# ==============================================================================
# 06_SCENIC_networks_control_energy.R
# ==============================================================================
# Bootstrapped TF co-activity networks, PageRank influence analysis, GRN
# subnetwork visualization, and LASSO projection of TF AUC scores onto MOFA2
# consensus factors (B matrix for control energy estimation).
#
# Corresponds to Methods: "SCENIC regulon inference and co-activity networks"
# and "Network control energy and transcription factor prioritization"
#   "Bootstrapped co-activity networks were constructed by resampling cells
#    with replacement and computing pairwise Spearman correlations between
#    regulon AUC values... Edges retained in ≥50/100 bootstraps at ρ ≥ 0.55
#    were used to define consensus networks... PageRank influence scores were
#    computed on directed networks (ρ ≥ 0.30) and differences in mean
#    bootstrapped PageRank (Δ influence = SCβ − β) were used to rank TFs by
#    differential network centrality."
#
# Input:
#   20260113_cells_plus_filtered_mi_SCENIC.rds
#     - cells.plus: integrated Seurat object with AUC assay
#   cells_plus_mi_consensus_grn_results.csv
#     - Consensus GRN from python/02_aggregate_consensus_grn.py
#   bootstrapped_DE_triage_table.csv
#     - DE triage results from scripts/04_differential_expression.R
#   MOFA2_consensus_factors_with_metadata.csv
#     - Consensus MOFA2 factor scores, from scripts/08_MOFA2.R
#     - Note: run scripts/08_MOFA2.R first to generate this file
#
# Output:
#   20260115_SCB_beta_bs.RData
#     - SCB_beta_bs: bootstrapped PageRank influence results (summary + all boots)
#   20260115_bootstrapped_page_rank_influence_scores_revised.csv
#     - summary_df: mean Δ influence and direction per TF
#   20260126_coactivity_networks_with_direction_revised.svg
#     - Co-activity network panels for SCβ, S6D1, and β-cells
#   20260326_tf_networks_with_direction_revised_NS.svg
#     - GRN subnetwork: top 13 TFs + PDX1 with targets colored by DE
#   TF_to_MOFA_lasso_matrix_consensus.csv
#     - B matrix: TF × MOFA factor LASSO coefficients
#     - Input to python/04_control_energy_bootstrap.py
#
# Dependencies:
#   tidyverse, Seurat, igraph, ggraph, ggrepel, patchwork, scico,
#   glmnet, EnsDb.Hsapiens.v86, ensembldb, clusterProfiler, org.Hs.eg.db, rrvgo
#
# Note on script ordering:
#   This script requires MOFA2_consensus_factors_with_metadata.csv from
#   scripts/08_MOFA2.R (Section 6, LASSO). If running in strict order, run
#   scripts/07 and 08 first, then return here for the LASSO section.
#   Sections 1-5 (networks) can be run independently of MOFA2.
# ==============================================================================

library("tidyverse")
library("Seurat")
library("igraph")
library("ggraph")
library("ggrepel")
library("patchwork")
library("scico")
library("glmnet")
library("EnsDb.Hsapiens.v86")
library("ensembldb")
library("clusterProfiler")
library("org.Hs.eg.db")
library("rrvgo")

# ── Colour palette ────────────────────────────────────────────────────────────
cluster_cols2 <- c(
  "#ff7db1","#e2ff50","#713be8","#669b00","#c734e9","#01a457","#e400b3",
  "#8cffb0","#010f92","#ffd26e","#015fde","#005a02","#ff6de1","#01deaf",
  "#9d0098","#4f6600","#a282ff","#867100","#170047","#b06100","#01b4e8",
  "#ff5662","#02cfcb","#ff3a85","#006c55","#ac004a","#baecff","#862200",
  "#006d95","#ff7d56","#003767","#ffcaa5","#33002d","#ffd6db","#101c00",
  "#f6bbff","#001e22","#ff9da8","#360c00","#5f3800","#001024","#d0f700"
)

# ═══════════════════════════════════════════════════════════════════════════════
# Section 1: Load data and prepare AUC matrix
# ═══════════════════════════════════════════════════════════════════════════════

cells.plus <- readRDS("20260113_cells_plus_filtered_mi_SCENIC.rds")

# AUC matrix: regulons (rows) × cells (columns)
auc_matrix <- cells.plus[["AUC"]]$counts

metadata <- cells.plus@meta.data
metadata$CellID <- rownames(metadata)

# Cell groups for co-activity network analysis
cells_SCB  <- metadata$CellID[metadata$celltype == "SCβ"]
cells_beta <- metadata$CellID[metadata$celltype == "beta"]
cells_S6D1 <- metadata$CellID[metadata$celltype == "S6D1"]

# Regulon network (TF → target gene list) from pySCENIC loom
# Used to restrict PageRank analysis to TF-only nodes
library("SCopeLoomR")
loom <- open_loom("20260115_cells_plus_mi_pyscenic_final.loom")
regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
regulons          <- regulonsToGeneLists(regulons_incidMat)
regulon_network   <- purrr::map(names(regulons), function(tf) {
  tibble(TF = tf, Target = regulons[[tf]])
}) %>% bind_rows()
close_loom(loom)

# ═══════════════════════════════════════════════════════════════════════════════
# Section 2: Bootstrapped co-activity networks
# ═══════════════════════════════════════════════════════════════════════════════
# For each cell group, build 100 bootstrap networks by resampling cells with
# replacement and computing pairwise Spearman correlations between regulon AUCs.
# Edges are retained at ρ > 0.55 (elbow of edge-density curve).
# Consensus edges: present in ≥50/100 bootstrap runs.

# ── Bootstrap co-activity network construction ────────────────────────────────
bootstrap_auc_networks <- function(auc_matrix, group_cells,
                                   n_boot = 100, cor_thresh = 0.55) {
  net_list <- vector("list", n_boot)
  for (i in seq_len(n_boot)) {
    sample_cells <- sample(group_cells, replace = TRUE)
    auc_sample   <- auc_matrix[, sample_cells]
    cor_mat      <- cor(t(auc_sample), method = "spearman")
    cor_long     <- as.data.frame(as.table(cor_mat)) %>%
      rename(from = Var1, to = Var2, correlation = Freq) %>%
      filter(from != to, abs(correlation) > cor_thresh)
    net_list[[i]] <- igraph::graph_from_data_frame(cor_long, directed = FALSE)
  }
  net_list
}

# ── Network property quantification ──────────────────────────────────────────
quantify_boot_networks <- function(net_list) {
  net_stats <- purrr::map(net_list, ~ data.frame(
    num_nodes  = vcount(.x),
    num_edges  = ecount(.x),
    avg_degree = mean(degree(.x)),
    density    = edge_density(.x),
    modularity = modularity(cluster_louvain(.x)),
    stringsAsFactors = FALSE
  )) %>% bind_rows()
  net_stats %>% summarise(across(everything(), list(mean = mean, sd = sd)))
}

# ── Edge frequency counter ────────────────────────────────────────────────────
# Counts how many bootstraps each edge appeared in.
# Pre-filters to n > 10 for memory efficiency (main threshold applied later).
get_edge_counts <- function(net_list) {
  edge_dfs <- purrr::map(net_list, ~ as_data_frame(.x, what = "edges") %>%
    mutate(edge = paste(pmin(from, to), pmax(from, to), sep = "_")))
  bind_rows(edge_dfs) %>%
    dplyr::count(edge) %>%
    separate(edge, into = c("from", "to"), sep = "_") %>%
    dplyr::filter(n > 10)
}

cat("Building co-activity bootstrap networks...\n")
net_list_SCB  <- bootstrap_auc_networks(auc_matrix, cells_SCB)
net_list_beta <- bootstrap_auc_networks(auc_matrix, cells_beta)
net_list_S6D1 <- bootstrap_auc_networks(auc_matrix, cells_S6D1)

# Network property summaries
boot_nets_SCB  <- quantify_boot_networks(net_list_SCB)
boot_nets_beta <- quantify_boot_networks(net_list_beta)
boot_nets_S6D1 <- quantify_boot_networks(net_list_S6D1)

cat("Network property summaries:\n")
print(bind_rows(SCB = boot_nets_SCB, beta = boot_nets_beta,
                S6D1 = boot_nets_S6D1, .id = "group"))

# Consensus edges (≥50/100 bootstraps)
edge_SCB  <- get_edge_counts(net_list_SCB)
edge_beta <- get_edge_counts(net_list_beta)
edge_S6D1 <- get_edge_counts(net_list_S6D1)

# ═══════════════════════════════════════════════════════════════════════════════
# Section 3: Bootstrapped PageRank Δ influence (SCβ vs β)
# ═══════════════════════════════════════════════════════════════════════════════
# For each bootstrap iteration, compute PageRank centrality in directed
# co-activity graphs (ρ ≥ 0.30) built separately for SCβ and β-cells.
# Δ influence = PageRank_SCβ − PageRank_β, summarized across 100 iterations.
# Note: directed graphs are built from symmetric correlation matrices; the
# direction assignment is arbitrary but consistent within each bootstrap run.
# Node coloring in the final figure reflects mean_delta (positive = SCβ-dominant).

bootstrap_influence_deltas_regulon <- function(auc_matrix, meta_data,
                                               celltypes_to_compare,
                                               n_iter = 100,
                                               cor_thresh = 0.6) {
  tf_list              <- unique(regulon_network$TF)
  delta_influence_list <- vector("list", n_iter)

  for (i in seq_len(n_iter)) {
    cat("Iteration", i, "\n")

    # Subsample 80% of cells per group (without replacement for stable estimates)
    sampled_cells <- meta_data %>%
      dplyr::filter(celltype %in% celltypes_to_compare) %>%
      group_by(celltype) %>%
      sample_frac(0.8) %>%
      pull(CellID)

    auc_sub  <- auc_matrix[, sampled_cells, drop = FALSE]
    meta_sub <- meta_data %>% dplyr::filter(CellID %in% sampled_cells)

    # Within-group correlation matrices → directed graphs
    cor_list <- lapply(celltypes_to_compare, function(ct) {
      cells <- meta_sub$CellID[meta_sub$celltype == ct]
      cor(t(auc_sub[, cells, drop = FALSE]),
          method = "spearman", use = "pairwise.complete.obs")
    })
    names(cor_list) <- celltypes_to_compare

    g_list <- lapply(cor_list, function(cormat) {
      cor_edges <- as.data.frame(as.table(cormat)) %>%
        dplyr::filter(Var1 != Var2, abs(Freq) > cor_thresh) %>%
        rename(from = Var1, to = Var2, weight = Freq)
      graph_from_data_frame(cor_edges, directed = TRUE)
    })

    # PageRank with edge weights = |correlation|
    influence_scores <- lapply(g_list, function(g) {
      score <- page_rank(g, directed = TRUE,
                         weights = abs(E(g)$weight))$vector
      tibble(TF = names(score), influence = score)
    })

    delta_influence_list[[i]] <- full_join(
      influence_scores[[1]], influence_scores[[2]],
      by = "TF", suffix = c("_1", "_2")
    ) %>% mutate(delta = influence_1 - influence_2)
  }

  all_deltas <- bind_rows(delta_influence_list, .id = "bootstrap")
  summary_df <- all_deltas %>%
    group_by(TF) %>%
    summarize(
      mean_delta = mean(delta,          na.rm = TRUE),
      sd_delta   = sd(delta,            na.rm = TRUE),
      n_top      = sum(rank(-abs(delta)) <= 10, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    arrange(desc(abs(mean_delta)))

  list(summary = summary_df, all_boots = all_deltas)
}

cat("Running bootstrapped PageRank influence analysis...\n")
SCB_beta_bs <- bootstrap_influence_deltas_regulon(
  auc_matrix, metadata, c("SCβ", "beta"), cor_thresh = 0.3
)
save(SCB_beta_bs, file = "20260115_SCB_beta_bs.RData")

# Summarize: add direction and stability columns
summary_df <- SCB_beta_bs$summary %>%
  mutate(
    stability          = 1 / sd_delta,
    influence_direction = case_when(
      mean_delta > 0 ~ "SCB",
      mean_delta < 0 ~ "Beta",
      TRUE           ~ "Equal"
    )
  ) %>%
  arrange(desc(abs(mean_delta))) %>%
  mutate(TF = str_replace(TF, "\\-\\(", "\\("))

write_csv(summary_df,
          file = "./for_manuscript/20260115_bootstrapped_page_rank_influence_scores_revised.csv")

cat("Top 10 TFs by |Δ influence|:\n")
print(head(summary_df, 10))

# ═══════════════════════════════════════════════════════════════════════════════
# Section 4: Consensus co-activity network visualization
# ═══════════════════════════════════════════════════════════════════════════════
# Build consensus graphs from edges present in ≥50/100 bootstraps.
# Node color = bootstrapped PageRank Δ influence direction (from summary_df).

threshold <- 50  # consensus threshold: ≥50/100 bootstraps

build_consensus_graph <- function(edge_df, summary_df, threshold) {
  edges <- edge_df %>%
    dplyr::filter(n >= threshold) %>%
    mutate(from = str_replace(from, "\\-\\(", "\\("),
           to   = str_replace(to,   "\\-\\(", "\\("))

  g <- graph_from_data_frame(edges, directed = FALSE)
  V(g)$direction <- summary_df$influence_direction[
    match(V(g)$name, summary_df$TF)
  ]
  V(g)$delta <- summary_df$mean_delta[match(V(g)$name, summary_df$TF)]
  g
}

g_SCB_consensus  <- build_consensus_graph(edge_SCB,  summary_df, threshold)
g_S6D1_consensus <- build_consensus_graph(edge_S6D1, summary_df, threshold)
g_beta_consensus <- build_consensus_graph(edge_beta,  summary_df, threshold)

# Fix regulon name formatting for S6D1 and beta (already done for SCB above)
V(g_S6D1_consensus)$name <- str_replace(V(g_S6D1_consensus)$name, "\\-\\(", "\\(")
V(g_beta_consensus)$name  <- str_replace(V(g_beta_consensus)$name,  "\\-\\(", "\\(")

direction_colors <- c("SCB" = "#A6D278", "Beta" = "#053059", "None" = "#B2F2FD")

plot_network <- function(g, title) {
  set.seed(42)
  layout <- ggraph::create_layout(g, layout = "fr")
  ggraph(layout) +
    geom_edge_link(aes(alpha = n), color = "grey70") +
    geom_node_point(aes(fill = direction), shape = 21, size = 5, color = "black") +
    scale_fill_manual(values = direction_colors, name = "Upregulated") +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    ggtitle(title) +
    guides(fill = guide_legend("Upregulated",
                               override.aes = list(shape = 21))) +
    theme_void()
}

p1 <- plot_network(g_S6D1_consensus, "S6D1 Coactivity Network")
p2 <- plot_network(g_SCB_consensus,  "SCβ Coactivity Network")
p3 <- plot_network(g_beta_consensus, "Beta Coactivity Network")

ggsave("./for_manuscript/20260126_coactivity_networks_with_direction_revised.svg",
       plot = p1 + p2 + p3,
       width = 45, height = 10, units = "cm", dpi = 600)

# ═══════════════════════════════════════════════════════════════════════════════
# Section 5: GRN subnetwork visualization (Fig 5)
# ═══════════════════════════════════════════════════════════════════════════════
# Build a directed subnetwork from the consensus GRN, restricted to the top 13
# TFs by |Δ influence| plus PDX1. Edges: importance > 95th percentile.
# Node color = Δ influence direction (TF nodes) or DE direction (target nodes).

adj_df <- read_csv("cells_plus_mi_consensus_grn_results.csv") %>%
  mutate(TF = paste0(TF, "(+)"))

importance_threshold <- quantile(adj_df$importance, 0.95)

top_tfs <- c(
  summary_df %>%
    top_n(13, abs(mean_delta)) %>%
    pull(TF) %>%
    str_replace("\\-\\(", "\\_\\("),
  "PDX1(+)"
)

reg_edges_filtered <- adj_df %>%
  dplyr::filter(TF %in% top_tfs & importance > importance_threshold)

tf_network <- graph_from_data_frame(reg_edges_filtered, directed = TRUE)

# Load DE triage results for target gene coloring
triage2 <- read_csv("bootstrapped_DE_triage_table.csv")
deg_labeled <- triage2 %>%
  dplyr::filter(triage_tier2 %in% c("Tier1A_strong_stable_JM_supported",
                                     "Tier2_strong_sensitive"))

# Annotate nodes with influence direction (TFs) and DE direction (targets)
influence_summary <- summary_df %>%
  dplyr::select(TF, mean_delta) %>%
  mutate(TF_clean = str_remove(TF, "-\\(\\+\\)"))

tf_list <- unique(summary_df$TF)
V(tf_network)$node_type <- ifelse(V(tf_network)$name %in% tf_list, "TF", "target")
V(tf_network)$logFC     <- deg_labeled$mean_logFC[
  match(V(tf_network)$name, deg_labeled$gene_name)
]

node_df <- as_data_frame(tf_network, what = "vertices") %>%
  mutate(
    is_tf     = node_type == "TF",
    tf_symbol = ifelse(is_tf, str_remove(name, "_\\(\\+\\)"), NA),
    mean_delta = ifelse(
      is_tf,
      influence_summary$mean_delta[match(tf_symbol, influence_summary$TF_clean)],
      NA
    ),
    color_group = case_when(
      is_tf & mean_delta > 0    ~ "Higher in SCβ",
      is_tf & mean_delta < 0    ~ "Higher in Beta",
      !is_tf & logFC > 0.25     ~ "Higher in SCβ",
      !is_tf & logFC < -0.25    ~ "Higher in Beta",
      TRUE                      ~ "NS"
    ),
    node_class = ifelse(is_tf, "TF", "Target")
  )

# Align node attributes to graph ordering
node_df <- node_df[match(V(tf_network)$name, node_df$name), ]
V(tf_network)$color_group <- node_df$color_group
V(tf_network)$node_class  <- node_df$node_class
V(tf_network)$label_flag  <- node_df$node_type == "TF"
V(tf_network)$name        <- str_remove(V(tf_network)$name, "_")

combined_palette <- c(
  "Higher in SCβ"  = "#A6D278",
  "Higher in Beta" = "#053059",
  "NS"             = "#C0CFDE"
)

set.seed(123)
g_fig5B <- ggraph(tf_network, layout = "stress") +
  geom_edge_link(aes(alpha = importance), color = "gray80") +
  geom_node_point(
    data = function(x) dplyr::filter(x, color_group == "NS"),
    aes(fill = color_group, shape = node_class), color = "black", size = 3
  ) +
  geom_node_point(
    data = function(x) dplyr::filter(x, color_group == "Higher in Beta"),
    aes(fill = color_group, shape = node_class), color = "black", size = 4
  ) +
  geom_node_point(
    data = function(x) dplyr::filter(x, color_group == "Higher in SCβ"),
    aes(fill = color_group, shape = node_class), color = "black", size = 4.5
  ) +
  geom_node_text(aes(label = ifelse(label_flag, name, "")),
                 repel = TRUE, size = 5, max.overlaps = 500) +
  scale_shape_manual(values = c("TF" = 21, "Target" = 22), name = "Node Type") +
  scale_fill_manual(values = combined_palette, name = "Expression/Influence") +
  guides(fill  = guide_legend("Expression/Influence",
                              override.aes = list(shape = 21)),
         alpha = "none") +
  theme_void()

ggsave("./for_manuscript/20260326_tf_networks_with_direction_revised_NS.svg",
       plot = g_fig5B, bg = "white",
       width = 20, height = 12, units = "cm", dpi = 600)

# ── Theme-faceted subnetworks (Supplementary) ─────────────────────────────────
# GO enrichment on subnetwork target genes, reduced with rrvgo,
# used to color target nodes by functional theme in supplementary figure.

target_genes  <- unique(reg_edges_filtered$target)
entrez_ids    <- bitr(target_genes, fromType = "SYMBOL",
                      toType = "ENTREZID", OrgDb = org.Hs.eg.db)

ego <- enrichGO(
  gene          = entrez_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "none",
  qvalueCutoff  = 1,
  readable      = TRUE
)

# Reduce GO terms using semantic similarity
simMatrix    <- calculateSimMatrix(ego$ID, orgdb = "org.Hs.eg.db",
                                   ont = "BP", method = "Rel")
term_scores  <- setNames(rep(1, length(ego$ID)), ego$ID)
reducedTerms <- reduceSimMatrix(simMatrix, term_scores,
                                threshold = 0.7, orgdb = "org.Hs.eg.db")

# Assign biological themes to reduced GO parent terms
reducedTerms <- reducedTerms %>%
  mutate(theme = case_when(
    str_detect(parentTerm,
               regex("mitochondria|oxidoreductase|respiratory|NADH|pyruvate|coenzyme|energy homeostasis|lipid",
                     ignore_case = TRUE)) ~ "Mitochondrial/Metabolic",
    str_detect(parentTerm,
               regex("vesicle|Golgi|transport|endosomal|tethering|plasma membrane organization|receptor transport",
                     ignore_case = TRUE)) ~ "Vesicle/Transport",
    str_detect(parentTerm,
               regex("differentiation|development|maturation|stem cell|organ regeneration|gland|morphogenesis",
                     ignore_case = TRUE)) ~ "Differentiation/Maturation",
    str_detect(parentTerm,
               regex("transcription|translation|RNA polymerase|regulation of gene expression",
                     ignore_case = TRUE)) ~ "Gene Expression Regulation",
    str_detect(parentTerm,
               regex("neurotransmitter|neuron|spine|nerve|recognition|synapse",
                     ignore_case = TRUE)) ~ "Neural/Neuroendocrine",
    TRUE ~ "Other"
  ))

# Build gene → theme lookup
gene_theme_df <- ego@result %>%
  dplyr::select(Description, geneID) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  left_join(reducedTerms %>% dplyr::select(term, theme),
            by = c("Description" = "term")) %>%
  distinct(geneID, theme) %>%
  rename(gene = geneID)

# Update node_df with theme, using name without underscore
node_df$name <- str_remove(node_df$name, "_")
node_df <- node_df %>%
  left_join(gene_theme_df, by = c("name" = "gene")) %>%
  mutate(theme = ifelse(is_tf, "TF", theme),
         theme = ifelse(is.na(theme), "Unannotated", theme))
V(tf_network)$theme <- node_df$theme[match(V(tf_network)$name, node_df$name)]

# Theme-faceted plots (Supplementary)
theme_list <- c("Differentiation/Maturation", "Mitochondrial/Metabolic",
                "Vesicle/Transport", "Neural/Neuroendocrine")

# Canonical gene labels retained across theme facets
deg_label <- c("INS","PDX1","MAFA","COX6B1","COX5B","COX7A2","COX7C",
               "MTOR","IAPP","CHGA","PTPRN","SLC25A6")
tf_label  <- node_df %>% dplyr::filter(is_tf) %>% pull(name)

theme_plots <- lapply(theme_list, function(th) {
  set.seed(42)
  keep_genes <- node_df %>%
    dplyr::filter(theme == th | is_tf) %>%
    pull(name)

  g_sub <- induced_subgraph(tf_network,
                             vids = V(tf_network)[name %in% keep_genes])
  V(g_sub)$color_group <- node_df$color_group[match(V(g_sub)$name, node_df$name)]
  V(g_sub)$node_class  <- node_df$node_class[match(V(g_sub)$name, node_df$name)]
  V(g_sub)$label_flag  <- V(g_sub)$name %in% unique(c(deg_label, tf_label))

  ggraph(g_sub, layout = "stress") +
    geom_edge_link(aes(alpha = importance), color = "gray60") +
    geom_node_point(
      data = function(x) dplyr::filter(x, color_group == "NS"),
      aes(fill = color_group, shape = node_class), color = "black", size = 3
    ) +
    geom_node_point(
      data = function(x) dplyr::filter(x, color_group == "Higher in Beta"),
      aes(fill = color_group, shape = node_class), color = "black", size = 4
    ) +
    geom_node_point(
      data = function(x) dplyr::filter(x, color_group == "Higher in SCβ"),
      aes(fill = color_group, shape = node_class), color = "black", size = 4.5
    ) +
    geom_node_text(aes(label = ifelse(label_flag, name, "")),
                   repel = TRUE, size = 5, max.overlaps = 30) +
    scale_fill_manual(
      values = c("Higher in SCβ"  = "#A6D278",
                 "Higher in Beta" = "#053059",
                 "NS"             = "#C0CFDE"),
      name = "Upregulated"
    ) +
    scale_shape_manual(values = c("TF" = 21, "Target" = 22),
                       name = "Node Type") +
    scale_alpha_continuous(range = c(0.1, 0.6), guide = "none") +
    ggtitle(th) +
    guides(fill = guide_legend("Influence",
                               override.aes = list(shape = 21))) +
    theme_void() +
    theme(plot.title = element_text(size = 14, hjust = 0.5),
          legend.position = "right")
})

ggsave(
  "./for_manuscript/20260127_tf_networks_with_direction_theme_facet_revised.svg",
  plot = wrap_plots(theme_plots, ncol = 2) + plot_layout(guides = "collect"),
  width = 30, height = 20, units = "cm", dpi = 600
)

# ═══════════════════════════════════════════════════════════════════════════════
# Section 6: LASSO projection of TF AUC scores onto MOFA2 consensus factors
# ═══════════════════════════════════════════════════════════════════════════════
# Regresses each MOFA2 consensus factor score on TF AUC values using LASSO
# (alpha = 1, lambda selected by 10-fold CV). The resulting coefficient matrix
# (TF × factors) is the B matrix used in python/04_control_energy_bootstrap.py.
#
# Note: Z_consensus (cells × MOFA factors) must be loaded from
# MOFA2_consensus_factors_with_metadata.csv — run scripts/08_MOFA2.R first.

Z_consensus <- read_csv("MOFA2_consensus_factors_with_metadata.csv") %>%
  dplyr::select(CellID, starts_with("Factor")) %>%
  column_to_rownames("CellID") %>%
  as.matrix()

# Align AUC matrix to MOFA cells; strip "-" from regulon names for LASSO
auc_t <- t(auc_matrix)  # cells × TFs
auc_t <- auc_t[rownames(Z_consensus), , drop = FALSE]
colnames(auc_t) <- str_remove(colnames(auc_t), pattern = "-")

X            <- as.matrix(auc_t)
Y            <- Z_consensus
factor_names <- colnames(Y)

# LASSO CV per factor
lasso_results <- map(seq_len(ncol(Y)), function(i) {
  y       <- Y[, i]
  cv_fit  <- cv.glmnet(X, y, alpha = 1)
  coefs   <- as.matrix(coef(cv_fit, s = cv_fit$lambda.min))
  coefs   <- coefs[-1, , drop = FALSE]  # remove intercept
  tibble(TF = rownames(coefs), Coefficient = coefs[, 1], Factor = factor_names[i])
}) %>% bind_rows()

lasso_results_filtered <- lasso_results %>%
  dplyr::filter(Coefficient != 0) %>%
  mutate(TF_clean = gsub("\\(\\+\\)", "", TF))

lasso_wide <- lasso_results_filtered %>%
  dplyr::select(TF_clean, Coefficient, Factor) %>%
  pivot_wider(names_from = Factor, values_from = Coefficient, values_fill = 0) %>%
  rename(TF = TF_clean)

write.csv(lasso_wide, "TF_to_MOFA_lasso_matrix_consensus.csv", row.names = FALSE)

cat("B matrix saved: TF_to_MOFA_lasso_matrix_consensus.csv\n")
cat("TFs with non-zero coefficients:", nrow(lasso_wide), "\n")
cat("Next step: run python/04_control_energy_bootstrap.py\n")
