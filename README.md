# Multi-modal comparison of primary and stem cell-derived β-cells nominates targets for maturation

**Maghera J, Ellis CE, Spigelman AF, Smith N, Sasaki S, Lynn FC, MacDonald PE**

---

## Overview

This repository contains the analytical code for a multi-omic, single-cell comparison of stem cell-derived β-cells (SCβ-cells) and primary human β-cells, integrating patch-clamp electrophysiology, scRNA-seq, and SCENIC regulatory network inference. The analysis identifies transcriptional-electrophysiological coupling differences between SCβ-cells and primary β-cells and nominates SREBP1 as a candidate maturation target.

---

## Repository structure

```
bash/
  01_pyscenic_grn_array.sh          # SLURM array: 100 pySCENIC GRN runs
  02_pyscenic_ctx_aucell.sh         # pySCENIC ctx + aucell on consensus GRN

python/
  01_h5ad_to_loom.py                # Convert Seurat h5ad to loom for pySCENIC
  02_aggregate_consensus_grn.py     # Aggregate 100 GRN runs → consensus GRN
                                    # (edges present in ≥90/100 runs retained)
  03_auc_extraction.py              # Extract AUC matrix from pySCENIC loom output
  04_control_energy_bootstrap.py    # Forced response simulation + bootstrap
                                    # ranking of TF control energy

R/
  01_data_ingestion_JM.R            # STARsolo counts → Seurat object,
                                    # metadata merging, decontX, Scrublet (JM study)
  02_integrate_patchseq_datasets.R  # Merge JM + external patch-seq studies,
                                    # SCTransform, cells.plus object
  03_imputation.R                   # MICE multiple imputation of missing
                                    # electrophysiology values (m=50)
  04_differential_expression.R      # Bootstrapped pseudobulk DE (glmGamPoi),
                                    # LOSO stability, triage tiers
  05_correlations.R                 # Bootstrapped Spearman correlations,
                                    # DerSimonian-Laird meta-analysis, Δz,
                                    # GSEA on Δz-ranked gene lists
  06_SCENIC_networks.R              # Co-activity networks, bootstrapped PageRank,
                                    # LASSO projection of TFs onto MOFA factors
  07_DIABLO.R                       # Pseudobulk prep, block.splsda (mixOmics),
                                    # consensus features, LME4 annotations,
                                    # cross-block correlation heatmap
  08_MOFA2.R                        # Single-cell MOFA2 across 50 imputations,
                                    # consensus separation factor, UMAP,
                                    # lollipop plot of factor loadings

data/
  patch_key.xlsx                    # Electrophysiology parameter pretty-name lookup
```

---

## Data availability

Raw sequencing data for the JM_patchSeq study are deposited at [GEO accession — to be added].

Previously published patch-seq datasets used in this analysis:

| Study | Accession |
|-------|-----------|
| Camunas-Soler et al. | GSE124742 |
| Dai et al. | GSE164875 |
| Dos Santos et al. | GSE270484 |
| HPAP | https://hpap.pmacs.upenn.edu/ |
| CryoPatchSeq | Human Cell Atlas: CryoPancreaticIsletCellPatchSeq |
| Krentz et al. | GSE120522 |

pySCENIC auxiliary files (TF list, rankings databases, motif annotations) are available from the [pySCENIC resources page](https://scenic.aertslab.org/).

---

## Dependencies

### R packages
```r
# Core
tidyverse, Seurat, EnsDb.Hsapiens.v86, ensembldb

# Quality control
celda        # decontX
# Scrublet is run via reticulate (Python)

# Imputation
mice, AER, truncnorm, survival

# Differential expression
glmGamPoi, lme4, lmerTest, broom

# Correlations
matrixStats, future, future.apply, fgsea, msigdbr

# Multi-omics integration
mixOmics    # DIABLO
MOFA2       # unsupervised integration
uwot        # UMAP on consensus factors

# Network analysis
igraph, ggraph, tidygraph

# Visualization
ggplot2, patchwork, ComplexHeatmap, circlize, scico,
pheatmap, ggtext, ggrepel, tidytext
```

### Python packages
```
# pySCENIC pipeline
pyscenic, loompy, scanpy, anndata, sceasy

# GRN aggregation
pandas, dask

# Control energy
numpy, networkx, control, sklearn
```

### Compute environment
pySCENIC was run using Apptainer (Singularity) with the official pySCENIC container (`pyscenic.sif`) on a SLURM cluster. GRN inference (100 runs) was parallelized using SLURM job arrays.

---

## Reproducing the analysis

Scripts are numbered to reflect execution order. Each script loads the saved output of the previous script rather than re-running upstream steps. Key intermediate objects are saved as `.rds` or `.RData` files.

**Recommended order:**

1. `R/01` → `R/02` → `R/03` (data preparation and imputation)
2. `bash/01` → `python/02` → `bash/02` → `python/03` (pySCENIC pipeline)
3. `R/04` → `R/05` (differential expression and correlations)
4. `R/06` → `python/04` (SCENIC networks and control energy)
5. `R/07` → `R/08` (DIABLO and MOFA2)

---

## Notes on multiple imputation

Electrophysiology data contains substantial missingness. Missing values were handled using MICE (`R/03`), producing 50 completed datasets (`imp_bundle.rds`). Downstream analyses either randomly sample one imputed dataset per bootstrap iteration (correlations, MOFA2) or run independently on all 50 datasets and extract consensus (DIABLO). See Methods for full details.

---

## AI use statement

Claude (Anthropic, claude.ai) was used to refactor analytical code into reusable functions and loops, to annotate and document R and Python scripts for reproducibility, and to assist in drafting the methods section based on author-provided code and analytical decisions. All scientific decisions, interpretations, and conclusions are the authors' own.

---

## Citation

Maghera J, Ellis CE,  Spigelman AF, Smith N, Sasaki S, Lynn FC, MacDonald PE. (2026) Multi-modal comparison of primary and stem cell-derived β-cells nominates targets for maturation. bioRxiv, doi.org/10.64898/2026.06.04.730032.

---

## Contact

Cara Ellis: cara.e.ellis@gmail.com
Jasmine Maghera: maghera@ualberta.ca
Patrick MacDonald: pem@ualberta.ca
