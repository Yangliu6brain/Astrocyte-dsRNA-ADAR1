# Astrocyte-specific dsRNA–ADAR1 signalling links innate immunity to complement activation in Parkinson's disease

This repository holds the R code for the snRNA-seq analyses and figures of the manuscript. It covers donor covariates, the eight-module scoring framework, cell-type coupling analyses and robustness checks.

- **Preprint:** doi: [10.64898/2026.10.06.757209](https://doi.org/10.64898/2026.10.06.757209)
- **Archived code:** Zenodo, doi: [10.5281/zenodo.21315999](https://doi.org/10.5281/zenodo.21315999)
- **Source data:** the Zenodo data deposit (`Bioinformatics data_v2.xlsx`, `Figure_Source_Data_v2.xlsx`; see its `README_v2.md`)

Funded by the Aligning Science Across Parkinson's [000497, 020505, 027228] through the Michael J. Fox Foundation for Parkinson's Research and the National Health and Medical Research Council of Australia to GMH [Investigator Grant 2034292].

------------------------------------------------------------------------

## Version 2 (2026-11): what changed

- **Covariate correction.** In version 1, donor age, sex and PMI for the GSE178146 and GSE202210 donors were assigned by position, so these donors received another donor's values. `Full_merge_RDS_script_with_covariates.R` now reads them from `data/donor_covariates.csv`, one row per donor keyed by `sample_id`. It validates the table and stops on any mismatch.
- **Reruns.** All covariate-adjusted analyses were rerun with the corrected values. The scripts here are the copies that produced the corrected results.
- **Drawing scripts.** The final drawing scripts for Figures 5, 6 and 7 have been added.

------------------------------------------------------------------------

## 1. Pipeline

```
GEO raw data -> per-dataset QC (Step 1) -> merge + integrate + annotate + covariates
                                                   |
                                                   v
                                data/CX_Astrocytes_withCov.rds, data/CX_Neurons_Total_withCov.rds
                                                   |
                                                   v
                       Script5_TE_Control_Module_Analysis.R  (8 module scores, DESeq2, Fig 6)
                                -> results/CX_TE_Control/CX_8Module_AllResults.rds
                                                   |
     +-------------------+------------------+------+-----------+------------------+-------------------+
     v                   v                  v                  v                  v                   v
  Script4 (Fig 5)   Script_A1A2 (Fig 7)  Script2 fGSEA     Figure4B source    Braak (Supp Fig 3)  Script5 immune panels,
                                             (Fig 8)      (Supp Fig 2)                           S1 sensitivity
                                                                                
```

The eight modules are ADAR1, ADAR2, SNCA, dsRNA_Processing, dsRNA_EV_Transport, dsRNA_Sensing, Innate_Immune and TE_Control.

------------------------------------------------------------------------

## 2. Repository structure

```
.
├── README.md
├── LICENSE
├── Astrocyte-dsRNA-ADAR1/                      # unchanged from the previous release
│   ├── GSE_202210_Step1_Seurat_UPDATED.R, GSE_178146_Step1_Seurat_UPDATED.R, GSE253462_Step1_Seurat_UPDATED.R   # per-dataset QC
│   └── Script1_Combined_8Module_DESeq2_Analysis_V1_Final.R
├── Full_merge_RDS_script_with_covariates.R     # merge, integration, annotation, covariates
├── data/
│   └── donor_covariates.csv                    # age, sex, PMI per donor (corrected)
├── analysis/                                   # statistics; each writes a timestamped folder under results/
│   ├── Script5_TE_Control_Module_Analysis.R
│   ├── Script4_Module_x_Immune_ADD_CI_value_PMI_SENSITIVITY_FIXED.R
│   ├── Script_A1A2_PvsCol_Complement_ADD CIvalue.R
│   ├── Script2_KEGG_v4_FIXED.R
│   ├── Figure4B_SourceData_CSV.R
│   ├── MS2_Figure9B_Braak_Analysis-2.R
│   ├── Script5_PD_vs_Control_Immune_Panels.R
│   ├── Script_S1_Module_Sensitivity_Analysis.R
│   ├── NewScript4_Module_x_Immune_Correlation-2.R   # stratified Fig 5d rho quoted in the Results text (see §6)
│   
└── figures/                                    # drawing scripts
    ├── Fig4A/   Figure4_SourceData_CSV.R, Figure4A_SourceData_v2_block.R, Fig4A_aligned_v2_8mod_shapeDx.R
    ├── Fig5/    F3_figure5_heatmap_v3.R, F10_figure5_assembled_v7.R
    ├── Fig6/    Script6_Figure6_TE_Control_CorrelationMatrix.R, F9_figure6_v5.R
    ├── Fig7/    F13_figure7_v6.R
    ├── SuppFig2/       Fig4B_aligned_v1_layout2.R
```

The per-dataset QC scripts (`GSE_202210_Step1_Seurat_UPDATED.R`, `GSE_178146_Step1_Seurat_UPDATED.R`, `GSE253462_Step1_Seurat_UPDATED.R`) and Script 1 (`Script1_Combined_8Module_DESeq2_Analysis_V1_Final.R`) are unchanged from the previous release and stay in the folder `Astrocyte-dsRNA-ADAR1/`; they do not use covariates.

`data/` (apart from `donor_covariates.csv`), `results/` and `figures/` outputs are not redistributed. Raw matrices come from GEO (§4), and every processed object is regenerated by the scripts.

------------------------------------------------------------------------

## 3. Software environment

- **Single-cell analysis:** Seurat 5, harmony and scDblFinder.
- **Statistics:** DESeq2, fgsea, msigdbr and ppcor.
- **Data handling:** dplyr, tidyr, tibble and readxl.
- **Output and figures:** openxlsx, ggplot2 4.0.2, patchwork 1.3.2, ragg, RColorBrewer and scales.


------------------------------------------------------------------------

## 4. Data

The analyses re-use three public snRNA-seq datasets of human post-mortem cortex, 43 donors in total:

| Accession | Region | Donors (Control / PD) | Notes |
|---|---|---|---|
| GSE178146 | Cingulate cortex | 14 (7 / 7) | |
| GSE253462 | Cingulate cortex | 17 (8 / 9) | Age given as ranges, so the midpoint is used. PMI not available, so it is mean-imputed. |
| GSE202210 | Prefrontal cortex | 12 (6 / 6) | |

Download the raw matrices from GEO into `data/raw/<ACCESSION>/`. The donor covariates are in `data/donor_covariates.csv`, taken from the GEO records and original publications.

------------------------------------------------------------------------

## 5. How to run

Run everything from the repository root. Analysis scripts write timestamped folders such as `results/CX_TE_Control_20261113_1032/`. **Rename each folder to its stable name**, shown in the table below, before running the scripts that read it.

| Step | Script | Reads | Output folder → rename to |
|---|---|---|---|
| 1 | `Astrocyte-dsRNA-ADAR1/GSE*_Step1_Seurat_UPDATED.R` (unchanged) | `data/raw/` | `data/GSE*_merged_step1.rds` |
| 2 | `Full_merge_RDS_script_with_covariates.R` | Step 1 objects, `data/donor_covariates.csv` | `data/CX_Astrocytes_withCov.rds`, `data/CX_Neurons_Total_withCov.rds` |
| 3 | `analysis/Script5_TE_Control_Module_Analysis.R` | `data/` | `CX_TE_Control_<date>` → `results/CX_TE_Control` |
| 4 | `analysis/Script4_Module_x_Immune_ADD_CI_value_PMI_SENSITIVITY_FIXED.R` | `data/`, `results/CX_TE_Control` | `CX_Module_x_Immune_Correlation_<date>` → `results/CX_Module_x_Immune_Correlation` |
| 4 | `analysis/NewScript4_Module_x_Immune_Correlation-2.R` (stratified ρ quoted in the Fig 5 Results text) | `data/`, `results/CX_TE_Control` | also named `CX_Module_x_Immune_Correlation_<date>`: keep it apart from Script4's folder, e.g. → `results/CX_Module_x_Immune_Correlation_stratified` (table: `FigD_Stratified_Stats.xlsx`) |
| 4 | `analysis/Script_A1A2_PvsCol_Complement_ADD CIvalue.R` | `data/` | `CX_A1A2_PD_vs_Control_<date>` → `results/CX_A1A2_PD_vs_Control` |
| 4 | `analysis/Script2_KEGG_v4_FIXED.R` | `data/`, `results/CX_TE_Control` | `CX_Unified_fGSEA_<date>` → `results/CX_Unified_fGSEA` |
| 4 | `analysis/Figure4B_SourceData_CSV.R` | `data/`, `results/CX_TE_Control` | `CX_Fig4B_SourceData_<date>` → `results/CX_Fig4B_SourceData` |
| 4 | `analysis/MS2_Figure9B_Braak_Analysis-2.R` | `data/` | `results/MS2_Figure9B_Braak_<date>` |
| 4 | `analysis/Script5_PD_vs_Control_Immune_Panels.R` | `data/` | `results/CX_PD_vs_Ctrl_ImmunePanels_<date>` |
| 4 | `analysis/Script_S1_Module_Sensitivity_Analysis.R` | `data/`, Script 1 output `results/CX_7Mod/CX_7Module_AllResults.rds` | `results/CX_Sensitivity_<date>` |

RB3 also compares against two runs made with the version 1 covariates: `results/baseline/CX_Unified_fGSEA_published`, `results/baseline/CX_Unified_fGSEA_v4_oldcov` and `results/baseline/ASAPextract`. These cannot be regenerated from the corrected covariate table. RB3 needs them only for its "before vs after" columns.

**Figures** (run after the analysis steps):

``` r
# Figure 5
source("figures/Fig5/F3_figure5_heatmap_v3.R")        # -> figures/Figure5_heatmap/ (read by F10)
source("figures/Fig5/F10_figure5_assembled_v7.R")     # -> figures/Figure5/
# Figure 6
source("figures/Fig6/Script6_Figure6_TE_Control_CorrelationMatrix.R")  # -> results/CX_Figure6_TE_Control_<date>; rename to results/CX_Figure6_TE_Control
source("figures/Fig6/F9_figure6_v5.R")                # -> figures/Figure6/
# Figure 7
source("figures/Fig7/F13_figure7_v6.R")               # -> figures/Figure7/
# Supplementary Figure 2
source("figures/SuppFig2/Fig4B_aligned_v1_layout2.R") # -> figures/SuppFig2/

```

F10, F9 and F13 check the printed values against the analysis output (`stopifnot`). They also write a `*_printed_values.csv` next to each figure, listing every number printed in it.

------------------------------------------------------------------------

## 6. Script-to-manuscript mapping

| Manuscript item | Analysis script | Drawing script(s) |
|---|---|---|
| Figure 4A | — | `Figure4_SourceData_CSV.R` (+ `Figure4A_SourceData_v2_block.R`), `Fig4A_aligned_v2_8mod_shapeDx.R` |
| Figure 5a–d | `Script4_Module_x_Immune_ADD_CI_value_PMI_SENSITIVITY_FIXED.R` |
| Figure 6a, b | `Script5_TE_Control_Module_Analysis.R` | `Script6_Figure6_TE_Control_CorrelationMatrix.R` |
| Figure 7a–c | `Script_A1A2_PvsCol_Complement_ADD CIvalue.R` |
| Supplementary Figure 2 | `Figure4B_SourceData_CSV.R` | `Fig4B_aligned_v1_layout2.R` |
| Supplementary Figure 3 (Braak) | `MS2_Figure9B_Braak_Analysis-2.R` | — |
| Deposit: DE_results, correlation_Spearman, correlation_fisherZ | `Script5_TE_Control_Module_Analysis.R` | — |
| Deposit: Immune_fisherZ | `Script4_…_FIXED.R` | — |
| Deposit: Immune_Condition | `Script5_PD_vs_Control_Immune_Panels.R` | — |
| Deposit: Sensitivity_analysis | `Script_S1_Module_Sensitivity_Analysis.R` | — |
| Deposit: braak_analysis | `MS2_Figure9B_Braak_Analysis-2.R` | — |
| Deposit: fgsea_results| `Script2_KEGG_v4_FIXED.R` | — |


------------------------------------------------------------------------

## 7. License

Code is released under the **MIT License** (see `LICENSE`), as required by the ASAP Open Science Policy.

------------------------------------------------------------------------

## 8. Citation

- Preprint: doi: [10.64898/2026.10.06.757209](https://doi.org/10.64898/2026.10.06.757209)
- Archived code: doi: [10.5281/zenodo.21315999](https://doi.org/10.5281/zenodo.21315999)

------------------------------------------------------------------------

## 9. ASAP compliance

This repository is one of the research outputs of the ASAP-funded study above. It is released with an open reuse license and archived on Zenodo with a citable DOI. Processed data, wet-lab protocols (protocols.io) and microscopy images (BioImage Archive, [S-BIAD3762](https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD3762)) are deposited separately; see the manuscript's Key Resources Table.

**Correspondence:** Yang Liu, Yang.liu6@sydney.edu.au
