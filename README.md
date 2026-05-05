# PeatClim v1.0

## Overview

This repository contains the code and configuration files used to build and evaluate **PeatClim v1.0**, a climate-driven machine-learning model for global peatland distribution designed for climate-model applications, especially in paleoclimate settings.

The repository also includes a selected derived product: the modern HTLTPeatland projection raster (`data_processed/projections/HTLTPeatland_EM01_GBM012_WC5min_EM01_GBM013_WC5min.tif`; see the relevant `jobid` entries in the config files for the corresponding model settings). This file is archived for reproducibility and convenience, and can also be regenerated from the provided code and required external inputs.

The workflow:
1. prepares modern peatland occurrence data from PEATMAP and Peat-DBase,
2. defines high- and low-temperature peatland groups (HTPeatland, LTPeatland),
3. prepares climate predictor layers,
4. trains and projects single and ensemble SDMs with R package `biomod2`,
5. integrates HT and LT predictions into a unified **HTLTPeatland** probability,
6. evaluates model performance and propagated uncertainty.

## Citation

If you use this repository, please cite:

Chen, L., Valdes, P., & Farnsworth, A. (2026). PeatClim: A climate-driven machine-learning model for predicting potential paleo-peatland distribution and its key climate controls (v1.0). Zenodo. https://doi.org/10.5281/zenodo.20039650

## Repository structure

```text
PeatClim_v1.0/
├── code/
├── config/
├── data_external/
├── data_processed/
├── outputs/
├── hpc/
├── README.md
└── sessionInfo_project.txt
```

## Main scripts

### Occurrence preparation

- `01a_prepare_PEATMAP_occurrence.R`
- `01b_prepare_Peat-DBase_occurrence.R`
- `01c_integrate_PEATMAP_Peat-DBase.R`
- `01d_define_plot_HT_LT_groups.R`

### Predictor preparation

PeatClim v1.0 uses bioclimatic variables from WorldClim v2.1 at 5 arc-min resolution as training predictors and as an example projection scenario. Other climate datasets, including WorldClim v2.1 aggregated to N48 resolution and simulated preindustrial (PI) climate data from the HadCM3BL model, are used only for projection tests.

- `02a_upscale_WorldClim_to_N48.R`
- `02b_prepare_PI_19_biovars.R`
- `03_calculate_bioclim_correlations.R`

### Model building, evaluation and projection

- `04_build_single_models.R`
- `05_project_single_models.R`
- `06_plot_TSS_across_algos_and_PAs.R`
- `07_plot_TSS_across_predictor_sets.R`
- `08_build_ensemble_models.R`
- `09_project_ensemble_models.R`

### HT/LT integration and uncertainty

- `10a_plot_schematic_weighted_probability.R`
- `10b_integrate_and_project_HTLTPeatland.R`
- `10c_evaluate_HTLTPeatland.R`
- `11_calculate_uncertainty_HTLTPeatland.R`

## External data

This archive does not redistribute large third-party raw datasets. These should be downloaded from their original repositories and placed in `data_external/` as described below.

### PEATMAP

Place the PEATMAP shapefiles under:

```
data_external/PEATMAP/
```

Dataset citation:  
Xu, J., Morris, P. J., Liu, J., and Holden, J. (2017). _PEATMAP: Refining estimates of global peatland distribution based on a meta-analysis._ University of Leeds. [Dataset]. [https://doi.org/10.5518/252](https://doi.org/10.5518/252)

Associated paper:  
Xu, J., Morris, P. J., Liu, J., and Holden, J. (2018). PEATMAP: Refining estimates of global peatland distribution based on a meta-analysis. _CATENA_, 160, 134–140.

### Peat-DBase

Place the Peat-DBase file under:

```
data_external/Peat-DBase/Peat_DBase_version_1_0_0_b.csv
```

Dataset citation:  
Skye, J. et al. (2025). _Peat-DBase v.1: A Compiled Database of Global Peat Depth Measurements (1.0.0b)._ Zenodo. [https://doi.org/10.5281/zenodo.15530644](https://doi.org/10.5281/zenodo.15530644)

Associated paper:  
Skye, J. et al. (2025). Peat-DBase v.1: a compiled database of global peat depth measurements. _Earth System Science Data_, 17, 7313–7330.

### WorldClim

Place the WorldClim v2.1 bioclimatic rasters under:

```
data_external/WorldClim/wc2.1_5m_bio/
```

Citation:  
Fick, S. E. and Hijmans, R. J. (2017). _WorldClim 2: new 1 km spatial resolution climate surfaces for global land areas_. _International Journal of Climatology_, 37(12), 4302–4315.

### Preindustrial climate inputs

Target land–sea masks and preindustrial climate inputs from HadCM3BL should be placed under:

```
data_external/HadCM3BL/
```

Information on the `tdezc` simulation can be found at https://www.paleo.bristol.ac.uk/ummodel/data/tdezc/standard_new_html/tdezc.html. 

General information on the BRIDGE simulation archive is available at: https://www.paleo.bristol.ac.uk/resources/simulations/.

For guidance on navigating the archive and accessing data, see: [Using_BRIDGE_webpages.pdf](https://www.paleo.bristol.ac.uk/ummodel/scripts/papers/Using_BRIDGE_webpages.pdf), especially Section 5, “Accessing this Data”.

Citation: 
Valdes, P. J. et al. (2017). _The BRIDGE HadCM3 family of climate models: HadCM3@Bristol v1.0_, _Geosci. Model Dev._, 10, 3715–3743, https://doi.org/10.5194/gmd-10-3715-2017, 2017.

## Running the workflow

Scripts use project-relative paths and should be run from the project root.

The model training workflow was executed on an HPC cluster using SLURM. An example launcher script is provided in `hpc/submit_r_job.sh`. Users on other systems may need to adapt scheduler directives, module-loading commands, resource requests, and account settings.

Example submission pattern:

```
sbatch hpc/submit_r_job.sh -r code/04_build_single_models.R -c config/config_build_single_models.csv -j <jobid>
```

## Software

This project uses R. R and packages details are provided in `sessionInfo_project.txt`.
