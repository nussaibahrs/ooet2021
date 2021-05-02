# The evolution of the latitudinal diversity gradient of Cenozoic marine plankton
*Nussaïbah B. Raja & Wolfgang Kiessling*

[![](https://img.shields.io/badge/doi-10.17605/OSF.IO/NGK54-orange.svg)](https://doi.org/10.17605/OSF.IO/NGK54)
[![](https://img.shields.io/github/languages/code-size/nussaibahrs/ooet_plankton.svg)](https://github.com/nussaibahrs/ooet_plankton)

  - [Description](#description)
  - [Requirements](#requirements)
  - [Setup](#setup)
  - [Scripts](#scripts)
  - [Troubleshooting](#troubleshooting)

## Description

This repository contains all **R scripts** (in `/scripts`) necessary to evaluate the diversification dynamics and selectivity in, and the dispersal of marine plankton over the last 66 million years. The **original data files** can be downloaded at http://doi.org/10.17605/OSF.IO/NGK54. They were not included in the repository due to file size limits. Please download the data files and replace the contents of the data folder. The outputs of the scripts are provided in the `/output` and figures in the `/figs` folder. Please cite the study as:

Raja NB, Kiessling W. 2021 Out of the extratropics: the evolution of the latitudinal diversity gradient of Cenozoic marine plankton. Proc. R. Soc. B 20210545. https://doi.org/10.1098/rspb.2021.0545 

## Requirements

This code was developed in `R 4.0.0`. It is therefore recommended to use the same or any more up-to-date version of R for reproducing the analyses in this study.

## Setup
You will need to either use the Rstudio project environment or set your working directory to the root of this folder. 

To install all required depdendencies (packages), run:

```
source(file.path("inst","dependencies.R"))
```

## Scripts
The `scripts/` folder contains all the code generated for the above mentioned study. The folder contains **5** R scripts and two additional folders, namely, `utils` which contains scripts for custom functions and paramters used for plotting and `supplementary` that contains additional scripts for the results in the Supplementary Online Material in the study (Available here: [ADD LINK]). The scripts are numbered in the order in which they should be run.
 
* **01-prepare_data.R:** This script prepared the data files to be used in the analyses from the original download from the Nepture Database. It calculates the first appearance datum and last appearance datum of each species. It also categorises each occurrence into the tropics and extratropics, and calculate when dispersal from one zone to the other may have happened. 

* **02-analysis_diversity_dynamics.R:** This script calculates the diversity of marine plankton over the last 66 million years per climate zone and hemisphere.

* **03-analysis_LDG_div_overall.R:** This script calculates the latitudinal diversity gradients for marine plankton over the last 66 million years, divided into climate phases. 

* **04a-preferences_migration.R:** This script calculates origination, extinction and diversitfication preferences, and dispersal proportions across climate zones over the last 66 million years, divided into climate phases.

* **04b-plot_prefs_and_migration.R:** This script generates the plots of origination, extinction and diversitfication preferences, and dispersal proportions across climate zones over the last 66 million years, divided into climate phases, as per the analyses in `04a-preferences_migration.R`.

## Troubleshooting
The issue tracker is the preferred channel for bug reports. You may also contact me [by email](mailto:nussaibah.raja.schoob@fau.de).
