# Imputation – Business Expectations (BE)

This folder contains the imputation analysis for BE.
Three imputation methods LOCF, MC and RF are evaluated under a common simulation framework.

## Overview of files

- `LOCF_BE.Rmd`  
  Implements Last Observation Carried Forward (LOCF) imputation for BE.

- `MC_BE.Rmd`  
  Implements homogeneous Markov Chain (MC) imputation for BE.

- `RF_BE.Rmd`  
  Implements Random Forest (RF) imputation for BE.

- `Calibration_BE.Rmd`  
  Produces calibration analyses and diagnostic plots for the RF imputation.

- `BE_Summary.Rmd`  
  Summarizes imputation performance and generates comparison results across methods.

## Execution order

The scripts should be executed in the following order:

1. `LOCF_BE.Rmd`
2. `MC_BE.Rmd`
3. `RF_BE.Rmd`
4. `Calibration_BE.Rmd`
5. `BE_Summary.Rmd`

The first three scripts generate imputed datasets and evaluation results, which are subsequently used by the calibration and summary scripts.

## Reproducibility

For reproducibility and comparability across imputation methods, the random seed should not be changed unless explicitly stated.

## Runtime considerations

When applied to the full firm-level panel dataset, the MC and RF imputation procedures may require substantial runtime. This behavior is expected given the size of the data and the repeated simulation design.

## Inputs and outputs

The scripts assume that cleaned input data are available prior to execution.
All generated figures, tables, and intermediate results are written to subfolders within the `outputs/` directory located in this folder.
