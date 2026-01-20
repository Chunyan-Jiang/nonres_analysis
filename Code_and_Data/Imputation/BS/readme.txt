# Imputation – Business Situation (BS)

This folder contains the imputation analysis for BS.
Three imputation methods LOCF, MC and RF are evaluated under a common simulation framework.

## Overview of files

- `LOCF_BS.Rmd`  
  Implements Last Observation Carried Forward (LOCF) imputation for BS.

- `MC_BS.Rmd`  
  Implements homogeneous Markov Chain (MC) imputation for BS.

- `RF_BS.Rmd`  
  Implements Random Forest (RF) imputation for BS.

- `Calibration_BS.Rmd`  
  Produces calibration analyses and diagnostic plots for the RF imputation.

- `BS_Summary.Rmd`  
  Summarizes imputation performance and generates comparison results across methods.

## Execution order

The scripts should be executed in the following order:

1. `LOCF_BS.Rmd`
2. `MC_BS.Rmd`
3. `RF_BS.Rmd`
4. `Calibration_BS.Rmd`
5. `BS_Summary.Rmd`

The first three scripts generate imputed datasets and evaluation results, which are subsequently used by the calibration and summary scripts.

## Reproducibility

For reproducibility and comparability across imputation methods, the random seed should not be changed unless explicitly stated.

## Runtime considerations

When applied to the full firm-level panel dataset, the MC and RF imputation procedures may require substantial runtime. This behavior is expected given the size of the data and the repeated simulation design.

## Inputs and outputs

The scripts assume that cleaned input data are available prior to execution.
All generated figures, tables, and intermediate results are written to subfolders within the `outputs/` directory located in this folder.
