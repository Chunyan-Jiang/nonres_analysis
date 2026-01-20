# Quantitive Analysis

This folder contains Quantitive Analysis for unit nonresponse.


## Overview of files

- `quan_crgam.Rmd`  
  Generate Model with state (in paper) and model with Westeast (optional)

- `quan_crgam_outputs.Rmd`  
  Generate the results of model (summary, gamcheck,concurvity, cluster-robust     covariance matrix, cluster-robust coefficient tests, forest plot and spline

- `quan_effective_G.Rmd`  
  Calculate effective G:serves as a diagnostic for assessing whether the effective number of independent clusters is sufficiently large for reliable cluster-robust inference.

- helper-folder
  The function we need in analysis

## Execution order

The scripts should be executed in the following order:

1. `quan_crgam.Rmd`
2. `quan_crgam_outputs.Rmd` 
3. `quan_effective_G.Rmd`


In quan_crgam_outputs.Rmd, cluster-robust covariance matrices are saved as .rds files in the outputs directory. These .rds files are subsequently read back into quan_crgam_outputs.Rmd and used for plotting.


## Inputs and outputs

The scripts assume that cleaned input data are available prior to execution.
All generated figures, tables, and intermediate results are written to subfolders within the `outputs/` directory located in this folder.
