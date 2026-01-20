# Imputation

This folder contains the imputation workflows used in the simulation study.

## Folder structure

- `BS/`  
  Contains scripts and outputs for imputation of BS
  This includes method-specific scripts (LOCF, MC, RF) simulation runs,
  calibration and summary scripts, as well as all corresponding output files
  (CSV results, figures, and intermediate RData objects).

- `BE/`  
  Contains scripts and outputs for imputation of BE
  The structure mirrors the `BS/` folder, with outcome-specific scripts and
  outputs organized in the same way.

- `helpers/`  
  Contains reusable helper functions used across BS and BE. These scripts define:
  - missingness simulation functions,
  - method-specific imputation functions (e.g. LOCF, MC, RF),
  - simulation and evaluation functions.
  
  The helper scripts are sourced by the method-specific scripts in `BS/` and
  `BE/` and are not intended to be executed standalone.

## Notes

Each script in `BS/` and `BE/` is designed to be runnable independently.
A common random seed is used to ensure reproducibility and comparability
across imputation methods.
Outputs are written to subfolders within the corresponding outcome directory.
