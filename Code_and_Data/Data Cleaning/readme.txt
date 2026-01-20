# Data Cleaning 
This folder contains the data cleaning process.

## Attention

The original data were provided by EBDC as test data. However, the test data differ from the real data in two aspects: (i) the variable vg_size is missing, and (ii) one variable is named inconsistently (vg_emplpast appears as vg_empl). To ensure structural consistency with the real data, we modified the original test data accordingly and generated a harmonized dataset, saved as data.dta in the Dataset directory. All subsequent analyses are based on this harmonized dataset (Dataset/data.dta).

## Overview of files

- `data_cleaning.Rmd`  
  
## Execution order

The scripts should be executed in the following order:

1. `data_cleaning.Rmd`

## Inputs and outputs

Input: Dataset/data.dta

Output:
- Dataset/unit_nonresponse_data.dta: from 1991 without participation number = 0
- Dataset/item_nonresponse_data.dta: from 1991
- Dataset/imputation_data.dta: from 2014
