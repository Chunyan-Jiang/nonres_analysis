PROJECT ELECTRONIC APPENDIX
===========================

This repository contains the electronic appendix for the project.
It provides the program code, test data, and outputs generated from the test data. 

1. Repository Structure
-----------------------

The repository is organized as follows:

Code_and_Data/
    Dataset/
        Data files used in the analyses.

    Data Cleaning/
        R scripts for cleaning and preparing the raw data.

    Descriptive Analysis/
        Code and outputs for descriptive analyses.

    Quantitative Analysis/
        Code and outputs for quantitative analyses.

    Imputation/
        Code and outputs for the imputation methods:
        - LOCF (Last Observation Carried Forward)
        - Markov Chain (MC)
        - Random Forest (RF)

report.pdf
    Final project report.

presentation_slides.pdf
    Presentation slides.


2. Code and Data
----------------

The folder "Code_and_Data/" contains all R scripts used to conduct the analyses,
create visualizations, and export results. Where relevant, outputs generated
from test data are included for verification purposes.


3. Execution Order
------------------

To reproduce the main results:

1) Run the scripts in the "Data_Cleaning/" folder to generate the cleaned dataset.

2) Run the scripts in the following folders as needed:
   - Descriptive_Analysis/
   - Quantitative_Analysis/
   - Imputation/

3) Within each subfolder, a README file documents the required execution order
   and any method-specific settings.

Note: Some imputation procedures (especially MC and RF) may require substantial
runtime due to repeated simulation.


4. Software Environment
-----------------------

The R session information used to run the analyses is documented in the file
"session_info.txt". 


5. Additional Notes
-------------------

- The code is written to be fully executable and reproducible.
- Random seeds are fixed within the scripts unless explicitly stated otherwise.
- Outputs may be overwritten if scripts are re-run.

For any questions regarding the repository or reproduction of results,
please refer to the documentation in the corresponding subfolders.
