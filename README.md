# Bioinformatic methods for stratification of obese patients and identification of cancer susceptibility biomarkers based on the analysis of gut microbiota

Scripts used to analyze the metagenomic data, build the random forest classifier and the patient similarity network, and perform the LEfSe analysis on CRC data.

## Table of contents
* [Data directory](#data-directory)
* [Taxonomic and functional profiling](#taxonomic-and-functional-profiling)
* [R scripts](#R-scripts)
* [CRC LEfSe](#CRC-LEfSe)

## Data directory
This directory contains the scripts used to preprocess the raw reads and for taxonomic and functional profiling.
- `samples_from_ENA.txt`: text file containing the project ID, the sample ID and the subject ID for all samples downloaded from the European Nucleotide Archive at EMBL-EBI
- `adjusted_abundances.txt`: relative abundances table for the 356 samples after batch effect adjustment

**NOTE:** the celiac disease cohort and the WGS CRC cohort are unavailable, as the data have not been published.

## Taxonomic and functional profiling
The `process_data_from_ena` directory contains:
- `metawrap_and_metaphlan.sh`: bash script used for preprocessing and for obtaining the relative abundances table
  - `metaphlan-script.pl`: script called by the `metawrap_and_metaphlan.sh` pipeline to run MetaPhlAn on the raw reads
  - `concatenator.sh`: helper script that concatenates reads from the same subject into a single file
- `run_humann.sh`: script for functional profiling 
  - `cat_files.pl`: script called by `run_humann.sh` to concatenate forward and reverse read files together
 
## R scripts
The `R_scripts` directory contains the R files for the alpha and beta diversity analyses, the random forest classifier, and the patient similarity network.
- `alpha_and_beta_diversity.R`: performs alpha and beta diversity analyses, as well as batch effect correction; and plots the corresponding figures
- for the **random forest**:
  - `randomForest_preprocessing.R`: filters species with low abundances, correlated features and performs CLR transformation on the data
  - `randomForest_tuning.R`: tries different configurations for the model and obtains hyperparameter tuning and ROC curves figures 
  - `randomForest_finalModel.R`: trains the final model, plots the ROC curves and the heatmap
  - `randomForest_validation.R`: validates the data on the CRC dataset (validation on the celiac disease cohort follows the same structure)

## CRC LEfSe
The `LEfSe_CRC` directory contains two files:
- `prepare_lefse_input.R`: R script to transform the relative abundances table into an input file appropriate for LEfSe
- `running_lefse.sh`: commands used for running LEfSe on the conda version 1.0.0 and for obtaining the corresponding figures
