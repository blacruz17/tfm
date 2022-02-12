library(curatedMetagenomicData)
library(tidyverse)
library(janitor)
library(readxl)

library(phyloseq)
library(ape)
library(mia)
library(vegan)

library(MMUPHin)

library(caret)
library(pROC)

# Dataset preparation -------
metaphlanToPhyloseq <- source('DATA/metaphlanToPhyloseq.R')$value

# We will retrieve the data as in other scripts. Note that we first
# keep all samples, and we will filter out controls in later steps:
feng <- sampleMetadata %>%
  filter(study_name == "FengQ_2015")  %>%
  returnSamples("relative_abundance", counts = FALSE)

# Fix metadata object:
feng.df <- data.frame(colData(feng))
feng.df <- feng.df %>%
  select(where(~ !all(is.na(.x)))) %>%
  clean_names() 

# adds supplementary data from the publication (PMID 25758642):
# gender: 1 = male, 2 = female
feng.excel <- read_xlsx('DATA/metadata-fengq.xlsx', 
                        na = 'n.a.', skip = 3) %>% 
  clean_names() %>%
  rename(subject_id = x1) %>%
  mutate(subject_id = paste('SID', 
                            subject_id, sep = ''))
keep.df <- keep.df[-9] # no treatment column here

## joins both data frames:
feng.df.2 <- feng.df %>%
  select(all_of(keep.df)) %>%
  left_join(feng.excel, 'subject_id') %>%
  relocate(study_name, subject_id, ncbi_accession, pmid)

rownames(feng.df.2) <- feng.df.2$subject_id
colnames(feng.df.2) <- gsub('_1_yes_0_no', '', colnames(feng.df.2))

# adds obesity_status, metab and group columns:
feng.df.2 <- feng.df.2 %>%
  mutate(obesity_status = if_else(bmi >= 30, 
                                  'obese', 'lean', 'missing')) %>%
  mutate(metab = 
           if_else(tg_mg_l <= 150 &
                     fasting_glucose_mg_l <= 100 &
                     ((gender == 'male' & hdl_mg_l > 40)|
                        (gender == 'female' & hdl_mg_l > 50)) &
                     fatty_liver_in_ultrasound == 0 &
                     diabetes == 0 &
                     hypertension == 0 &
                     met_s == 0,
                   'healthy', 'unhealthy', 'missing')) %>%
  mutate(group = case_when(
    metab == 'healthy' & obesity_status == 'lean' ~ 'MHNO',
    metab == 'unhealthy' & obesity_status == 'lean' ~ 'MUNO',
    metab == 'healthy' & obesity_status == 'obese' ~ 'MHO',
    metab == 'unhealthy' & obesity_status == 'obese' ~ 'MUO'
  )) %>%
  mutate_at(vars(study_condition, obesity_status, 
                 metab, group), as.factor)

# now back to the treeSummarizedExperiment:
identical(rownames(DataFrame(feng.df.2)), rownames(colData(feng)))
colData(feng) <- DataFrame(feng.df.2)

# Preprocessing -----
pseq.feng <- makePhyloseqFromTreeSummarizedExperiment(feng,
                                                      "relative_abundance")
pseq.feng <- tax_glom(pseq.feng, 'Species')

# set abundances to numbers between 0-1:
otu_table <- pseq.feng@otu_table %>% 
              apply(2, function(x) x / 100)

# fix the rownames:
otu_table <- t(data.frame(t(otu_table)) %>% clean_names())
rownames(otu_table) <- gsub('k_.*_s_', 's_', rownames(otu_table))

# import model and training dataset:
# we need the training dataset to know which variables are missing
# in this dataset, we will add these variables and fill them
# with 0s (that will be converted to pseudocounts)
model <- readRDS('DATA/randomForest/bestmodel.rds')
trainDescr <- readRDS('DATA/randomForest/trainDescr_bestmodel.rds')

# detect missing variables:
missing <- t(trainDescr[1:ncol(otu_table),
                      !colnames(trainDescr) %in% rownames(otu_table)])
# fill with 0s:
missing[ , ] <- 0

otu_table <- rbind(otu_table, missing)

feng_pseq <- metaphlanToPhyloseq(otu_table,
                                         metadat = pseq.feng@sam_data)
# CLR-transform (automatically adds pseudocounts):
feng_pseq_clr <- microbiome::transform(feng_pseq, "clr")

# Prepare data for predictions -------
totu <- as.data.frame(t(feng_pseq_clr@otu_table))
totu$SID <- rownames(totu)

tags <- data.frame(feng_pseq_clr@sam_data)
tags$SID <- rownames(tags)

tags <- tags %>% 
  select(SID, group, study_condition, age_yrs, gender, bmi)

totu <- totu %>% left_join(tags, by = 'SID')
rownames(totu) <- totu$SID

data.RF <- totu %>% 
  select(-c(SID)) %>%
  rename(class = group) %>%
  clean_names() 

########## filters out controls:
data.RF <- data.RF %>%
            filter(study_condition != 'control')

class.RF <- data.RF %>%
  select(class, study_condition)

data.RF <- data.RF %>%
  select(-c(class)) %>%
  select(-c(study_condition, age_yrs, gender, bmi))

# Predict labels on the dataset ------------
class.RF$bi <- factor(gsub('O|NO', '', class.RF$class))

predicted <- predict(model, data.RF, type = 'prob') %>%
  mutate('class'= names(.)[apply(., 1, which.max)],
         'real' = class.RF$bi)

predicted_roc <- roc(class.RF$bi, predicted$MU)

# save the ROC curve for plotting:
saveRDS(predicted_roc, 'DATA/randomForest/roc_feng.rds')