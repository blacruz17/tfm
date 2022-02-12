library(curatedMetagenomicData)
library(dplyr)
library(tidyr)

library(ggplot2)
library(ggpubr)
library(RColorBrewer)

library(janitor)
library(readxl)

library(phyloseq)
library(ape)
library(mia)
library(vegan)

library(rstatix)
library(MMUPHin)
library(scater)

metaphlanToPhyloseq <- source('DATA/metaphlanToPhyloseq.R')$value

# Preparing the data -------
## T2D data from MetaHit ----
# relative abundances table:
abundances.t2d <- read.table('DATA/abundances_t2d.txt',
                             row.names = 1)
# metadata:
sample.t2d <- read.table('DATA/sample_t2d.txt', row.names = 1)

pseq.t2d <- metaphlanToPhyloseq(abundances.t2d, metadat = sample.t2d,
                                simplenames = FALSE)

random_tree <- rtree(ntaxa(pseq.t2d), rooted = T, 
                     tip.label = taxa_names(pseq.t2d))

pseq.t2d <- merge_phyloseq(pseq.t2d, random_tree)

## Healthy data from MetaHit ----
# relative abundances:
abundances.mh <- read.table('DATA/abundances_mh.txt',
                            row.names = 1)
# metadata:
sample.mh <- read.table('DATA/sample_mh.txt', row.names = 1)

pseq.mh <- metaphlanToPhyloseq(abundances.mh, metadat = sample.mh,
                               simplenames = FALSE)

random_tree <- rtree(ntaxa(pseq.mh), rooted = T, 
                     tip.label = taxa_names(pseq.mh))

pseq.mh <- merge_phyloseq(pseq.mh, random_tree)

## Data from curatedMetagenomicData ----
### Karlsson ----
karlsson <- sampleMetadata %>%
  filter(study_name == "KarlssonFH_2013")  %>%
  returnSamples("relative_abundance", counts = FALSE)

## metadata fix:
## turns to dataframe, cleans names and removes empty columns:
karl.df <- data.frame(colData(karlsson)) %>% 
  clean_names() %>%
  select(where(~ !all(is.na(.x))))

## imports metadata from the publication
## the suppdata-karlsson file can be obtained from the online
## version of the paper (PMID 23719380)
karl.metadata.1 <- read_xlsx('DATA/suppdata-karlsson.xlsx',
                             sheet = 'metadata1') %>%
  clean_names()

karl.metadata.2 <- read_xlsx('DATA/suppdata-karlsson.xlsx',
                             sheet = 'metadata2') %>%
  clean_names()

## joins and generates subject_id column
## for joining with karl.df:
karl.metadata <- karl.metadata.1 %>% 
  left_join(karl.metadata.2, by = "sample_id") %>%
  mutate(subject_id = paste0('S', sample_id)) %>%
  select(-sample_id)

# converts subject_id to rownames while maintaining the column:
rownames(karl.metadata) <- karl.metadata$subject_id

# columns from curatedMetagenomicData that we will keep:
keep.df <- c('age_category', 'gender', 'ncbi_accession', 'pmid', 'subject_id',
             'sequencing_platform', 'study_condition','study_name','treatment')

# joins, selects columns of interest,
# renames bmi column, relocates some columns:
karl.df.2 <- karl.df %>%
  select(all_of(keep.df)) %>%
  left_join(karl.metadata, 'subject_id') %>%
  rename(bmi = bmi_kg_m2) %>% 
  relocate(study_name, subject_id, ncbi_accession, pmid)

rownames(karl.df.2) <- karl.df.2$subject_id

# adds obesity_status, metab and group (MUO/MHO/MUNO/MHNO) columns:
karl.df.2 <- karl.df.2 %>%
  select(where(~ !all(is.na(.x)))) %>%
  clean_names() %>%
  mutate(obesity_status = if_else(bmi >= 30, 
                                  'obese', 'lean', 'missing')) %>%
  mutate(metab = if_else(study_condition == 'control' &
                           triglycerides_mmol_l <= 1.7 &
                           fasting_glucose_mmol_l <= 6.1 &
                           hdl_mmol_l > 1.3,
                         'healthy', 'unhealthy', 'missing')) %>%
  mutate(group = case_when(
    metab == 'healthy' & obesity_status == 'lean' ~ 'MHNO',
    metab == 'unhealthy' & obesity_status == 'lean' ~ 'MUNO',
    metab == 'healthy' & obesity_status == 'obese' ~ 'MHO',
    metab == 'unhealthy' & obesity_status == 'obese' ~ 'MUO'
  ))

# there is one NA / 'missing':
karl.df.2 %>% 
  filter(metab == 'missing') %>%
  select(subject_id, study_condition, bmi,
         triglycerides_mmol_l, fasting_glucose_mmol_l, hdl_mmol_l, 
         obesity_status, metab, group)

# we can see that except for glucose, she looks healthy, so:
karl.df.2[karl.df.2$metab == 'missing', ] # row S58
karl.df.2['S58', 'metab'] <- 'healthy'
karl.df.2['S58', 'group'] <- 'MHO'

# converts variables to factors:
# if we did this earlier we would get 'missing' as a level
# and we do not want that
karl.df.2 <- karl.df.2 %>%
  mutate_at(vars(study_condition, obesity_status, 
                 metab, group), as.factor)

# checks everything is OK (should be TRUE):
identical(rownames(DataFrame(karl.df.2)), rownames(colData(karlsson)))
# now back to the treeSummarizedExperiment:
colData(karlsson) <- DataFrame(karl.df.2)

### Feng ----
feng <- sampleMetadata %>%
  filter(study_name == "FengQ_2015")  %>%
  filter(study_condition == "control") %>%
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

# Joins datasets  -----
## Coerce all data to phyloseq objects:
pseq.feng <- makePhyloseqFromTreeSummarizedExperiment(feng,
                                                      "relative_abundance")
pseq.karl <- makePhyloseqFromTreeSummarizedExperiment(karlsson,
                                                      "relative_abundance")

# agglomerate to species level:
pseq.feng <- tax_glom(pseq.feng, 'Species')
pseq.karl <- tax_glom(pseq.karl, 'Species')
pseq.mh <- tax_glom(pseq.mh, 'Species')
pseq.t2d <- tax_glom(pseq.t2d, 'Species')

pseq.mh@sam_data$study_name <- 'MetaHit_MH'
pseq.t2d@sam_data$study_name <- 'MetaHit_T2D'

### merging ---
merged_otu_table <- t(dada2::mergeSequenceTables(t(pseq.feng@otu_table), 
                                                 t(pseq.karl@otu_table),
                                                 t(pseq.mh@otu_table),
                                                 t(pseq.t2d@otu_table)))

merged_metadata <- rbind(data.frame(pseq.feng@sam_data) %>% 
                           rename(age = age_yrs) %>%
                           select(group, metab, obesity_status, 
                                  study_name, age, gender, bmi),
                         data.frame(pseq.karl@sam_data) %>% 
                           rename(age = age_years) %>%
                           select(group, metab, obesity_status, 
                                  study_name, age, gender, bmi),
                         data.frame(pseq.mh@sam_data) %>% 
                           select(group, metab, obesity_status, 
                                  study_name, age, gender, bmi),
                         data.frame(pseq.t2d@sam_data) %>% 
                           select(group, metab, obesity_status, 
                                  study_name, age, gender, bmi)) %>%
  mutate(sex = case_when(
    gender == "M" ~ "male",
    gender == "F" ~ "female",
    gender == 'male' ~'male',
    gender == 'female' ~'female'
  ))

merged_metadata <- merged_metadata %>% 
  mutate_at(vars(group, metab, obesity_status,
                 study_name, sex), as.factor)

# batch effect adjustment:
merged_otu <- merged_otu_table %>%
                apply(2, function(x) x / 100)

fit_adjust_batch <- adjust_batch(feature_abd = merged_otu,
                                 batch = "study_name",
                                 covariates = "group",
                                 data = merged_metadata,
                                 control = list(verbose = FALSE))

adjusted_otu_table <- fit_adjust_batch$feature_abd_adj

# Preprocessing the relative abundances table ----------- 
## removes samples whose relative abundance is < 0.05
## in more than half of the samples:
a <- apply(adjusted_otu_table, 1, function(x) x < 0.0005)
b <- apply(a, 2, function(x) sum(x) > dim(adjusted_otu_table)[2]/2)
filtered_otu_table <- adjusted_otu_table[b, ]

## CLR transformation:
filtered_phyloseq <- metaphlanToPhyloseq(filtered_otu_table,
                                         metadat = merged_metadata,
                                         simplenames = FALSE)
## this function adds pseudocounts automatically:
filtered_phyloseq.clr <- microbiome::transform(filtered_phyloseq, "clr")

# Obtaining class labels ----------------------
totu <- as.data.frame(t(filtered_phyloseq.clr@otu_table))
totu$SID <- rownames(totu)

tags <- merged_metadata %>%
        mutate(SID = rownames(merged_metadata)) %>%
        select(SID, group, age, sex, bmi)

totu <- totu %>% 
          left_join(tags, by = 'SID')
rownames(totu) <- totu$SID

data.RF <- totu %>% 
              select(-c(SID)) %>%
              rename(class = group) %>% 
              clean_names() 

class.RF <- data.RF %>%
              select(class)

data.RF <- data.RF %>%
            select(-c(class))

data.RF_nometa <- data.RF %>%
                  select(-c(age, sex, bmi))

## save input files for the Random Forest:
write.table(class.RF, file = 'DATA/randomForest/classRF.tsv',
            quote = FALSE, sep = '\t')
write.table(data.RF_nometa, file = 'DATA/randomForest/dataRF.tsv',
            quote = FALSE, sep = '\t')
write.table(data.RF, file = 'DATA/randomForest/dataRF_metadata.tsv',
            quote = FALSE, sep = '\t')
