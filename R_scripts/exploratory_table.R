library(curatedMetagenomicData)
library(dplyr)
library(tidyr)
library(gtsummary)
library(gt)

library(janitor)
library(readxl)

library(phyloseq)
library(ape)
library(mia)

metaphlanToPhyloseq <- source('DATA/metaphlanToPhyloseq.R')$value

# Data preparation -------------------
## This first section is the same as in the alpha_and_beta_diversity.R script
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

## get dataframes:
df.feng <- data.frame(pseq.feng@sam_data)
df.karl <- data.frame(pseq.karl@sam_data)
df.t2d <- data.frame(pseq.t2d@sam_data)
df.mh <- data.frame(pseq.mh@sam_data)

# Preparing dataframes ------
## some aesthetic fixes so that the final table looks a bit better:
df.feng <- df.feng %>% 
              rename(age = age_yrs) %>%
              mutate(sex = if_else(gender == 'male', 'Male', 'Female')) %>%
              mutate_at(vars(sex), as.factor) %>%
              mutate(fasting_glucose_mmol_l = fasting_glucose_mg_l * 0.0555,
                     hdl_mmol_l = hdl_mg_l / 38.67 ,
                     triglycerides_mmol_l = tg_mg_l / 88.57 )

df.karl <- df.karl %>% 
              rename(age = age_years) %>%
              mutate(sex = if_else(gender == 'male', 'Male', 'Female')) %>%
              mutate_at(vars(sex), as.factor)

df.t2d <- df.t2d %>%
              rename(country = cohort_country) %>%
              mutate(sex = if_else(gender == "M", "Male", "Female"),
                     study_name = "MetaHit_T2D")

df.mh <- df.mh %>%
              mutate(sex = if_else(gender == "male", "Male", "Female"),
                     study_name = "MetaHit_Healthy")

df.metahit <- df.mh %>%
              bind_rows(df.t2d)

## Summary table ----
df.all <- df.feng %>%
          mutate(country = "Austria") %>%
          bind_rows(df.karl %>% mutate(country = "Sweden"), 
                    df.metahit %>% mutate(country = "Denmark"))

save.metadata <- df.all %>% 
                  select(study_name, group, metab, obesity_status, 
                         sex, age,  bmi, tg_mg_l, study_condition, 
                         hdl_mg_l, fasting_glucose_mg_l) %>% 
                  tibble::rownames_to_column('sample_id')
write.table(save.metadata, 'DATA/metadata.txt', sep = '\t',
            row.names = FALSE, col.names = TRUE, quote = FALSE)
			
mytable <- df.all %>%
            select(group, age, sex, bmi,
                   hdl_mmol_l, fasting_glucose_mmol_l, triglycerides_mmol_l) %>%
            tbl_summary(by = group,
                        label = list(age ~ 'Age (years)',
                                     sex ~ 'Sex',
                                     bmi ~ 'BMI (kg/m<sup>2</sup>)',
                                     hdl_mmol_l ~ 'HDL (mmol/L)',
                                     fasting_glucose_mmol_l ~ 'Fasting glucose (mmol/L)',
                                     triglycerides_mmol_l ~ 'Triglycerides (mmol/L)'),
                        missing = 'no') %>%
            add_n() %>%
            add_p() %>%
            add_q() %>%
            modify_header(label ~ '**Variable**') %>%
            modify_spanning_header(list(c("stat_1", "stat_3") ~ "**Lean**",
                                         c("stat_2", "stat_4") ~ "**Obese**")) %>%
            bold_labels()

mytable

# gtsave(as_gt(mytable) %>% 
#              gt::fmt_markdown(columns = vars(label)),
#        "plots/table-exploratory.png", zoom = 1)