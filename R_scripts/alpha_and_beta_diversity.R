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
library(MMUPHin)

library(rstatix)
library(scater)

# Data preparation -------------------
# the metaphlanToPhyloseq script can be obtained at:
# https://github.com/waldronlab/presentations/blob/master/Waldron_2016-06-07_EPIC/metaphlanToPhyloseq.R
metaphlanToPhyloseq <- source('DATA/metaphlanToPhyloseq.R')$value

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

# agglomerate to species level:
pseq.feng <- tax_glom(pseq.feng, 'Species')
pseq.karl <- tax_glom(pseq.karl, 'Species')
pseq.mh <- tax_glom(pseq.mh, 'Species')
pseq.t2d <- tax_glom(pseq.t2d, 'Species')

pseq.mh@sam_data$study_name <- 'MetaHit_MH'
pseq.t2d@sam_data$study_name <- 'MetaHit_T2D'

# approximate to counts:
# this step is necessary to calculate alpha diversities
pseq.feng <- transform_sample_counts(pseq.feng, function(x) 1E6 * x)
pseq.karl <- transform_sample_counts(pseq.karl, function(x) 1E6 * x)
pseq.t2d <- transform_sample_counts(pseq.t2d, function(x) 1E6 * x)
pseq.mh <- transform_sample_counts(pseq.mh, function(x) 1E6 * x)

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

merged_phyloseq <- metaphlanToPhyloseq(merged_otu_table,
                                       metadat = merged_metadata,
                                       simplenames = FALSE)

### we can visualize this with a PCA plot:
## helper function:
pca_plot <- function(pseq) {
  pseq_ord <- ordinate(pseq, 'RDA')
  p <- plot_ordination(pseq, pseq_ord, 
                       type = "samples", color = "study_name") + 
    geom_point(size = 2) +
    stat_ellipse(aes(group = study_name), linetype = 2) +
    theme_bw() +
    theme(plot.title = element_text(face = 'bold', size = 16),
          axis.title = element_text(size = 14),
          legend.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 12))
  return(p)
}

p_before <- pca_plot(merged_phyloseq)

## MMUPHin batch effect correction -----------
# converts to numbers 0-1:
merged_otu <- merged_otu_table %>% 
                  apply(2, function(x) round(x))

# batch effect adjustment:
fit_adjust_batch <- adjust_batch(feature_abd = merged_otu,
                                 batch = "study_name",
                                 covariates = "group",
                                 data = merged_metadata,
                                 control = list(verbose = FALSE))

adjusted_otu_table <- fit_adjust_batch$feature_abd_adj

adjusted_phyloseq <- metaphlanToPhyloseq(adjusted_otu_table,
                                         metadat = merged_metadata,
                                         simplenames = FALSE)

random_tree <- rtree(ntaxa(adjusted_phyloseq), rooted = T, 
                     tip.label = taxa_names(adjusted_phyloseq))

adjusted_phyloseq <- merge_phyloseq(adjusted_phyloseq, random_tree)

p_after <- pca_plot(adjusted_phyloseq)

### PCA plot ----
ggarrange(p_before + ggtitle('a'), p_after + ggtitle('b'), 
          common.legend = TRUE, legend = 'bottom')

### PERMANOVA ----
# this will allow us to see if the effect of study_name was reduced:
D_before <- vegdist(t(merged_otu))
D_after <- vegdist(t(adjusted_otu_table))

set.seed(1) # for reproducibility
fit_adonis_before <- adonis(D_before ~ study_name, 
                            data = merged_metadata)
fit_adonis_after <- adonis(D_after ~ study_name, 
                           data = merged_metadata)
print(fit_adonis_before)
print(fit_adonis_after)

# the treeSummarizedExperiment object will let us 
# calculate alpha diversity measures with mia:
adjusted_phyloseq.tse <- makeTreeSummarizedExperimentFromphyloseq(adjusted_phyloseq)

# Alpha diversity ---------------
## Study by study ------
mh <- makeTreeSummarizedExperimentFromphyloseq(pseq.mh)
t2d <- makeTreeSummarizedExperimentFromphyloseq(pseq.t2d)
feng <- makeTreeSummarizedExperimentFromphyloseq(pseq.feng)
karlsson <- makeTreeSummarizedExperimentFromphyloseq(pseq.karl)

# calculates richness, diversity, evenness and dominance:
tse.list <- list("mh" = mh, "t2d" = t2d, 
                 "karl" = karlsson, "feng" = feng)

tse.list <- lapply(names(tse.list), function(x){
                    tse.list[[x]] <- tse.list[[x]] %>%
                                estimateRichness(abund_values = "counts",
                                                 index = "observed",
                                                 name = "observed") %>%
                                estimateRichness(abund_values = "counts",
                                                 index = "chao1",
                                                 name = "chao1") %>%
                                estimateDiversity(abund_values = "counts", 
                                                  index = "shannon", 
                                                  name = "shannon") %>%
                                estimateDiversity(abund_values = "counts", 
                                                  index = "gini_simpson", 
                                                  name = "simpson") %>%
                                estimateEvenness(abund_values = "counts", 
                                                 index = "pielou", 
                                                 name = "pielou") %>%
                                estimateDominance(abund_values = "counts", 
                                                  index = "relative", 
                                                  name = "relative") %>%
                                estimateDiversity(abund_values = "counts", 
                                                  index = "log_modulo_skewness", 
                                                  name = "log_modulo_skewness")
                  })

names(tse.list) <- c('MetaHit Healthy', 'MetaHit T2D', 
                     'Karlsson FH (2013)', 'Feng Q (2015)')

group.palette <- brewer.pal(4, "Set2")
names(group.palette) <- c('MHNO', 'MUNO', 'MUO', 'MHO')

# helper function to build plots:
plotter <- function(study_name){
  
  df <- as.data.frame(colData(tse.list[[study_name]])) %>%
    select(chao1, shannon, simpson, 
           pielou, relative, log_modulo_skewness,
           group, metab, obesity_status)
  
  # Kruskal-Wallis test:
  stat.test <- df %>%
    kruskal_test(chao1 ~ group) %>%
    bind_rows(kruskal_test(df, shannon ~ group),
              kruskal_test(df, simpson ~ group),
              kruskal_test(df, pielou ~ group),
              kruskal_test(df, relative ~ group),
              kruskal_test(df, log_modulo_skewness ~ group)) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance() %>%
    mutate(p.adj.signif = gsub('ns', '', p.adj.signif))
  
  # Plots:
  counter <<- 0
  plist <- lapply(df[ ,c("chao1", "shannon", "simpson",
                         "pielou", "relative", "log_modulo_skewness")], 
                  function(a) 
                  {counter <<- counter + 1
                  ggplot(df, aes(x = group, y = a)) +
                    geom_boxplot(aes(fill = group), 
                                 alpha=.5,
                                 outlier.shape = NA) +
                    scale_fill_manual(values = group.palette) +
                    geom_jitter(width = 0.2,
                                aes(colour = group), size = 1.5) +
                    scale_color_manual(values = group.palette) +
                    geom_text(x = 1, y = min(a),
                              label = stat.test$p.adj.signif[counter]) +
                    ylab(gsub('_', ' ', colnames(df)[counter])) +
                    xlab(NULL) +
                    theme_bw() +
                    theme(axis.text = element_text(size = 12),
                          axis.title = element_text(size = 14),
                          legend.text = element_text(size = 12),
                          legend.title = element_text(size = 14),
                          axis.text.x = element_text(angle = 45, hjust=1))}
  ) 
  return(plist)
}

# gets plots for all studies:
plots <- lapply((names(tse.list)), 
                function(i) {tmp <- plotter(i)
                ggarrange(plotlist = tmp,
                          common.legend = TRUE, legend = "right",
                          align = "hv") %>%
                  annotate_figure(top = paste0(i,': Alpha diversity'))
                })

# there are no significant differences:
plots[[1]]
plots[[2]]
plots[[3]]
plots[[4]]

## With merged + adjusted OTU table -----
adjusted_phyloseq.tse <- adjusted_phyloseq.tse %>%
                            estimateRichness(abund_values = "counts",
                                             index = "observed",
                                             name = "observed") %>%
                            estimateRichness(abund_values = "counts",
                                             index = "chao1",
                                             name = "chao1") %>%
                            estimateDiversity(abund_values = "counts", 
                                              index = "shannon", 
                                              name = "shannon") %>%
                            estimateDiversity(abund_values = "counts", 
                                              index = "gini_simpson", 
                                              name = "simpson") %>%
                            estimateEvenness(abund_values = "counts", 
                                             index = "pielou", 
                                             name = "pielou") %>%
                            estimateDominance(abund_values = "counts", 
                                              index = "relative", 
                                              name = "relative") %>%
                            estimateDiversity(abund_values = "counts", 
                                              index = "log_modulo_skewness", 
                                              name = "log_modulo_skewness")

## runs kruskal-wallis test:
df <- as.data.frame(colData(adjusted_phyloseq.tse)) %>%
                    select(chao1, shannon, simpson,
                           pielou, relative, log_modulo_skewness,
                           group, metab, obesity_status, study_name)

kruskal.test <- df %>%
  kruskal_test(chao1 ~ group) %>%
  bind_rows(kruskal_test(df, shannon ~ group),
            kruskal_test(df, simpson ~ group),
            kruskal_test(df, pielou ~ group),
            kruskal_test(df, relative ~ group),
            kruskal_test(df, log_modulo_skewness ~ group)) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>%
  mutate(p.adj.signif = gsub('ns', '', p.adj.signif))

## gets plots:
counter <- 0
plist <- lapply(df[ ,c("chao1", "shannon", "simpson",
                       "pielou", "relative", "log_modulo_skewness")], 
                function(a) 
                {counter <<- counter + 1
                ggplot(df, aes(x = group, y = a)) +
                  geom_boxplot(aes(fill = group), 
                               alpha=.5,
                               outlier.shape = NA) +
                  scale_fill_manual(values = group.palette) +
                  geom_jitter(width = 0.2,
                              aes(colour = group), size = 1.5) +
                  scale_color_manual(values = group.palette) +
                  
                  ylab(gsub('_', ' ', colnames(df)[counter])) +
                  xlab(NULL) +
                  theme_bw() +
                  theme(axis.text = element_text(size = 12),
                        axis.title = element_text(size = 14),
                        legend.text = element_text(size = 12),
                        legend.title = element_text(size = 14),
                        axis.text.x = element_text(angle = 45, hjust=1))}
) 

ggarrange(plotlist = plist,
          common.legend = TRUE, legend = "right",
          align = "hv")

# chao and pielou panels do not need modifications:
chao_panel <- plist[[1]]
pielou_panel <- plist[[4]]

# runs Wilcoxon on the remaining measurements:
wil.test <- bind_rows(df %>%
                        pairwise_wilcox_test(shannon ~ group) ,
                      df %>%
                        pairwise_wilcox_test(simpson ~ group),
                      df %>%
                        pairwise_wilcox_test(relative ~ group),
                      df %>%
                        pairwise_wilcox_test(log_modulo_skewness ~ group)) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

# plots significances:
signif.p <- plist[-c(1, 4)]
tests <- wil.test[-c(1, 4)]

signif.list <- list(
  signif.p[[1]] + stat_pvalue_manual(
    wil.test %>% filter(.y. == names(signif.p)[[1]]), 
    label = "p.adj.signif", 
    hide.ns = TRUE,
    tip.length = 0.01,
    y.position = 4,
    bracket.shorten = 0.05,
    size = 5, face = 'bold'),
  
  signif.p[[2]] + stat_pvalue_manual(
    wil.test %>% filter(.y. == names(signif.p)[[2]]), 
    label = "p.adj.signif", 
    hide.ns = TRUE,
    tip.length = 0.01,
    y.position = 1,
    bracket.shorten = 0.05,
    size = 5, face = 'bold'),
  
  signif.p[[3]] + stat_pvalue_manual(
    wil.test %>% filter(.y. == names(signif.p)[[3]]), 
    label = "p.adj.signif", 
    hide.ns = TRUE,
    tip.length = 0.01,
    y.position = 1,
    bracket.shorten = 0.05,
    size = 5, face = 'bold'),
  
  signif.p[[4]] + stat_pvalue_manual(
    wil.test %>% filter(.y. == names(signif.p)[[4]]), 
    label = "p.adj.signif", 
    hide.ns = TRUE,
    tip.length = 0.01,
    y.position = 2.07,
    bracket.shorten = 0.05,
    size = 5, face = 'bold')
)

# some modifications for the final figure:
chao_panel <- chao_panel +
              theme(legend.title = element_blank(),
                    axis.text.x = element_blank())
shannon_panel <- signif.list[[1]] +
                 theme(axis.text.x = element_blank())
simpson_panel <- signif.list[[2]] +
                 theme(axis.text.x = element_blank())
relative_panel <- signif.list[[3]]
log_panel <- signif.list[[4]]

ggarrange(chao_panel, shannon_panel, simpson_panel,
          pielou_panel, relative_panel, log_panel, 
          common.legend = TRUE, legend = 'right',
          align = 'v', heights = c(0.45, 0.55))

## Obesity and metabolic health separately ----
## rewrites helper function for this section:
plotter <- function(study_df, var, var.palette, hide.ns = FALSE){
  
  # Kruskal-Wallis:
  stat.test <- study_df %>%
    kruskal_test(chao1 ~ var) %>%
    bind_rows(kruskal_test(study_df, shannon ~ var),
              kruskal_test(study_df, simpson ~ var),
              kruskal_test(study_df, pielou ~ var),
              kruskal_test(study_df, relative ~ var),
              kruskal_test(study_df, log_modulo_skewness ~ var)) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  
  ## PLOTS
  counter <<- 0
  if (hide.ns == TRUE){stat.test$p.adj.signif <- gsub('ns', '', stat.test$p.adj.signif)}
  plist <- lapply(study_df[ ,c("chao1", "shannon", "simpson",
                               "pielou", "relative", "log_modulo_skewness")], 
                  function(a) 
                  {counter <<- counter + 1
                  ggplot(study_df, aes(x = var, y = a)) +
                    geom_boxplot(aes(fill = var), 
                                 alpha=.5,
                                 outlier.shape = NA) +
                    scale_fill_manual(values = var.palette) +
                    geom_jitter(width = 0.2,
                                aes(colour = var), size = 1.5) +
                    scale_color_manual(values = var.palette) +
                    geom_text(x = 1, y = min(a),
                              label = stat.test$p.adj.signif[counter]) +
                    ylab(gsub('_', ' ', colnames(study_df)[counter])) +
                    xlab(NULL) +
                    theme_bw() +
                    theme(legend.title = element_blank()) +
                    theme(plot.title = element_text(face = 'bold', size = 16),
                          axis.title = element_text(size = 14),
                          axis.text = element_text(size = 12),
                          legend.text = element_text(size = 12),
                          axis.text.x = element_text(angle = 45, hjust=1))}
  ) 
  return(plist)
}

df_all <- as.data.frame(colData(adjusted_phyloseq.tse)) %>%
  select(chao1, shannon, simpson,
         pielou, relative, log_modulo_skewness,
         group, metab, obesity_status, study_name)

### Obesity -----
obesity.palette <- brewer.pal(4, "Paired")[3:4]
names(obesity.palette) <- c('lean', 'obese')

all_plots_obesity <- plotter(df_all, df_all$obesity_status, obesity.palette,
                             hide.ns = TRUE)

ob_panel <- ggarrange(plotlist = all_plots_obesity,
                      common.legend = TRUE, legend = "bottom",
                      align = "hv")
### Metabolic health -----
metab.palette <- brewer.pal(4, "Paired")[1:2]
names(metab.palette) <- c('healthy', 'unhealthy')

all_plots_metab <- plotter(df_all, df_all$metab, metab.palette, hide.ns = TRUE)

metab_panel <- ggarrange(plotlist = all_plots_metab,
                          common.legend = TRUE, legend = "bottom",
                          align = "hv")

ggarrange(ob_panel + 
            ggtitle('a') + 
            theme(plot.title = element_text(size = 16, 
                                            face = 'bold', hjust = 0.1)),
          metab_panel + 
            ggtitle('b') + 
            theme(plot.title = element_text(size = 16, 
                                            face = 'bold', hjust = 0.1)),
          nrow = 2, hjust = 'hv')

# Beta diversities --------------
## Helper functions ------
# Calculates distance matrix and PCoA:
distances <- function(study_pseq) {
  for( i in dist_methods ){
    # Calculate distance matrix
    iDist <- distance(study_pseq, method = i)
    dlist[[i]] <- iDist
    
    # Calculate ordination
    iPCoA  <- ordinate(study_pseq, "PCoA", distance = iDist)
    pcoa_list[[i]] <- iPCoA
  }
  to_return <- list('dlist' = dlist, 'pcoa_list' = pcoa_list)
  return(to_return)
}

# Performs PERMANOVA:
adonis_calculator <- function(dlist, study_pseq) {
  
  results.group <- lapply(names(dlist), 
                        function(x) {
                          z <- adonis(dlist[[x]] ~ group, 
                                      data = data.frame(sample_data(study_pseq)))
                          return(as.data.frame(z$aov.tab))
                        })
  names(results.group) <- names(dlist)
  return(results.group)
}

# Makes plots:
plotter_beta <- function(dist_methods, study_pseq, pcoa_list){
  for( i in dist_methods ){
    p <- NULL
    
    p <- plot_ordination(study_pseq, pcoa_list[[i]], 
                         color = "group") + 
          geom_point(size = 2) +
          stat_ellipse(aes(group = group), linetype = 2) +
          theme_bw() +
          theme(plot.title = element_text(face = 'bold', size = 16),
                axis.title = element_text(size = 14),
                legend.title = element_text(size = 14),
                axis.text = element_text(size = 12),
                legend.text = element_text(size = 12))
    
    p <- p + ggtitle(i) 
    p <- p + scale_colour_brewer(type="qual", palette="Set2")
    
    plist[[i]] <- p
  }
  
  return(plist)
}

## Calculations -----
# gets distance methods of interest:
dist_methods <- unlist(distanceMethodList)[c(1, 2, 8, 10)]
dist_methods

# to do this individually for each study:
dlist <- vector("list", length(dist_methods)) # distance matrix
names(dlist) <- dist_methods
pcoa_list <- dlist # PCoA
plist <- dlist # plots
# example: with Feng dataset
feng.beta <- distances(pseq.feng)
feng.beta.adonis <- adonis_calculator(feng.beta$dlist, pseq.feng)
feng.beta.plot <- plotter_beta(dist_methods, pseq.feng, feng.beta$pcoa_list)

feng.beta.adonis
ggarrange(plotlist = feng.beta.plot,
          common.legend = TRUE, legend = "right")
# and so on...

# focusing on the merged dataset:
dlist <- vector("list", length(dist_methods)) # distances
names(dlist) <- dist_methods
pcoa_list <- dlist # PCoA
plist <- dlist # plots

# Calculate distances and ordination:
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- distance(adjusted_phyloseq, method = i)
  dlist[[i]] <- iDist
  
  # Calculate ordination
  iPCoA  <- ordinate(adjusted_phyloseq, "PCoA", distance = iDist)
  pcoa_list[[i]] <- iPCoA
}

# permanova:
results.group <- lapply(names(dlist), 
                      function(x) {
                        z <- adonis(dlist[[x]] ~ group, 
                                    data = 
                                      data.frame(sample_data(adjusted_phyloseq)))
                        return(as.data.frame(z$aov.tab))
                      })
names(results.group) <- names(dlist)
results.group

# plots:
plist <- plotter_beta(dist_methods, adjusted_phyloseq, pcoa_list)
ggarrange(plotlist = plist,
          common.legend = TRUE, legend = "right")

ggarrange(plist[[1]] + ggtitle('a'),
          plist[[2]] + ggtitle('b'),
          plist[[3]] + ggtitle('c'),
          plist[[4]] + ggtitle('d'),
          common.legend = TRUE, legend = "bottom",
          align = 'hv')

# Saving data ------
tab.save <- adjusted_otu_table %>%
  apply(1, function(x) x/1E8)
  # returns to relative abundance values
  # and forces them to be between 0-1

# go back to relative abundances and force to 0-1 values:
adjusted_phyloseq <- transform_sample_counts(adjusted_phyloseq, 
                                             function(x) x/1E8)
# save phyloseq object:
saveRDS(adjusted_phyloseq, 'DATA/adjusted_phyloseq.RDS')
# save relative abundances table:
write.table(tab.save, 'DATA/adjusted_abundances.txt',
            quote = FALSE, sep = '\t', row.names = TRUE, col.names = TRUE)