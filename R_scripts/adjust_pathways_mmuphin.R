library(tidyverse)
library(janitor)
library(MMUPHin)

library(phyloseq)
library(scater)
library(vegan)
library(mia)

library(ggpubr)

metaphlanToPhyloseq <- source('DATA/metaphlanToPhyloseq.R')$value

# Pathway abundances ----
# load abundances table obtained with HUMAnN:
abundances <- read_tsv('./DATA/pathabundance_unstratified_head-fix.tsv')

abundances <- as.data.frame(abundances) %>%
                    column_to_rownames(var = 'MetaCyc_pathway') %>%
                    select(sort(names(.)))

# load metadata file:
metadata <- read_tsv('./DATA/metadata_runacc.tsv') %>%
            arrange(sample_id)
metadata <- column_to_rownames(metadata, var = 'sample_id')

pseq_abundances <- metaphlanToPhyloseq(abundances,
                                       metadat = metadata)

pca_plot <- function(pseq) {
      pseq_ord <- ordinate(pseq, 'PCA')
      p <- plot_ordination(pseq, pseq_ord, 
                                type="samples", color="study_name") + 
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
p_before_abundance <- pca_plot(pseq_abundances)

# adjust for batch effect:
fit_adjust_batch <- adjust_batch(feature_abd = abundances,
                                 batch = "study_name",
                                 covariates = "group",
                                 data = metadata,
                                 control = list(verbose = FALSE))

to_write <- as.data.frame(fit_adjust_batch$feature_abd_adj) %>%
              rownames_to_column('MetaCyc_pathway')

write.table(to_write,
            'DATA/adjusted_pathabundance.tsv',
            quote = FALSE, row.names = FALSE, sep = '\t')

# to see the PCA plot:
abundances_table <- fit_adjust_batch$feature_abd_adj
pseq_abundances_after <- metaphlanToPhyloseq(abundances_table,
                                      metadat = metadata)

p_after_abundance <- pca_plot(pseq_abundances_after)

# Pathway coverages ----
# load coverages table obtained with HUMAnN:
coverages <- read_tsv('DATApathcoverage_unstratified_head-fix.tsv')

coverages <- as.data.frame(coverages) %>%
  column_to_rownames(var = 'MetaCyc_pathway') %>%
  select(sort(names(.)))

pseq_coverages <- metaphlanToPhyloseq(coverages,
                                       metadat = metadata)

p_before_coverage <- pca_plot(pseq_coverages)

# adjust for batch effect:
fit_adjust_batch <- adjust_batch(feature_abd = coverages,
                                 batch = "study_name",
                                 covariates = "group",
                                 data = metadata,
                                 control = list(verbose = FALSE))


coverages_table <- fit_adjust_batch$feature_abd_adj
pseq_coverages_after <- metaphlanToPhyloseq(coverages_table,
                                       metadat = metadata)

p_after_coverage <- pca_plot(pseq_coverages_after)

to_write <- as.data.frame(fit_adjust_batch$feature_abd_adj) %>%
  rownames_to_column('MetaCyc_pathway')

write.table(to_write,
            'DATA/adjusted_pathcoverage.tsv',
            quote = FALSE, row.names = FALSE, sep = '\t')

# we can take a look at the batch effect correction now:
ggarrange(p_before_abundance + ggtitle('a'), 
          p_after_abundance + ggtitle(' '),
          p_before_coverage + ggtitle('b'), 
          p_after_coverage + ggtitle(' '),
          common.legend = TRUE, legend = 'right',
          align = 'hv')

# Adonis(PERMANOVA) ----
## Abundances ----- 
D_before <- vegdist(t(abundances))
D_after <- vegdist(t(abundances_table))

set.seed(1)
fit_adonis_before <- adonis(D_before ~ study_name, 
                            data = metadata)
fit_adonis_after <- adonis(D_after ~ study_name, 
                           data = metadata)
print(fit_adonis_before)
print(fit_adonis_after)

## Coverages ----
D_before <- vegdist(t(coverages))
D_after <- vegdist(t(coverages_table))

set.seed(1)
fit_adonis_before <- adonis(D_before ~ study_name, 
                            data = metadata)
fit_adonis_after <- adonis(D_after ~ study_name, 
                           data = metadata)
print(fit_adonis_before)
print(fit_adonis_after)
