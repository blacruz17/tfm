library(caret)
library(pROC)
library(ggpubr)
library(tidyverse)
library(viridis)

# Loading data ----
class.RF <- read.table('DATA/randomForest/classRF.tsv',
                       header = TRUE, sep = '\t',# skip = 1,
                       row.names = 1)
## classes for binary and multiclass predictors:
class.RF$class.binary <- gsub('O|NO', '', class.RF$class)
class.binary <- factor(class.RF$class.binary)
class.multi <- factor(class.RF$class)

## relative abundances table:
data.RF <- read.table('DATA/randomForest/dataRF.tsv',
                      header = TRUE, sep = '\t',# skip = 1,
                      row.names = 1)
# keeps long names:
long.names <- colnames(data.RF)
# shortens names:
colnames(data.RF) <- gsub('k_.*_s_', 's_', colnames(data.RF))

# Preprocessing -----
# same as in the tuning script, removes correlated variables:
descrCorr <- cor(data.RF)
highCorr <- findCorrelation(descrCorr, 0.90)
data.RF.uncor <- data.RF[ , -highCorr]

# Classifier ------
## Train/Test split ----
set.seed(505)
inTrain <- createDataPartition(class.multi, 
                               p = 3/4, list = FALSE)

trainDescr <- data.RF.uncor[inTrain, ]
testDescr <- data.RF.uncor[-inTrain, ]
trainClass.bi <- class.binary[inTrain]
testClass.bi <- class.binary[-inTrain]

# sanity check:
dim(trainDescr)[1] == length(trainClass.bi)
dim(testDescr)[1] == length(testClass.bi)

## Model training -----
# In this case we will only define one value
# for each hyperparameter:
mtry.grid <- expand.grid(.mtry = floor(sqrt(ncol(trainDescr))))
ntrees <- 1000

set.seed(505)
seeds <- vector(mode = "list", length = 51)
for(i in 1:50) seeds[[i]] <- sample.int(1000, 22)
## For the last model:
seeds[[51]] <- sample.int(1000, 1)


fitControl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 5,
  savePredictions = TRUE,
  classProbs = TRUE,
  summaryFunction = multiClassSummary,
  verboseIter = TRUE,
  search = "grid",
  seeds = seeds
  )

# train the model:
fit <- train(trainDescr, trainClass.bi,
               method = 'rf',
               tuneGrid = mtry.grid,
               trControl = fitControl,
               ntree = ntrees)

# These two objects will be necessary to validate our model
# in accompanying scripts:
saveRDS(fit, 'DATA/randomForest/bestmodel.rds')
saveRDS(trainDescr, 'DATA/randomForest/trainDescr_bestmodel.rds')

# ROC Curves -------------------------
predictions <- predict(fit, testDescr, type = 'prob')
roc_curve <- roc(testClass.bi, predictions$MU)

# compare with the random predictor:
set.seed(505)
random_curve <- roc(sample(testClass.bi), predictions$MU)
roc.list.shuf <- list('Model' = roc_curve, 'Shuffled' = random_curve)
p.shuf.ci <- roc.plot(roc.list.shuf, 
                      shuffle = TRUE, title.in = '', p.ci = TRUE)

p.shuf.ci

# Variable importance -----
imp.bi <- varImp(fit, scale = FALSE)

imp.bi.df <- imp.bi$importance %>% 
              as.data.frame() %>%
              tibble::rownames_to_column() %>%
              rename(feature = rowname) %>%
              arrange(desc(Overall))

## Abundance testing -------
clr.abundances <- data.frame(t(data.RF))
clr.subset <- clr.abundances[imp.bi.df[1:20, 'feature'], ]

clr.mu <- clr.subset %>% 
  select(all_of(mu.samples)) %>%
  mutate(mean = rowMeans(.),
         type = 'MU') %>%
  rownames_to_column('feature')

clr.mh <- clr.subset %>% 
  select(all_of(mh.samples))%>%
  mutate(mean = rowMeans(.),
         type = 'MH') %>%
  rownames_to_column('feature')

# join both datasets, some aesthetic fixes, 
# wilcoxon test and p-value correction:
both.clr <- left_join(clr.mu %>% 
                        select(feature, mean) %>% 
                        rename(mean.mu = mean), 
                      clr.mh %>% 
                        select(feature, mean) %>% 
                        rename(mean.mh = mean),
                      by = 'feature') %>%
  select(mean.mu, mean.mh, feature) %>%
  mutate(feature = gsub('k_.*_s_', '', feature)) %>%
  column_to_rownames('feature') %>%
  mutate(p = apply(clr.subset, 1, 
                   function(x) wilcox.test(x[mu.samples], 
                                           x[mh.samples])$p.value)) %>%
  rstatix::adjust_pvalue(method = 'BH') %>%
  rstatix::add_significance()

# keep significant results and add a column
# indicating the group with the greater mean:
greater.clr <- both.clr %>% 
  filter(p.adj.signif != 'ns') %>%
  mutate(greater = if_else(
    mean.mu - mean.mh > 0, 'MU', 'MH')) 

## Heatmap ------
library(gplots)
library(Heatplus)
library(RColorBrewer)

adjusted_phyloseq <- readRDS('DATA/adjusted_phyloseq.rds')
clr_phyloseq <- microbiome::transform(adjusted_phyloseq, 'clr')

# prepares the input:
otus <- data.frame(clr_phyloseq@otu_table)
rownames(otus) <- tolower(gsub('k_.*s__', 's_', rownames(otus)))

markers <- rownames(greater.clr[greater.clr$greater == 'MU', ])
otus_mu <- otus[markers, ]

selection <- imp.bi.df[1:20, 'feature']
otus_imp <- otus[selection, ]

t_otus_imp_id <- data.frame(t(otus_imp)) %>% 
                  tibble::rownames_to_column(('ID'))
metadata <- data.frame(adjusted_phyloseq@sam_data) %>%
                  tibble::rownames_to_column('ID')

otus_and_metadata <- left_join(t_otus_imp_id, 
                               metadata, 
                               by = 'ID') %>% 
                      arrange(desc(metab))

otus_by_metab <- otus_and_metadata %>% 
                  column_to_rownames('ID') %>%
                  select(-colnames(adjusted_phyloseq@sam_data))

# now we will introduce some changes to correctly 
# format feature names in the plot:
df.to.melt <- data.frame(t(otus_by_metab))

df.to.melt$feature <- rownames(df.to.melt)
df.to.melt$feature <- gsub('^s_', '', df.to.melt$feature)
df.to.melt$feature <- gsub('_', ' ', df.to.melt$feature)
df.to.melt$feature <- Hmisc::capitalize(df.to.melt$feature)

# We want features (columns) to be ordered by their importance:
row_order <- imp.bi.df[1:20, 'feature']
row_order <- gsub('^s_', '', row_order) # repeat aesthetic name fixes
row_order <- gsub('_', ' ', row_order)
row_order <- Hmisc::capitalize(row_order)

df.to.melt$feature <- factor(df.to.melt$feature, levels = row_order)

# melt the dataframe:
df.melted <- melt(df.to.melt, id = c('feature'))
df.melted$variable <- factor(df.melted$variable, 
                            levels = colnames(df.to.melt))

## Create the heatmap:
hm <- df.melted %>%  
        ggplot(aes(x = feature, y = variable, fill = value)) +
          geom_tile() + 
          scale_fill_distiller(name = "",
                               palette = "PuBuGn", direction = 1,
                               na.value = "transparent") +
          scale_colour_manual(values = bluepalette) +
          scale_x_discrete(breaks = unique(df.melted$feature),
                           labels = unique(df.melted$feature)) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, 
                                           vjust = 0.5, hjust = 1)) +
          theme(legend.position = "bottom", legend.direction = "horizontal",
                legend.title = element_text(size = 15), 
                legend.key.size = unit(1,"cm"),
                legend.text = element_text(size = 10)) +
          guides(fill = guide_colorbar(title.position = "top", 
                                       title.hjust = 0.5))

hm

# to check this is doing what it should 
# (we want MH patients on top and MU in the bottom)
# gets labels from the plot:
hm_labels_patients <- ggplot_build(hm)$layout$panel_params[[1]]$y$get_labels()
# checks this vector is equal to the colnames of the sorted df:
all(colnames(df.to.melt)[1:356] == hm_labels_patients)
# TRUE -- great!

# clean heatmap:
hm.clean <- hm +
            theme(axis.title.y = element_blank(), 
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x = element_text(size = 12, face = 'italic'),
                  legend.position = "right",
                  legend.direction = "vertical")

# Now we will create a feature importance barplot 
# specifically for this plot.
# First, we need the feature importance data:
a <- imp.bi$importance %>% 
      as.data.frame() %>%
      tibble::rownames_to_column() %>%
      arrange(desc(Overall)) %>%
      filter(rowname %in% imp.bi.df[1:20, 'feature']) %>%
      mutate(rowname = gsub('k_.*_s_', '', rowname)) %>%
      mutate(rowname = forcats::fct_inorder(rowname)) %>%
      rename(feature = rowname) %>%
      mutate(type = "binary")

# we qill colour feature importance bars according to MH/MU status:
pal <- rep(viridis(1), 20)
idx <- which(a$feature %in% 
             rownames(greater.clr %>% filter(greater == 'MU')))
pal[idx] <- viridis(2)[2]

# create the barplot:
p3 <- a %>%
      arrange(desc(Overall)) %>%
      ggplot(aes(x = feature, 
                 y = Overall))+
      geom_col(width = 0.5, fill = pal_2) +
      labs(y = 'Importance', x = 'Features') +
      theme_minimal() +
      ggtitle('a') +
      theme(plot.title = element_text(face = 'bold', size = 16),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size = 14),
            axis.text.y = element_text(size = 12))

# get grob widths right to combine the heatmap + barplot:
gA <- ggplotGrob(p3)
gB <- ggplotGrob(hm.clean +
                   theme(axis.text.x = element_text(face = 'italic',
                                                    size = 12)))

gB <- ggplotGrob(hm.clean)

maxWidth <- grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth)
gB$widths[2:5] <- as.list(maxWidth)

# Next, a sidebar to tell if patients are MU or MH:
metab_df_to_melt <- otus_and_metadata[ ,c('ID', 'metab')]
row_order <- hm_labels_patients

# convert to factors to get orders right:
metab_df_to_melt$ID <- factor(metab_df_to_melt$ID,
                              levels = row_order)
metab_df_to_melt$metab <- factor(metab_df_to_melt$metab, 
                                 levels = c('healthy', 'unhealthy'))

# melt dataframe - and convert to factors again:
metab_df_melted <- melt(metab_df_to_melt, id = c('metab'))
metab_df_melted$variable <- metab_df_melted$value
metab_df_melted$variable <- factor(metab_df_melted$variable, 
                                   levels = metab_df_to_melt$ID)
# convert to numeric data:
metab_df_melted$value <- 1
metab_df_melted$value[metab_df_melted$metab == 'healthy'] <- 0

# create the sidebar:
bar <- metab_df_melted %>% 
          ggplot(aes(x = 1, y = variable, fill = value)) +
          scale_fill_viridis() +
          geom_tile() +
          theme_minimal() +
          theme(axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank(),
                legend.position = 'none')

## Put all the plots together:
# but first, a few more fixes to the heatmap:
hm.clean <- hm.clean + 
            ggtitle('b') +
            theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(), 
                  plot.title = element_blank())

bar_hm <- ggarrange(bar, hm.clean, widths = c(0.05, 1))

panel_bar_hm <- bar_hm + 
                ggtitle('b') + 
                theme(plot.title = element_text(face = 'bold', 
                                                size = 16, 
                                                hjust = 0.06))
# and, finally, get our plot:
gridExtra::grid.arrange(gA, panel_bar_hm,
                        ncol = 1, heights = c(20, 70))

## Build the ROC curve comparing with validation data: ----
# To run this last part, it is necessary to obtain the final model,
# export it, use it for validation in the celiac disease and CRC
# cohorts and then import the resulting ROC curves.

# helper function for plotting:
roc.plot <- function(roc.list, title.in, p.ci = FALSE,
                     shuffle = TRUE){
  # generates palette
  if (shuffle == TRUE){
    lengthroc <- length(roc.list)/2
    pal <- rep(viridis::viridis(lengthroc), 2)
    aes_ggroc <- c("linetype", "colour")
    line_types <- c(rep("solid", lengthroc), rep("twodash", lengthroc))
  }
  else {
    pal <- viridis::viridis(length(roc.list))
    aes_ggroc <- "colour"
    line_types <- "solid"
  }
  
  # extracts auc
  roc.list %>% 
    map(~tibble(AUC = .x$auc)) %>% 
    bind_rows(.id = "name") -> data.auc
  
  # generates labels
  data.auc %>% 
    mutate(label_long=paste0(name,", AUC = ",paste(round(AUC,2))),
           label_AUC=paste0("AUC = ",paste(round(AUC,2)))) -> data.labels
  
  names(roc.list) <- data.labels$label_long
  
  p <- ggroc(c(roc.list),
             size = 1.2, legacy.axes = TRUE, aes = aes_ggroc) + 
    geom_line(size = 1.2) +
    labs(x = "False Positive Rate", 
         y = "True Positive Rate",
         title = title.in) +
    scale_color_manual(values = pal) +
    scale_linetype_manual(values = line_types) +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                 color = "black", linetype = "dotted") +
    theme_bw() +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.title = element_blank(),
          legend.text = element_text(size = 12)) +
    coord_equal()
  
  if (p.ci == TRUE){
    ## plots + confidence intervals
    ci.list <- lapply(roc.list, ci.se, specificities = seq(0, 1, l = 25))
    
    dat.ci.list <- lapply(ci.list, function(ciobj) 
      data.frame(x = as.numeric(rownames(ciobj)),
                 upper =  1- ciobj[, 1],
                 lower = 1- ciobj[, 3]))
    
    for(i in 1:length(roc.list)) {
      p <- p + geom_ribbon(
        data = 1 - dat.ci.list[[i]],
        aes(x = x, ymin = lower, ymax = upper),
        fill = pal[i],
        alpha = 0.2,
        inherit.aes = FALSE) 
    }}
  return(p)
}

# Validate the data and import the ROC curves:
# note: the roc_celiac file is not available
roc_celiac <- readRDS('DATA/randomForest/roc_celiac.rds')
roc_feng <- readRDS('DATA/randomForest/roc_feng.rds')

rocs <- list('Holdout' = roc_curve, 
             'Celiac disease' = roc_celiac, 
             'CRC' = roc_feng)

rocs_plot <- roc.plot(rocs, '', 
                      p.ci = FALSE, shuffle = FALSE)
# get our nice figure:
ggarrange(p.shuf.ci + 
            ggtitle('a') + 
            theme(plot.title = element_text(size = 16, face = 'bold')),
          rocs_plot + 
            ggtitle('b') + 
            theme(plot.title = element_text(size = 16, face = 'bold')),
          nrow = 2, align = 'hv')