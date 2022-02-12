library(caret)
library(MLmetrics)
library(pROC)
library(ggpubr)
library(tidyverse)

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

###### NOTE ######
## To work with the model with metadata, one can use this same
## script, but with three small changes:
# 1 - load the file containing metadata information
data.RF <- read.table('DATA/randomForest/dataRF_metadata.tsv',
                      header = TRUE, sep = '\t',# skip = 1,
                      row.names = 1)
# 2 - remove a sample with age == 0:
age.0 <- rownames(metadata.RF[metadata.RF$age == 0, ])
metadata.RF <- metadata.RF[metadata.RF$age != 0, ]
class.RF <- subset(class.RF, !rownames(class.RF) %in% age.0)

# 3 - before training the multiclass model, we would also 
#     have to remove the BMI from the training set:
trainDescr <- data.RF.uncor[inTrain, names(metadata.RF.uncor) != "bmi"]
testDescr <- data.RF.uncor[-inTrain, names(metadata.RF.uncor) != "bmi"]
# the rest of the script is exactly the same
##################

# Preprocessing -----
# checks for zero and near-zero variance variables:
data.nzv <- nearZeroVar(data.RF, saveMetrics = TRUE)
sum(data.nzv$zeroVar) ## 0
sum(data.nzv$nzv) ## 0
# there are none, no need to remove any then

# removes correlated variables:
descrCorr <- cor(data.RF)
highCorr <- findCorrelation(descrCorr, 0.90)
length(highCorr) # 117 variables
data.RF.uncor <- data.RF[ , -highCorr]

dim(data.RF.uncor) # 356  543

# Multiclass classifier ------
## Train/Test split ----
# creates split based on the class labels
set.seed(505)
inTrain <- createDataPartition(class.multi, 
                               p = 3/4, list = FALSE)
# selects the corresponding samples:
trainDescr <- data.RF.uncor[inTrain, ]
testDescr <- data.RF.uncor[-inTrain, ]
trainClass <- class.multi[inTrain]
testClass <- class.multi[-inTrain]

# these should be TRUE:
dim(trainDescr)[1] == length(trainClass)
dim(testDescr)[1] == length(testClass)
# OK!

## Model training -----
# creates grid of mtry values:
mtry.grid <- expand.grid(.mtry = floor(seq(sqrt(ncol(trainDescr)), 
                                           ncol(trainDescr), 
                                           length.out = 5)))
# creates sequence of ntree values:
ntrees <- seq(200, 2000, 400)

# sets seeds for reproducibility:
set.seed(505)
seeds <- vector(mode = "list", length = 51)
for(i in 1:50) seeds[[i]] <- sample.int(1000, 22)
## For the last model:
seeds[[51]] <- sample.int(1000, 1)

### Without subsampling --------
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

originals <- list()

# train with different ntree parameters:
for (ntree in ntrees){print(ntree)
  set.seed(505)
  fit <- train(trainDescr, trainClass,
               method = 'rf',
               tuneGrid = mtry.grid,
               trControl = fitControl,
               ntree = ntree)
  key <- toString(ntree)
  originals[[key]] <- fit
}

### upsampling --------------
fitControl$sampling <- "up"
upsampled <- list()
for (ntree in ntrees){
  set.seed(505)
  fit <- train(trainDescr, trainClass,
               method = 'rf',
               tuneGrid = mtry.grid,
               trControl = fitControl,
               ntree = ntree)
  key <- toString(ntree)
  upsampled[[key]] <- fit
}

### downsampling --------------
fitControl$sampling <- "down"
downsampled <- list()
for (ntree in ntrees){
  set.seed(505)
  fit <- train(trainDescr, trainClass,
               method = 'rf',
               tuneGrid = mtry.grid,
               trControl = fitControl,
               ntree = ntree)
  key <- toString(ntree)
  downsampled[[key]] <- fit
}

# Binary classifier --------
## Data preparation
## Train/Test split ----
trainClass.bi <- class.binary[inTrain]
testClass.bi <- class.binary[-inTrain]

# sanity check:
dim(trainDescr)[1] == length(trainClass.bi)
dim(testDescr)[1] == length(testClass.bi)
# OK!

## Model training ------
### Without subsampling --------
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

originals.bi <- list()

for (ntree in ntrees){print(ntree)
  set.seed(505)
  fit <- train(trainDescr, trainClass.bi,
               method = 'rf',
               tuneGrid = mtry.grid,
               trControl = fitControl,
               ntree = ntree)
  key <- toString(ntree)
  originals.bi[[key]] <- fit
}

### upsampling --------------
fitControl$sampling <- "up"
upsampled.bi <- list()
for (ntree in ntrees){
  set.seed(505)
  fit <- train(trainDescr, trainClass.bi,
               method = 'rf',
               tuneGrid = mtry.grid,
               trControl = fitControl,
               ntree = ntree)
  key <- toString(ntree)
  upsampled.bi[[key]] <- fit
}

### downsampling --------------
fitControl$sampling <- "down"
downsampled.bi <- list()
for (ntree in ntrees){
  set.seed(505)
  fit <- train(trainDescr, trainClass.bi,
               method = 'rf',
               tuneGrid = mtry.grid,
               trControl = fitControl,
               ntree = ntree)
  key <- toString(ntree)
  downsampled.bi[[key]] <- fit
}

# Comparing hyperparameters -------
# helper function for plotting:
hyper.plot <- function(modellist, title.in){
  
  a <- lapply(modellist, function(x) {
    acc <- x$results$Accuracy
  })
  
  a.df <- data.frame(a)
  rownames(a.df) <- mtry.grid$.mtry
  colnames(a.df) <- gsub('X', '', colnames(a.df))
  a.df$ntree <- colnames(a.df)
  
  data_ggp <- data.frame(x = as.numeric(rownames(a.df)),                            # Reshape data frame
                         y = c(a.df$`200`, a.df$`600`, 
                               a.df$`1000`, a.df$`1400`,
                               a.df$`1800`),
                         ntrees = c(rep("200", nrow(a.df)),
                                   rep("600", nrow(a.df)),
                                   rep("1000", nrow(a.df)),
                                   rep("1400", nrow(a.df)),
                                   rep("1800", nrow(a.df))))
  # this line makes sure the order of ntrees in the legend is right:
  data_ggp$ntrees <- factor(data_ggp$ntrees, 
                           levels = c('200', '600', '1000', '1400', '1800'))
  
  best_ggp_x <- data_ggp[which.max(data_ggp$y), 'x']
  best_ggp_y <- max(data_ggp$y)
  
  ggp <- ggplot(data_ggp, aes(x, y, col = ntrees)) +
    geom_line(size = 1) +
    geom_point(aes(best_ggp_x, best_ggp_y),
               size = 3, shape = 17, color = 'black') +
    labs(y = 'Accuracy', x = 'mtry') +
    scale_x_continuous(name="mtry",
                       breaks = mtry.grid$.mtry,
                       labels = mtry.grid$.mtry) +
    ggtitle(title.in) +
    theme_bw()
  return(ggp)           
}

# helper function to get the best model:
get_best_model <- function(modellist){
  acclist <- list()
  for (i in names(modellist)){
    key <- i
    acc <- modellist[[i]]$results$Accuracy
    acclist[[i]] <- max(acc)
  }
  idx_best <- names(which.max(acclist))
  return(modellist[idx_best])
}

# helper function for plotting roc curves:
roc.plot <- function(roc.list, title.in, p.ci = FALSE,
                     shuffle = TRUE){
  # generates palette
  # if we want to add a shuffled version, we want
  # original and random to have the same colour,
  # originals to have solid lines and randoms to have dashed lines
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
  # generates plot:
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
    theme(legend.title = element_blank()) +
    coord_equal()
  # to add confidence intervals:
  if (p.ci == TRUE){
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

# helper function to build multiclass ROC curves:
roccer.multi <- function(model, data, tags){
  result.predicted.prob <- predict(model, data, type="prob") 
  result.roc <- multiclass.roc(tags,
                               result.predicted.prob)
  return(result.roc)
}

# ROC Curves -------------------------
# helper function to calculate roc curves:
roccer <- function(model, data, tags){
  result.predicted.prob <- predict(model, data, type="prob") # Prediction
  result.roc <- roc(tags, result.predicted.prob$MU)
  return(result.roc)
}

# shuffled versions of the class labels:
set.seed(505)
shuf.bi <- sample(testClass.bi)
set.seed(505)
shuf.multi <- sample(testClass)

## Binary model ------
# takes a look at the hyperparameter tuning process:
og.bi.p <- hyper.plot(originals.bi, 'Binary model (original)')
up.bi.p <- hyper.plot(upsampled.bi, 'Binary model (upsampled)')
down.bi.p <- hyper.plot(downsampled.bi, 'Binary model (downsampled)')

best_binary_models <- c(get_best_model(originals.bi),
                        get_best_model(upsampled.bi),
                        get_best_model(downsampled.bi))

names(best_binary_models) <- c("Original", "Upsampled", "Downsampled")

# calculates roc curves and generates plots:
roc.list.bi <- lapply(best_binary_models, 
                      function(x) roccer(x, testDescr, testClass.bi))
roc.list.bi.s <- lapply(best_binary_models, 
                             function(x) roccer(x, testDescr, shuf.bi))

my.roc <- c(roc.list.bi, roc.list.bi.s)
names(my.roc)[4:6] <- paste(names(my.roc)[4:6], '+ shuffle')

## first panel of ROC curves figure:
roc_a <- roc.plot(my.roc, 'a',
         p.ci = FALSE, shuffle = TRUE) + 
          theme(plot.title = element_text(face="bold", size = 16),
                axis.title = element_text(size = 14),
                axis.text = element_text(size = 12),
                legend.text = element_text(size = 12))
roc_a


## Multiclass model ----
## Hyperparameters:
og.p <- hyper.plot(originals, 'Multiclass model (original)')
up.p <- hyper.plot(upsampled, 'Multiclass model (upsampled)')
down.p <- hyper.plot(downsampled, 'Multiclass model (downsampled)')

ggarrange(og.p, up.p, down.p,
          common.legend = TRUE,
          legend = 'right')

# get them all together for the figure:
ggarrange(og.bi.p, up.bi.p, down.bi.p,
          og.p, up.p, down.p,
          common.legend = TRUE,
          legend = 'right')

# this could look better...
# aesthetic fixes for the figure:
og.bi.p <- hyper.plot(originals.bi, 'a') + 
            theme(axis.title.x = element_blank(),
                  plot.title = element_text(face = "bold", size = 16),
                  axis.title.y = element_text(size = 14),
                  axis.text = element_text(size = 12),
                  legend.text = element_text(size = 12),
                  legend.title = element_text(size = 14)) +
            scale_y_continuous(labels = scales::label_number(accuracy = 0.001)) +
            labs(color = 'ntree')
og.bi.p
up.bi.p <- hyper.plot(upsampled.bi, 'b') + 
            theme(axis.title.x = element_blank(),
                  plot.title = element_text(face = "bold", size = 16),
                  axis.title.y = element_text(size = 14),
                  axis.text = element_text(size = 12))
down.bi.p <- hyper.plot(downsampled.bi, 'c') +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.001)) +
              theme(plot.title = element_text(face = "bold", size = 16),
                    axis.title.y = element_text(size = 14),
                    axis.title.x = element_text(size = 14),
                    axis.text = element_text(size = 12))

og.p <- hyper.plot(originals, '') + 
          theme(axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text = element_text(size = 12))
up.p <- hyper.plot(upsampled, '') + 
          theme(axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text = element_text(size = 12))
down.p <- hyper.plot(downsampled, '') +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.001)) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 12))

# This is a bit prettier now:
ggarrange(og.bi.p, og.p, 
          up.bi.p, up.p, 
          down.bi.p, down.p,
          ncol = 2, nrow = 3,
          common.legend = TRUE,
          legend = 'right', 
          align = 'hv')

# get the best models:
best_multiclass_models <- c(get_best_model(originals),
                            get_best_model(upsampled),
                            get_best_model(downsampled))

names(best_multiclass_models) <- c("Original", "Upsampled", "Downsampled")

## ROC curves:
roc.list.multi <- lapply(best_multiclass_models, 
                         function(x) roccer.multi(x, testDescr, testClass))
roc.list.multi.s <- lapply(best_multiclass_models,
                           function(x) roccer.multi(x, testDescr, shuf.multi))

## in multiclass models, two ROC curves are created 
## for each class combination, using each class as control. 
## We will keep only the first one:
rs <- roc.list.multi$Downsampled[['rocs']]
roc.multi <- list()
for (i in 1:length(rs)){
  key <- names(rs)[i]
  roc.multi[[key]] <- rs[[i]][[1]]
}

## get AUCs:
for (i in 1:length(roc.multi)){
  roc.multi[[i]]$auc <- auc(roc.multi[[i]])
}

# panel b for ROCs figure:
roc_b <- roc.plot(roc.multi, 'b',
         p.ci = FALSE, shuffle = FALSE) + 
          theme(plot.title = element_text(face="bold", size = 16),
                axis.title = element_text(size = 14),
                axis.text = element_text(size = 12),
                legend.text = element_text(size = 12))
roc_b

## on the shuffled version:
rs.s <- c(roc.list.multi$Original[['rocs']], 
          roc.list.multi.s$Original[['rocs']])
roc.multi.s <- list()

for (i in 1:length(rs.s)){
  key <- names(rs.s)[i]
  if (i > 6){key <- paste(key, '+ shuffle')}
  roc.multi.s[[key]] <- rs.s[[i]][[1]]
}

for (i in 1:length(roc.multi.s)){
  roc.multi.s[[i]]$auc <- auc(roc.multi.s[[i]])
}

# panel c for the ROCs figure:
roc_c <- roc.plot(roc.multi.s[c(1, 6, 7, 12)], 'c',
         p.ci = FALSE, shuffle = TRUE) + 
          theme(plot.title = element_text(face="bold", size = 16),
                axis.title = element_text(size = 14),
                axis.text = element_text(size = 12),
                legend.text = element_text(size = 12))
roc_c

# get them all together:
ggarrange(roc_a, roc_b, roc_c, 
          nrow = 3, ncol = 1, 
          align = 'hv')
# great!

# Variable importance plots -----
# for the multiclass model:
imp.multi <- lapply(best_multiclass_models, function(x) varImp(x, scale = FALSE))
chosen.model <- 'Original'

imp.multi.og <- imp.multi[[chosen.model]]$importance %>% 
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  rename(feature = rowname) %>%
  arrange(desc(Overall))

features_to_plot <- 20
# We will add a lot of fies because we want each figure
# to look exactly as it should (italics for species and so on)
df <- imp.multi[[chosen.model]]$importance %>% 
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  arrange(Overall) %>%
  filter(rowname %in% imp.multi.og[1:features_to_plot, 'feature']) %>%
  mutate(rowname = gsub('k_.*_s_', '', rowname)) %>%
  mutate(rowname = gsub('^s_', '', rowname)) %>%
  mutate(rowname = gsub('_cag_', '*_CAG_', rowname, fixed = TRUE)) %>%
  mutate(rowname = gsub('_sp*_CAG', '*_sp._CAG', rowname, fixed = TRUE)) %>%
  mutate(rowname = gsub('_', ' ', rowname)) %>%
  mutate(rowname = Hmisc::capitalize(rowname)) %>%
  mutate(rowname = paste0('*', rowname)) %>%
  mutate(rowname = if_else(grepl('CAG', rowname),
                           rowname, 
                           paste0(rowname, '*'))) %>%
  mutate(rowname = gsub('*Age*', 'Age', rowname, fixed = TRUE)) %>%
  mutate(rowname = gsub('*Bmi*', 'BMI', rowname, fixed = TRUE)) %>%
  mutate(rowname = fct_inorder(rowname)) %>%
  rename(feature = rowname) %>%
  mutate(type = "multiclass")
# generate barplot:
p <- df %>%
  arrange(desc(Overall)) %>%
  ggplot(aes(x = feature, 
             y = Overall, 
             fill = type))+
  geom_col(width = 0.7, position = 'dodge') +
  scale_fill_manual(values = '#21908CFF') +
  coord_flip() +
  labs(y = 'Importance', x = 'Features',
       fill = 'Model') +
  theme_bw()  + 
  ggtitle('a') +
  theme(plot.title = element_text(face = 'bold', size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'none') +
  theme(axis.text.y = ggtext::element_markdown())

# add a bit of space on the right side:
p + 
  theme(plot.margin = margin(r = 10, unit = "pt"))