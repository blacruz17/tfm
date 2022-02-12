library(tidyverse)
library(RColorBrewer)
library(janitor)
library(phyloseq)
library(igraph)
library(rstatix)
library(ggpubr)
library(gt)

# Data preparation ----
## import phyloseq object and clr-transform:
adjusted_phyloseq <- readRDS('DATA/adjusted_phyloseq.RDS')
ps_clr <- microbiome::transform(adjusted_phyloseq, 'clr')

# generate distance matrix based on Euclidean distances:
clr_dist_matrix <- phyloseq::distance(ps_clr, method = "euclidean") 

# we also need some metadata:
metadata <- data.frame(ps_clr@sam_data)
metadata$sample_id <- rownames(metadata)

# now, we only want some relative abundances:
otus.to.add <- ps_clr@otu_table
rownames(otus.to.add) <- gsub('k__.*s__', '', rownames(otus.to.add))
otus.to.keep <- c('Clostridium_leptum',
                  'Gordonibacter_pamelaeae',
                  'Collinsella_intestinalis',
                  'Eggerthella_lenta')

otus.to.add <- data.frame(t(otus.to.add[otus.to.keep, ])) %>%
                  rownames_to_column('sample_id')

metadata <- left_join(metadata, otus.to.add, by = 'sample_id')

# let's add alpha diversity data:
## approximate to counts multiplying by 1 million:
adjusted_phyloseq_counts <- transform_sample_counts(adjusted_phyloseq, 
                                                    function(x) round(1E6 * x))
alpha <- estimate_richness(adjusted_phyloseq_counts) %>%
  rownames_to_column('sample_id') %>%
  select(where(~ !any(is.na(.x)))) # some alpha measures have NA values

metadata <- left_join(metadata, alpha, by = 'sample_id') 

rownames(metadata) <- metadata$sample_id

# Network generation ----
# helper function that we will need to add edge weights:
flatten_beta_div <- function(beta_div){
  mx <- as.matrix(beta_div)
  ut <- upper.tri(mx)
  rows <- rownames(mx)[row(mx)[ut]]
  cols <- rownames(mx)[col(mx)[ut]]
  cor <- mx[ut]
  flat <- data.frame(rows, cols, beta = cor)
  return(flat)
}

# this function will test different parameters for our network:
try_networks <- function(pseq, dist_matrix, max_distances, 
                         metadata){
  results <- list()
  for (d in max_distances){
    # creates network from a phyloseq object and a distance matrix:
    ig <- make_network(pseq, distance = dist_matrix, max.dist = d)
    # adds aitchison distances as edge weights:
    df <- as_data_frame(ig, what = c('edges', 'vertices', 'both'))
    clr <- flatten_beta_div(dist_matrix)
    df$edge <- apply(df, 1, function(x) paste(x[1], x[2], sep = '--'))
    clr$edge <- apply(clr, 1, 
                      function(x) paste(x[1], x[2], sep = '--'))
    clr_2 <- clr %>% filter(edge %in% df$edge)
    
    weights <- clr_2[, 'beta']
    E(ig)$weight <- weights
    # adds metadata as vertex attributes:
    metadata.ig <- metadata %>%
      filter(sample_id %in% V(ig)$name) 
    for (n in colnames(metadata.ig)[2:length(colnames(metadata.ig))]){
      ig <- set_vertex_attr(ig, n,
                            index = metadata.ig$sample_id,
                            value = metadata.ig[[n]])
    }
    # creates clusters:
    com <- cluster_fast_greedy(ig)
    # returns information to results list:
    key <- toString(d)
    results[[key]]$network <- ig
    results[[key]]$communities <- com
  }
  return(results)
}

seq_to_try <- c(110, seq(80, 50, -5))
# build different networks with the clr-transformed abundances
# and the Euclidean distance matrix (= Aitchison distances)
tries <- try_networks(ps_clr, clr_dist_matrix, seq_to_try, metadata)

# we could take a look at the resulting networks:
sizelist <- vector(length = length(seq_to_try))
names(sizelist) <- seq_to_try
for (t in names(tries)){sizelist[[t]] <- gsize(tries[[t]]$network)}

orderlist <- vector(length = length(seq_to_try))
names(orderlist) <- seq_to_try
for (t in names(tries)){orderlist[[t]] <- gorder(tries[[t]]$network)}

# and the communities within then:
commnumber <- vector(length = length(seq_to_try))
names(commnumber) <- seq_to_try
for (t in names(tries)){commnumber[[t]] <- length(sizes(tries[[t]]$communities))}

maxcomm <- vector(length = length(seq_to_try))
names(maxcomm) <- seq_to_try
for (t in names(tries)) {maxcomm[[t]] <- max(sizes(tries[[t]]$communities))}

# Summarize this and build a nice table:
df_tries <- data.frame(rbind(sizelist, orderlist, commnumber, maxcomm),
                       row.names = c('Size', 'Order', 
                                     'No. of communities', 
                                     'Largest community')) 

colnames(df_tries) <- c('Whole', seq(80, 50, -5))

gt_tries <- gt(df_tries %>% rownames_to_column()) %>% 
            tab_style(style = cell_text(weight = 'bold'), 
                      locations = cells_column_labels()) %>% 
            tab_style(style = cell_text(weight = 'bold'), 
                      locations = cells_stub())
# gtsave(gt_tries,
#        "plots/network_configuration_table.png",
#        zoom = 1)

# We have chosen 65 as cutoff. Let's set this then:
ig <- tries[['65']]$network
ig_comms <- tries[['65']]$communities

# Differential cluster analysis -----
# This dataframe will help us with our analyses:
metadata.ig <- metadata %>%
                filter(sample_id %in% V(ig)$name) %>%
                mutate(membership = ig_comms$membership) 

## sanity check - did we do this right?
all(rownames(filter(metadata.ig, membership == 3)) == 
      communities(ig_comms)[['3']])
# try with any number 1-5. should be TRUE

## Calculating differences between communities ----
# Chi-square test for metabolic health status:
chisq_metab <- chisq.test(xtabs( ~ metab + membership, metadata.ig))
## pairwise:
chisq_pairs <- pairwise_prop_test(xtabs( ~ metab + membership, metadata.ig),
                                  p.adjust.method = 'BH')
chisq_pairs %>% filter(p.adj.signif != 'ns')
## keep this p-value for later:
p_metab_12 <- as.numeric(chisq_pairs[1, 4])

## other variables:
chisq.test(xtabs( ~ sex + membership, metadata.ig)) # ns
chisq.test(xtabs( ~ obesity_status + membership, metadata.ig)) # ns
chisq.test(xtabs( ~ study_name + membership, metadata.ig))

# kruskal-wallis test for quantitative variables:
kruskal <- apply(select(metadata.ig, -c('sample_id', 'group',
                                        'membership', 'metab', 'sex',
                                        'obesity_status', 'study_name')), 2,
                  function(x) kruskal.test(x, metadata.ig$membership)$p.value) 
kruskal <- as.data.frame(kruskal) %>%
            adjust_pvalue(method = 'BH', 
                                   p.col = 'kruskal', output.col = 'adj.p') %>%
            add_significance(p.col = 'adj.p', output.col = 'p.signif')

# which are different?
signif <- rownames(kruskal[kruskal$p.signif != 'ns', ])

# pairwise tests (wilcoxon):
pairs <- apply(select(metadata.ig, signif), 2,
                function(x) pairwise.wilcox.test(x, metadata.ig$membership,
                                                  p.adjust.method = 'BH'))
# this is the same thing, a bit clumsier code but 
# it makes things more easily readable later
# we only do this for significantly different variables:
c_leptum <- metadata.ig %>% 
              pairwise_wilcox_test(Clostridium_leptum ~ membership)
g_pamelaeae <- metadata.ig %>% 
    pairwise_wilcox_test(Gordonibacter_pamelaeae ~ membership)

c_intestinalis <- metadata.ig %>% 
  pairwise_wilcox_test(Collinsella_intestinalis ~ membership)

e_lenta <- metadata.ig %>% 
  pairwise_wilcox_test(Eggerthella_lenta ~ membership)

chao1 <- metadata.ig %>% 
  pairwise_wilcox_test(Chao1 ~ membership)

fisher <- metadata.ig %>% 
  pairwise_wilcox_test(fisher ~ membership)

age <- metadata.ig %>% 
  pairwise_wilcox_test(age ~ membership)

# get all the results together:
all <- rbind(c_leptum, g_pamelaeae, e_lenta, c_intestinalis,
             chao1, fisher, age) %>%
        adjust_pvalue(method = 'BH')
        
table.clusters <- all %>% filter(p.adj.signif != 'ns')
table.clusters

## Obtaining boxplots for each biomarker ------
# keep only samples from communities 1 and 2
# and information for the 4 biomarkers:
markers_communities <- metadata.ig %>%
                        filter(membership == 1 | membership == 2) %>%
                        select(Clostridium_leptum, Gordonibacter_pamelaeae,
                               Eggerthella_lenta, Collinsella_intestinalis,
                               membership)
colnames(markers_communities) <- c('C. leptum', 'G. pamelaeae', 
                                   'E. lenta', 'C. intestinalis', 'membership')
# melt the data frame and plot:
markers_communities_long <- reshape2::melt(markers_communities, 
                                           id.vars = 'membership')
pal <- viridis::viridis(2, direction = -1)
panel_boxplots <- ggplot(markers_communities_long, 
                          aes(x = membership, y = value, 
                          group = as.factor(membership))) + 
                  geom_boxplot(aes(fill = as.factor(membership)), alpha = .7,
                                outlier.shape = NA) + 
                  scale_fill_manual(values = pal, name = 'Community') +
                  geom_jitter(width = 0.2, 
                              aes(colour = as.factor(membership)), size = 1.5) +
                  facet_wrap(~ factor(variable, 
                             levels = c('C. leptum', 'G. pamelaeae',
                                        'E. lenta', 'C. intestinalis')), ncol = 4) +
                  scale_color_manual(values = pal,
                                     name = 'Community') +
                  ggtitle('a') +
                  ylab('CLR-abundances') +
                  scale_x_discrete(name ="Community", 
                                   limits=c("1", "2")) +
                  theme_bw() +
                  theme(plot.title = element_text(size = 16, face = 'bold'),
                        strip.text = element_text(size = 12, face = 'italic'),
                        axis.text = element_text(size = 12),
                        axis.title = element_text(size = 14),
                        legend.position = 'none')
# to add p-values:
stat_boxplots <- all %>% filter(group1 == '1' & group2 == '2')
stat_boxplots <- stat_boxplots[1:4, ]
stat_boxplots$variable <- c('C. leptum', 'G. pamelaeae', 
                            'E. lenta', 'C. intestinalis')
stat_boxplots$.y. <- 'value'
stat_boxplots$y.position <- 11

panel_boxplots <- panel_boxplots + 
                    stat_pvalue_manual(stat_boxplots, 
                                       label = "p.adj.signif",
                                       hide.ns = T)

panel_boxplots

## Obtaining a summary table -----
comm_table <- all %>% 
                filter(group1 == '1' & group2 == '2',
                       .y. %in% otus.to.keep) %>% 
                select(.y., p.adj)
comm_table <- rbind(c('Metabolic_category', 
                      formatC(p_metab_12, digits = 3)), 
                    comm_table)

# great, now we need to add mean values:
means <- markers_communities %>%
          group_by(membership) %>%
          summarise_at(vars(colnames(markers_communities)[1:4]), list(mean))

comm_table$comm2 <- as.numeric(means[2, ])
comm_table$comm1 <- as.numeric(means[1, ])

# and now, the % of MU nodes:
MU_1 <- sum(metadata.ig[metadata.ig$membership == 1, 'metab'] == 'unhealthy')
MU_2 <- sum(metadata.ig[metadata.ig$membership == 2, 'metab'] == 'unhealthy')

perc_MU_1 <- round(100 * MU_1/sum(metadata.ig$membership == 1), 2)
perc_MU_2 <- round(100 * MU_1/sum(metadata.ig$membership == 2), 2)

comm_table[1, 'comm1'] <- perc_MU_1
comm_table[1, 'comm2'] <- perc_MU_2

# a few aesthetic fixes:
comm_table <- comm_table %>%
              mutate(comm1 = round(comm1, 2),
                     comm2 = round(comm2, 2)) %>%
              relocate(.y., comm1, comm2, p.adj) %>%
              rename('Feature' = .y.,
                     'Community 1' = comm1,
                     'Community 2' = comm2,
                     'p-value' = p.adj) %>% 
              mutate(Feature = c('MU nodes (%)',
                             colnames(markers_communities)[1:4]))

# and now we can get our nice table:
panel_comm_table <- ggtexttable(comm_table, 
                                rows = NULL) %>% 
                    table_cell_font(row = 3:6, column = 1, 
                                    face = 'bold.italic', size = 11) %>%
                    table_cell_font(row = 2, column = 1,
                                    face = 'bold', size = 11) %>%
                    tab_add_title(text = 'b', face = "bold", size = 16)

ggarrange(panel_boxplots, panel_comm_table, nrow = 2)

# Export to Cytoscape -----
clr <- flatten_beta_div(clr_dist_matrix)
clr$edge <- apply(clr, 1, 
                  function(x) paste(x[1], x[2], sep = '-'))

ig.df <- as_data_frame(ig, what = c('edges', 'vertices', 'both'))
ig.df$edge <- apply(ig.df, 1, function(x) paste(x[1], x[2], sep = '-'))

clr_2 <- clr %>% filter(edge %in% ig.df$edge)

# write.table(clr_2, 'DATA/cytoscape/aitchison-65.csv',
#               sep = ',', quote = FALSE, row.names = FALSE)

# do the same thing for the 5 communities:
# helper function:
export_cytoscape <- function(comm, graph, dist_matrix, path){
  g <- induced.subgraph(graph, comm)
  df <- as_data_frame(g, what = c('edges', 'vertices', 'both'))
  df$edge <- apply(df, 1, function(x) paste(x[1], x[2], sep = '-'))
  
  tmp <- flatten_beta_div(dist_matrix)
  tmp$edge <- apply(tmp, 1, function(x) paste(x[1], x[2], sep = '-'))
  
  tmp2 <- tmp %>% filter(edge %in% df$edge)
  
  write.table(tmp2, path,
              sep = ',', quote = FALSE, row.names = FALSE)
}

# create the files:
# for (i in 1:length(ig_comms)) {
#   filepath <- paste0('DATA/cytoscape/cluster_', i, '.csv')
#   community <- ig_comms[[i]]
#   export_cytoscape(community, ig, clr_dist_matrix, filepath)
# }

# and the metadata file:
# write.table(metadata.ig, file = 'DATA/cytoscape/metadata.csv',
#               quote = FALSE, sep = ',', row.names = FALSE, col.names = TRUE)

# Network analysis -----
## Degree distribution ----
degs <- degree(ig)

# counts the frequencies of each degree
deg_hist <- as.data.frame(table(degs))

# converts the first column to numbers:
deg_hist$degs <- as.numeric(levels(deg_hist$degs))[deg_hist$degs]

# plot:
panel_degree <- ggplot(deg_hist, aes(x = degs, y = Freq)) +
                        geom_point(size = 2, col = 'grey50') +
                        scale_x_continuous("Degree", trans = 'log10') +
                        scale_y_continuous("Frequency", trans = 'log10') +
                        ggtitle("a") +
                        theme_bw() +
                        theme(axis.text = element_text(size = 12),
                              axis.title = element_text(size = 14),
                              plot.title = element_text(size = 16, 
                                                        face = "bold"))

## Betweenness centrality -----
bet <- betweenness(ig)
bet_hist <- as.data.frame(table(bet))
bet_hist$bet <- as.numeric(levels(bet_hist$bet))[bet_hist$bet]

betweenness_df <- data.frame(cbind(degree = degs, 
                                   betweenness = betweenness(ig)))
betweenness_df <- betweenness_df[order(desc(degs)), ]

panel_betweenness <- ggplot(bet_hist, aes(bet)) +
                      geom_histogram(binwidth = 100) +
                      scale_x_continuous("Betweennes centrality") +
                      scale_y_continuous("Frequency") +
                      ggtitle("b") +
                      theme_bw() +
                      theme(axis.text = element_text(size = 12),
                            axis.title = element_text(size = 14),
                            plot.title = element_text(size = 16, face = "bold"))

## other measures ------
measures_ig <- rbind(gorder(ig),
                      gsize(ig),
                      max(degs),
                      mean(degs),
                      edge_density(ig),
                      length(get_diameter(ig)),
                      mean_distance(ig), # average path length
                      transitivity(ig)) # clustering coefficient

measures_ig <- data.frame(measures_ig)
rownames(measures_ig) <- c('Order', 'Size', 'Max degree', 'Mean degree',
                           'Graph density', 'Diameter', 
                           'Average path length', 'Clustering coefficient')
colnames(measures_ig) <- 'Network'
measures_ig$Network <- round(measures_ig$Network, 2)

## Comparison with a random network -----
# build random graph with the same order and size:
random_g <- erdos.renyi.game(n = gorder(ig), 
                             p.or.m = gsize(ig), 
                             type = "gnm")

# Degree distribution:
degs <- degree(random_g)
deg_hist <- as.data.frame(table(degs))

deg_hist$degs <- as.numeric(levels(deg_hist$degs))[deg_hist$degs]

panel_degree_random <- ggplot(deg_hist, aes(x = degs, y = Freq)) +
  geom_point(size = 2, col = 'grey50') +
  scale_x_continuous("Degree", trans = 'log10') +
  scale_y_continuous("Frequency", trans = 'log10') +
  ggtitle("a") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"))

# Betweenness centrality distribution:
bet <- betweenness(random_g)
bet_hist <- as.data.frame(table(bet))

bet_hist$bet <- as.numeric(levels(bet_hist$bet))[bet_hist$bet]


betweenness_df <- data.frame(cbind(degree = degs, 
                                   betweenness = betweenness(ig)))
betweenness_df <- betweenness_df[order(desc(degs)), ]

panel_betweenness_random <- ggplot(bet_hist, aes(bet)) +
  geom_histogram(binwidth = 10) +
  scale_x_continuous("Betweennes centrality") +
  scale_y_continuous("Frequency") +
  ggtitle("b") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"))

ggpubr::ggarrange(panel_degree_random, panel_betweenness_random)

# other measures:
measures_random <- rbind(gorder(random_g),
                         gsize(random_g),
                         max(degs),
                         mean(degs),
                         edge_density(random_g), 
                         length(get_diameter(random_g)),
                         mean_distance(random_g),
                         transitivity(random_g))
measures_random <- data.frame(measures_random)
rownames(measures_random) <- c('Order', 'Size', 'Max degree', 'Mean degree',
                               'Graph density', 'Diameter', 
                               'Average path length', 'Clustering coefficient')

colnames(measures_random) <- 'Whole'
# add to the table:
measures_ig$Random <- round(measures_random$Whole, 2)
measures_ig

# Assortativity analysis ----
ast_clep <- assortativity(ig, V(ig)$Clostridium_leptum)
ast_gpam <- assortativity(ig, V(ig)$Gordonibacter_pamelaeae)
ast_elen <- assortativity(ig, V(ig)$Eggerthella_lenta)
ast_cint <- assortativity(ig, V(ig)$Collinsella_intestinalis)
ast_chao1 <- assortativity(ig, V(ig)$Chao1)
ast_shannon <- assortativity(ig, V(ig)$Shannon)
ast_simpson <- assortativity(ig, V(ig)$Simpson)
ast_study <- assortativity_nominal(ig, V(ig)$study_name)
ast_group <- assortativity_nominal(ig, V(ig)$group)
ast_metab <- assortativity_nominal(ig, V(ig)$metab)
ast_obese <- assortativity_nominal(ig, V(ig)$obesity_status)

ast.df <- data.frame(rbind(ast_clep, ast_gpam, ast_elen, ast_cint,
                           ast_chao1, ast_shannon, ast_simpson, 
                           ast_study, ast_group, ast_metab, ast_obese))
colnames(ast.df) <- 'Assortativity'
rownames(ast.df) <- c('Clostridium leptum', 'Gordonibacter pamelaeae',
                      'Collinsella intestinalis', 'Eggerthella lenta',
                      'Chao1', 'Shannon', 'Simpson',
                      'Study name', 'Group', 
                      'Metabolic category', 'Obesity status')

ast.df$Assortativity <- round(ast.df$Assortativity, 2)

## Network randomizations ------
randomizations <-  data.frame(matrix(NA, nrow = nrow(ast.df), ncol = 100),
                              row.names = c('Clostridium leptum', 
                                            'Gordonibacter pamelaeae',
                                            'Collinsella intestinalis', 
                                            'Eggerthella lenta',
                                            'Chao1', 'Shannon', 'Simpson',
                                            'Study name', 'Group', 
                                            'Metabolic category', 
                                            'Obesity status'))

# generate 100 random networks and calculate assortativities:
for (i in seq(1:100)) {
  randomized <- rewire(ig, each_edge(prob = 0.01))
  # calculate assortativities:
  ast_clep <- assortativity(randomized, V(randomized)$Clostridium_leptum)
  ast_gpam <- assortativity(randomized, V(randomized)$Gordonibacter_pamelaeae)
  ast_elen <- assortativity(randomized, V(randomized)$Eggerthella_lenta)
  ast_cint <- assortativity(randomized, V(randomized)$Collinsella_intestinalis)
  ast_chao1 <- assortativity(randomized, V(randomized)$Chao1)
  ast_shannon <- assortativity(randomized, V(randomized)$Shannon)
  ast_simpson <- assortativity(randomized, V(randomized)$Simpson)
  ast_study <- assortativity_nominal(randomized, V(randomized)$study_name)
  ast_group <- assortativity_nominal(randomized, V(randomized)$group)
  ast_metab <- assortativity_nominal(randomized, V(randomized)$metab)
  ast_obese <- assortativity_nominal(randomized, V(randomized)$obesity_status)
  col <- i
  
  randomizations[ , i] <- c(ast_clep, ast_gpam, ast_elen, ast_cint,
                            ast_chao1, ast_shannon, ast_simpson, 
                            ast_study, ast_group, ast_metab, ast_obese)
}
randomizations

# obtain p-values from the randomizations: 
pvals <- vector()
for (i in seq(1:nrow(ast.df))) {
  pvals <- c(pvals, sum(abs(randomizations[i, ]) >= abs(ast.df[i, 1]))/100)
}
names(pvals) <- rownames(ast.df)
pvals
# adjust via Benjamini-Hochberg:
adjs <- as.data.frame(pvals) %>% 
          adjust_pvalue(p.col = 'pvals', method = 'BH')

# we could get only the features with significant assortativities:
ast.keep <- data.frame(ast.df[adjs$pvals.adj < 0.05, ])
rownames(ast.keep) <- rownames(ast.df)[adjs$pvals.adj < 0.05]
colnames(ast.keep) <- 'Assortativity'

# get panel for the figure:
panel_assortativity <- ggtexttable(ast.df %>% rownames_to_column('Feature'), 
                            rows = NULL) %>% 
                        table_cell_font(row = 2:5, column = 1, 
                                        face = 'italic', size = 11) %>% 
                        table_cell_font(row = c(6, 9, 10, 11), column = 1, 
                                        face = 'bold', size = 11) %>%
                        table_cell_font(row = 4, column = 1, 
                                        face = 'bold.italic', size = 11) %>%
                        tab_add_title(text = 'c', face = "bold", size = 16)

# build the figure:
panel_deg_bet <- ggarrange(panel_degree, panel_betweenness, nrow = 2)
ggarrange(panel_deg_bet, NULL, panel_assortativity, 
          widths = c(1, 0.05, 1), nrow = 1)

# Cluster analysis ----
# We still need to get our measures for the 5 clusters:
cluster_analyser <- function(comm, graph){
  g <- induced.subgraph(graph, comm)
  deg_tmp <- degree(g)
  measures_g <- c(gorder(g),
                  gsize(g),
                  max(deg_tmp),
                  mean(deg_tmp), 
                  edge_density(g), 
                  length(get_diameter(g)),
                  mean_distance(g), 
                  transitivity(g))
  return(round(measures_g, 2))
  }

# do it for all clusters and add it to the table:
for (i in 1:length(ig_comms)) {
  community <- ig_comms[[i]]
  m <- cluster_analyser(community, ig)
  measures_ig <- cbind(measures_ig, m)
  }

colnames(measures_ig) <- c('Network', 'Random',
                           seq(1:length(ig_comms)))

## nice display as a gt table:
measures_ig <- measures_ig %>% rownames_to_column()
measures_gt <- gt(measures_ig) %>% 
                tab_spanner(label = 'Communities',
                            columns = 4:8)  %>%
                tab_style(style = cell_text(weight = 'bold'), 
                          locations = list(cells_column_labels(),
                                           cells_stub(),
                                           cells_column_spanners())) %>% 
                fmt_number(columns = everything(),
                           rows = c(1:3, 6), decimals = 0)

# gtsave(measures_gt,
#        "plots/exploratory-network-65.png",
#        zoom = 1)