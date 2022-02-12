library(tidyverse)
library(phyloseq)

# load taxonomy table and patient metadata:
# NOTE: unpublished data - not available 
taxa <- read_tsv("DATA/taxonomy.tsv")
metadata <- read_csv("DATA/data_to_lefse.csv")

# Convert classes to numbers
# we could use the classes instead but we want the colors
# to be assigned in a specific way: red for CRC, green 
# for controls, blue for adenoma
metadata$Study_Condition <- gsub("CRC", "1", metadata$Study_Condition)
metadata$Study_Condition <- gsub("Control", "2",metadata$Study_Condition)
metadata$Study_Condition <- gsub("adenoma", "3",metadata$Study_Condition)
# convert to factor:
metadata$Study_Condition <- factor(metadata$Study_Condition, 
                                   levels=c("1","2","3"))

# Load relative abundances table:
abundance <- read_tsv("DATA/vect_CRC.txt")
# Filter Homo sapiens MSP:
taxa <- taxa %>% filter(species != 'Homo sapiens')

# A few fixes to format the table as LEfSe requires:
taxa <- taxa %>% unite(kingdom, phylum, class, order, 
                       family, genus, species, 
                       sep="|")
colnames(taxa) <- c("MSP","OTU")
table <- taxa %>% left_join(abundance) %>% select(-MSP) %>%
         gather(-OTU, key="Sample_ID",value="Abundance")
table <- table %>% right_join(metadata)
table$OTU <- gsub(" / ", "_",table$OTU)
table$OTU <- gsub(" & ", "_",table$OTU)
table$OTU <- gsub(" ", "_",table$OTU)

lefse <- table %>% 
          pivot_wider(names_from = OTU, values_from = Abundance, 
                      values_fn = mean)

lefse <- as.data.frame(t(lefse))
lefse[is.na(lefse)] <- 0 # Converts NAs to 0s
lefse <- rownames_to_column(lefse,var = "var")

# Write output table:
write_tsv(lefse,"lefse.txt", col_names = F)
