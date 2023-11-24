## Packages
library(here)
source(here("source", "libraries.R"))
library(phyloseq)
library(reshape2)
##### Load data #### 

## ....Abundance table ####

#Load abundance table 
ab_table_raw <- read.csv(file = here("data","MesocosmExp2_AbTable_raw.tsv"), 
                         dec = ".", sep = "\t", header = T, row.names = 1, comment.char = "")
#20078 objs in 299 vars

# Remove X from sample names 
names(ab_table_raw) <- sub("X", "", names(ab_table_raw))
rownames(ab_table_raw) <- ab_table_raw$Sequence # replace row.names code by Sequence

# Remove sequence and taxonomy columns to keep only ab. data 
ab_table <- ab_table_raw[,!names(ab_table_raw) %in% c("Sequence", "taxonomy")]

# The negative controls reads found in common with the rest of the dataset 
# represented a low amount of reads. We decided that contamination wasn't an issue.
# Further analysis were done by removing the negative controls
No_ctrl <- !grepl("Ctrl", names(ab_table))
ab_table <- ab_table[,No_ctrl]
names(ab_table)

## ....Taxonomy ####
# Put the taxonomy information in another table 
# (drop= F trick to keep the row.names) 

taxo <- ab_table_raw[,"taxonomy", drop=F] 

# Split each taxonomy level into different columns 
taxo_split <- with(taxo, cbind(row.names(taxo),
                               colsplit(taxo$taxonomy,
                                        pattern = "\\;", names = c('Domain', 'Phylum','Class','Order','Family','Genus','Species'))))

# Change row.names to ASV sequence and remove ASV sequence from dataframe
row.names(taxo_split) <- taxo_split$`row.names(taxo)`
taxo_split_clean <- taxo_split[,-1]

# Remove all the p__, c__ etc
taxo_split_clean <- mutate(taxo_split_clean, across(where(is.character), str_remove, ".*__"))

# Empty string to NA
taxo_split_clean[taxo_split_clean == ''] <- NA

# Check for Chloroplast, Mitochondria and Non-bacterial sequences
# None in this dataset

###.... Metadata ####

#Load metadeta table 
meta.raw <- read.csv(file = here("data","metadata.tsv"), dec = ".", header = T, row.names = 1, sep = "\t", comment.char = "") #load
str(meta.raw) #Check structure

# Change characters to factors 
meta.raw[sapply(meta.raw, is.character)] <- lapply(meta.raw[sapply(meta.raw, is.character)], as.factor) 
# Did it work? Check with str(meta)

# Change in rownames
rownames(meta.raw) <- stringr::str_replace_all(rownames(meta.raw), '[-]', '.')

# Remove the ctrl samples
meta <- subset(meta.raw,!grepl("Ctrl",rownames(meta.raw))) 
rownames(meta)

# Drop unused factor levels while keeping the rownames 
meta[] <- as.data.frame(lapply(meta, function(x) if(is.factor(x)) droplevels(x) else x)) 
levels(meta$Plant.type)# Sanity check

# Reorder Time in chronological order
meta$Time<- factor(meta$Time, levels=c('D0', 'D10', 'D16', 'D35', 'D62'))
levels(meta$Time)

# Reorder Sample Type
meta$Sample.type<- factor(meta$Sample.type, levels=c('Roots', 'Rhizosphere', 'Sediment'))
levels(meta$Sample.type)

# Homogeneisation needed for T°
levels(meta$Temperature)
meta$Temperature <- sub("/","",meta$Temperature)
meta$Temperature <- as.factor(meta$Temperature)
levels(meta$Temperature)

#### Phyloseq object####

ASV <- otu_table(ab_table, taxa_are_rows = TRUE)
SDAT <- sample_data(meta)
TAX <- tax_table(as.matrix(taxo_split_clean))

ps_intermediate<- merge_phyloseq(ASV,SDAT,TAX)
ps_intermediate # 20078 taxa in 283 samples


ps = prune_samples(sample_sums(ps_intermediate) > 500, ps_intermediate)
setdiff(sample_names(ps_intermediate),sample_names(ps))
# Removes one rhizosphere sample with less than 500 reads

ps <- prune_taxa(taxa_sums(ps)>0, ps)
ps
# Removes 15 ASVs with 0 reads

# Save object 
save(ps, file = here("RData","ps_obj.RData"))
