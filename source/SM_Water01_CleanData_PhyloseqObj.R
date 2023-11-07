## Packages
library(here)
source(here("source", "libraries.R"))
library(phyloseq)
library(reshape2)
##### Load data #### 

## ....Abundance table ####
ab_table_raw <- read.csv(file = here("Data/Water","Mesocosm_Water_18S_oct2023.ASV_table_norarefaction_dnNA.tsv"), dec = ".", sep = "\t", header = T, row.names = 1, comment.char = "") #20078 objs in 299 vars

ab_table <- ab_table_raw[,!names(ab_table_raw) %in% c("taxonomy")]# Remove sequence and taxonomy columns to keep only ab. data 
names(ab_table)

## ....Taxonomy ####
# Put the taxonomy information in another table 
# (drop= F trick to keep the row.names) 

taxo <- ab_table_raw[,"taxonomy", drop=F] 

# Split into different columns 
taxo_split <- with(taxo, cbind(row.names(taxo), colsplit(taxo$taxonomy, pattern = "\\;", names = c('Domain', 'Phylum','Class','Order','Family','Genus','Species'))))

# Change row.names to ASV sequence and remove ASV sequence from dataframe
row.names(taxo_split) <- taxo_split$`row.names(taxo)`
taxo_split_clean <- taxo_split[,-1]

# Remove all the p__, c__ etc
taxo_split_clean <- mutate(taxo_split_clean, across(where(is.character), str_remove, ".*__"))

# Empty string to NA
taxo_split_clean[taxo_split_clean == ''] <- NA


## Fill NAs with last available taxon rank
tax_clean_filled <- taxo_split_clean |> 
    t() |> 
    as.data.frame() |> 
    fill(everything(), .direction = "down") |> 
    t() |> 
    as.data.frame()


###.... Metadata ####
meta.raw <- read.csv(file = here("Data/Water","metadata.tsv"), dec = ".", header = T, row.names = 1, sep = "\t", comment.char = "") #load
str(meta.raw)
meta.raw[sapply(meta.raw, is.character)] <- lapply(meta.raw[sapply(meta.raw, is.character)], as.factor) # did it work? Check with str(meta)

#replace _ by . to match sample names from abundance dataset
rownames(meta.raw) <- stringr::str_replace_all(rownames(meta.raw), '[-]', '.')

meta <- meta.raw


# Reorder Time in chronological order
levels(meta$Time)
meta$Time<- factor(meta$Time, levels=c('D-0', 'D8', 'D15', 'D34', 'D61',"D77"))
# Renq;ing time to have coherent dates with other sampling (sed,roots, rhizo)
levels(meta$Time) <- list(D0="D-0", D10="D8", D16 = "D15", D35= "D34", D62='D61',D77="D77")

#### Phyloseq object####

ASV <- otu_table(ab_table, taxa_are_rows = TRUE)
SDAT <- sample_data(meta)
TAX <- tax_table(as.matrix(taxo_split_clean))

ps_intermediate<- merge_phyloseq(ASV,SDAT,TAX)
ps_intermediate # 850 taxa in 141 samples


ps = prune_samples(sample_sums(ps_intermediate) > 500, ps_intermediate)
setdiff(sample_names(ps_intermediate),sample_names(ps))
# Removes one sample (M7G.Jan.4.2022) 
ab_table$M7G.Jan.4.2022 #Sample has only 5 reads

ps <- prune_taxa(taxa_sums(ps)>0, ps)
ps
getwd()
# Save object 
save(ps, file = here("RData","ps_water_obj.RData"))
