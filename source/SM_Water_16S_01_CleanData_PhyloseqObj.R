## Packages
library(here)
source(here("source", "libraries.R"))
library(phyloseq)
library(reshape2)
##### Load data #### 

## ....Abundance table ####
ab_table_raw <- read.csv(file = here("Data/Water/16S/","Mesocosm_Water_16S_oct2023.ASV_table_norarefaction_dnNA.tsv"), dec = ".", sep = "\t", header = T, row.names = 1, comment.char = "")
# 2462 obs in 145 vars

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
meta.raw <- read.csv(file = here("Data/Water/18S","metadata.tsv"), dec = ".", header = T, row.names = 1, sep = "\t", comment.char = "") #load
str(meta.raw)
meta.raw[sapply(meta.raw, is.character)] <- lapply(meta.raw[sapply(meta.raw, is.character)], as.factor) # did it work? Check with str(meta)

# Assign alias as rownames because some sample have .tar in there names
rownames(meta.raw) <- meta.raw$alias
#replace _ by . to match sample names from abundance dataset
rownames(meta.raw) <- stringr::str_replace_all(rownames(meta.raw), '[_]', '.')
names(ab_table)

(rownames(meta.raw)%in%names(ab_table))

meta <- meta.raw

# Reorder Time in chronological order
levels(meta$Time)
meta$Time<- factor(meta$Time, levels=c('D-0', 'D8', 'D15', 'D34', 'D61',"D77"))
# Renq;ing time to have coherent dates with other sampling (sed,roots, rhizo)
levels(meta$Time) <- list(D0="D-0", D10="D8", D16 = "D15", D35= "D34", D62='D61',D77="D77")

#


#### Phyloseq object####

ASV <- otu_table(ab_table, taxa_are_rows = TRUE)
SDAT <- sample_data(meta)
TAX <- tax_table(as.matrix(tax_clean_filled))



ps_intermediate<- merge_phyloseq(ASV,SDAT,TAX)
ps_intermediate # 2462 taxa in 144 samples

plot_bar(ps_intermediate)
# 5 samples have a low read abundance 

sort(sample_sums(ps_intermediate)) # see the number of reads per sample
ps = prune_samples(sample_sums(ps_intermediate) > 500, ps_intermediate)
setdiff(sample_names(ps_intermediate),sample_names(ps))
# Removes the 5 samples
# M2G.30.11.21 ; M10E.Jan.31.2022  ; M12E.10.12.21 ; M9E.434.Feb.16.2022 ; M7E.432.Feb.16.2022


# Save object 
save(ps, file = here("RData","ps_16S_water_obj.RData"))
