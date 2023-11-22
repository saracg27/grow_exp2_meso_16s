### Packages#####
library(here)
source(here("source", "libraries.R"))
library(phyloseq)
library(viridis)
library(metacoder)
library(dplyr)
library(ggpp)

#### Data ####

load(here("Rdata","ps_obj.RData"))

ps_sed <-  subset_samples(ps, Sample.type =="Sediment")
ps_sed <- prune_taxa(taxa_sums(ps_sed)>0, ps_sed)
ps_sed
# 120 samples 18851 taxa

## ______ #####
### DIFFERENTIAL ABUNDANCE ####
###____####
###Metacoder ####

####Metacoder Object####
# Taxonomic agglomeration
ps_glom_order <- tax_glom(ps_sed,taxrank="Order")
ps_glom_class <- tax_glom(ps_sed,taxrank="Class")

# Choose which agglomeration object you want ot use
ps_obj <- ps_glom_class

# Relative abundance transformation
ps.relab <- transform_sample_counts(ps_obj,function(x) x/sum(x))

# Transform the phyloseq object in a metacoder object
meta_obj <- parse_phyloseq(ps.relab) 

####Heat tree #####

### Sample Type ###
meta_obj$data$tax_abund <- calc_taxon_abund(meta_obj, "otu_table")


meta_obj$data$sample_data$Time <- factor(as.character(meta_obj$data$sample_data$Time), 
                                         levels = c("D0","D10","D16","D35","D62"),
                                         ordered = TRUE)


meta_obj$data$diff_table <- compare_groups(meta_obj, data = "tax_abund",
                                           cols = meta_obj$data$sample_data$sample_id,
                                           groups = meta_obj$data$sample_data$Time)# What columns of sample data to use
# Compares the taxonomic difference based on ASV abundance

# Adjust the p_value with FDR correction
meta_obj$data$diff_table$wilcox_FDR_p_value <- p.adjust(meta_obj$data$diff_table$wilcox_p_value,
                                                        method = "fdr")

# See the range of FDR p_values
range(meta_obj$data$diff_table$wilcox_FDR_p_value, finite = TRUE) 

# Set the log2MedianRatio to 0 if the corrected p values are bigger
# than 0.05. It's a trick for plotting purposes as you don't want to 
# highlight unsignificant results. 
meta_obj$data$diff_table$log2_median_ratio[meta_obj$data$diff_table$wilcox_FDR_p_value > 0.05] <- 0

# Only label significant differences
per_taxon_fold_changes <- obs(meta_obj, data = 'diff_table', value = 'log2_median_ratio')
per_taxon_max_change <- unlist(lapply(per_taxon_fold_changes, function(tax_changes) max(abs(tax_changes))))
meta_obj_simp <- filter_taxa(meta_obj, per_taxon_max_change !=0, supertaxa = TRUE, reassign_obs = c(diff_table = FALSE))


# Plot the tree
heat_tree_matrix (meta_obj_simp,
                  data="diff_table",
                  node_size = n_obs, # n_obs is a function that calculates, in this case, the number of ASV per taxon
                  node_label = taxon_names,
                  node_color = log2_median_ratio, # A column from `obj$data$diff_table`
                  node_color_range = diverging_palette(), # The built-in palette for diverging data                 
                  node_color_interval = c(-10, 10), # The range of `log2_median_ratio` to display
                  edge_color_interval = c(-10,10), # The range of `log2_median_ratio` to display
                  node_size_axis_label = "Number of orders",
                  node_color_axis_label = "Log2 ratio median proportions",
                  repel_labels = TRUE,
                  key_size = 0.8,
                  label_small_trees = F,
                  layout = "davidson-harel", # The primary layout algorithm
                  initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                  output_file = here("Results","Figures","Sed_Time_DiffTree.pdf")) # Saves the plot as a pdf file

### Sanity check of metacoder results ###
Diff_table <- as.data.frame(meta_obj$data$diff_table)
tax <- meta_obj$data$tax_data
tax_order <- tax[,c("taxon_id","Order")] 
# Merge taxonomy table and diff. abundance results
Diff_table_tax <- merge(tax_order,Diff_table,'taxon_id')

#Subset D0 samples
ps_D0 <- subset_samples(ps_sed,Time=="D0")
ps_D0 <- prune_taxa(taxa_sums(ps_D0)>0, ps_D0)
# 3735 Taxa in 24 samples 

#Subset D62 samples
ps_D62 <- subset_samples(ps_sed,Time=="D62")
ps_D62 <- prune_taxa(taxa_sums(ps_D62)>0, ps_D62)
# 13322 Taxa in 24 samples


### According to the tree,
# Burkholderiales should be more abundant in D0 than D62

ps_D0_burk <- subset_taxa(ps_D0,Order=="Burkholderiales")
ps_D0_burk <- prune_taxa(taxa_sums(ps_D0_burk)>0, ps_D0_burk)
# 352 Burk in 24 samples in D0
Ab_burk_D0 <- sum(rowSums(otu_table(ps_D0_burk))) # 214 050 reads Burk in D0

ps_D62_burk <- subset_taxa(ps_D62,Order=="Burkholderiales")
ps_D62_burk <- prune_taxa(taxa_sums(ps_D62_burk)>0, ps_D62_burk)
# 589 Burk in 24 samples in D62
Ab_burk_D62 <- sum(rowSums(otu_table(ps_D62_burk))) # 78 883 reads Burk in D62

Diff_Burk <- subset(Diff_table_tax, Order=="Burkholderiales" & wilcox_FDR_p_value<0.05)

#### Alicyclobacillales should be only be found in D66
ps_D0_alic <- subset_taxa(ps_D0,Order=="Alicyclobacillales")
ps_D0_alic <- prune_taxa(taxa_sums(ps_D0_alic)>0, ps_D0_alic)
# 3 alic in 24 samples in D0

Ab_alic_D0 <- sum(rowSums(otu_table(ps_D0_alic))) # 140 reads Alic in D0

check_alic <- as.data.frame(otu_table(ps_D0_alic))
check_alic_sum <- as.data.frame(colSums(check_alic))
median(check_alic_sum$`colSums(check_alic)`)


ps_D62_alic <- subset_taxa(ps_D62,Order=="Alicyclobacillales")
ps_D62_alic <- prune_taxa(taxa_sums(ps_D62_alic)>0, ps_D62_alic)
# 18 alic in 24 samples in D62
Ab_alic_D62 <- sum(rowSums(otu_table(ps_D62_alic))) # 1728 reads Alic in D62

check_alic62 <- as.data.frame(otu_table(ps_D62_alic))
check_alic_sum62 <- as.data.frame(colSums(check_alic62))
median(check_alic_sum62$`colSums(check_alic62)`)

Diff_alic <- subset(Diff_table_tax, Order=="Alicyclobacillales" & wilcox_FDR_p_value<0.05)
### So -inf / inf means that the median abundance of one of the groups is 0
#### Doesn't mean that there are no seq of this taxa in this group



### Test with low read sample removal
ps_sed_sub <- prune_samples(sample_sums(ps_sed)>10000,ps_sed)
ps_sed_sub <- prune_taxa(taxa_sums(ps_sed_sub)>0, ps_sed_sub)
# Loss of 5 D0 samples (8,15,27,30,31) and 13 ASVs only present in those samples


ps_glom_ord <- tax_glom(ps_sed_sub,taxrank="Order")

ps.relab <- transform_sample_counts(ps_glom_ord,function(x) x/sum(x))

levels(sample_data(ps.relab)$Time)
# First convert abundance to relative abundance

meta_obj <- parse_phyloseq(ps.relab) # transforms the phyloseq object in a metacoder object


### Sample Type ###
meta_obj$data$tax_abund <- calc_taxon_abund(meta_obj, "otu_table")


meta_obj$data$sample_data$Time <- factor(as.character(meta_obj$data$sample_data$Time), 
                                         levels = c("D0","D10","D16","D35","D62"),
                                         ordered = TRUE)


meta_obj$data$diff_table <- compare_groups(meta_obj, data = "tax_abund",
                                           cols = meta_obj$data$sample_data$sample_id,
                                           groups = meta_obj$data$sample_data$Time)# What columns of sample data to use
# Compares the taxonomic difference based on ASV abundance


meta_obj$data$diff_table$wilcox_FDR_p_value <- p.adjust(meta_obj$data$diff_table$wilcox_p_value,
                                                        method = "fdr")
range(meta_obj$data$diff_table$wilcox_FDR_p_value, finite = TRUE) 

meta_obj$data$diff_table$log2_median_ratio[meta_obj$data$diff_table$wilcox_FDR_p_value > 0.05] <- 0
# Sets the log2MedianRatio to 0 if the corrected p values are bigger than 0.05

heat_tree_matrix (meta_obj,
                  data="diff_table",
                  node_size = n_obs, # n_obs is a function that calculates, in this case, the number of ASV per taxon
                  node_label = taxon_names,
                  node_color = log2_median_ratio, # A column from `obj$data$diff_table`
                  node_color_range = diverging_palette(), # The built-in palette for diverging data                 
                  node_color_interval = c(-10, 10), # The range of `log2_median_ratio` to display
                  edge_color_interval = c(-10,10), # The range of `log2_median_ratio` to display
                  node_size_axis_label = "Number of orders",
                  node_color_axis_label = "Log2 ratio median proportions",
                  repel_labels = TRUE,
                  label_small_trees = T,
                  layout = "davidson-harel", # The primary layout algorithm
                  initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                  output_file = here("Results","Figures","Sed_Time_DiffTree_no5.pdf")) # Saves the plot as a pdf file


# Some changes exemple D62-D0 comparison on Proteobacteria 
# Overall remains quite similar. 




#### Potential NA degraders #####
NA_list <- c("Pseudomonas",
             "Rhodoferax",
             "Rhodococcus",
             "Acidovorax",
             "Alcaligenes",
             "Acinetobacter",
             "Kurthia",
             "Methonabrevibacter",
             "Methanolinea",
             "Methanoregula")

## In sediments
ps_NA_list <- subset_taxa(ps_sed,Genus%in%NA_list)

ps.relab <- transform_sample_counts(ps_NA_list,function(x) x/sum(x))
# First convert abundance to relative abundance

meta_obj <- parse_phyloseq(ps.relab) # transforms the phyloseq object in a metacoder object

####__Diff heat tree ##

### Time ###
meta_obj$data$tax_abund <- calc_taxon_abund(meta_obj, "otu_table")


meta_obj$data$sample_data$Time <- factor(as.character(meta_obj$data$sample_data$Time), 
                                         levels = c("D0","D10","D16","D35","D62"),
                                         ordered = TRUE)


meta_obj$data$diff_table <- compare_groups(meta_obj, data = "tax_abund",
                                           cols = meta_obj$data$sample_data$sample_id,
                                           groups = meta_obj$data$sample_data$Time)# What columns of sample data to use
# Compares the taxonomic difference based on ASV abundance

meta_obj$data$diff_table$wilcox_FDR_p_value <- p.adjust(meta_obj$data$diff_table$wilcox_p_value,
                                                        method = "fdr")
range(meta_obj$data$diff_table$wilcox_FDR_p_value, finite = TRUE) 

test <- as.data.frame(meta_obj$data$diff_table)

meta_obj$data$diff_table$log2_median_ratio[meta_obj$data$diff_table$wilcox_FDR_p_value > 0.05] <- 0
# Sets the log2MedianRatio to 0 if the corrected p values are bigger than 0.05

heat_tree_matrix (meta_obj,
                  data="diff_table",
                  node_size = n_obs, # n_obs is a function that calculates, in this case, the number of ASV per taxon
                  node_label = taxon_names,
                  node_color = log2_median_ratio, # A column from `obj$data$diff_table`
                  node_color_range = diverging_palette(), # The built-in palette for diverging data                  node_color_trans = "linear", # The default is scaled by circle area
                  node_color_interval = c(-10, 10), # The range of `log2_median_ratio` to display
                  edge_color_interval = c(-10,10), # The range of `log2_median_ratio` to display
                  node_size_axis_label = "Number of orders",
                  node_color_axis_label = "Log2 ratio median proportions",
                  repel_labels = TRUE,
                  layout = "davidson-harel", # The primary layout algorithm
                  initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                  output_file = here("Results","Figures","NA_taxaSed_Time_DiffTree.pdf")) # Saves the plot as a pdf file



#### Sample Type Comparison ####

ps.relab <- transform_sample_counts(ps_NA_list,function(x) x/sum(x))
# First convert abundance to relative abundance

meta_obj <- parse_phyloseq(ps.relab) # transforms the phyloseq object in a metacoder object
#### Diff heat tree ##

meta_obj$data$tax_abund <- calc_taxon_abund(meta_obj, "otu_table")


meta_obj$data$diff_table <- compare_groups(meta_obj, data = "tax_abund",
                                           cols = meta_obj$data$sample_data$sample_id,
                                           groups = meta_obj$data$sample_data$Sample.type)# What columns of sample data to use
# Compares the taxonomic difference based on ASV abundance

meta_obj$data$diff_table$wilcox_FDR_p_value <- p.adjust(meta_obj$data$diff_table$wilcox_p_value,
                                                        method = "fdr")
range(meta_obj$data$diff_table$wilcox_FDR_p_value, finite = TRUE) 

test <- as.data.frame(meta_obj$data$diff_table)

meta_obj$data$diff_table$log2_median_ratio[meta_obj$data$diff_table$wilcox_FDR_p_value > 0.05] <- 0
# Sets the log2MedianRatio to 0 if the corrected p values are bigger than 0.05

heat_tree_matrix (meta_obj,
                  data="diff_table",
                  node_size = n_obs, # n_obs is a function that calculates, in this case, the number of ASV per taxon
                  node_label = taxon_names,
                  node_color = log2_median_ratio, # A column from `obj$data$diff_table`
                  node_color_range = diverging_palette(), # The built-in palette for diverging data                  node_color_trans = "linear", # The default is scaled by circle area
                  node_color_interval = c(-10, 10), # The range of `log2_median_ratio` to display
                  edge_color_interval = c(-10,10), # The range of `log2_median_ratio` to display
                  node_size_axis_label = "Number of orders",
                  node_color_axis_label = "Log2 ratio median proportions",
                  repel_labels = TRUE,
                  layout = "davidson-harel", # The primary layout algorithm
                  initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                  output_file = here("Results","Figures","NA_taxa_Compartments_DiffTree.pdf")) # Saves the plot as a pdf file


####_______#####
###Differential ab. on D0-D62 #####
##### Deseq 2 #####
library(DESeq2)

# Subset d0 and d62 samples
ps_sed_D0D62 <- subset_samples(ps_sed,Time=="D0"|Time=="D62")
ps_sed_D0D62 <- prune_taxa(taxa_sums(ps_sed_D0D62)>0, ps_sed_D0D62)

# Phyloseq to deseq object conversion
diagdds = phyloseq_to_deseq2(ps_sed_D0D62, ~ Time)

# Estimate size factor for normalization
diagdds = estimateSizeFactors(diagdds, type = "poscounts")

# Differential abundance analysis with Wald test and local fit
diagdds = DESeq(diagdds,
                test = "Wald",
                fitType="local")
plotDispEsts(diagdds)

# Results of the test, FDR accepted = 5% 
res = results(diagdds, alpha=0.05) 
plotMA(res)

# Ordered and removes the ASVs which have padj = NAs 
res2 = res[order(res$padj, na.last=NA), ] 
plotMA(res2)

# Select significant ASVs
sigtab = res2[(res2$padj < 0.05), ] 
sigtab <- as.data.frame(sigtab)
sigtab <- cbind(rownames(sigtab),sigtab)
colnames(sigtab)[1] <- "OTU"

#Transform phyloseq object in long format 
ps_long <- psmelt(ps_sed_D0D62)
# Remove unrelevant columns
ps_long <- select(ps_long,-c("CFL_ID_Sediments",
                             "Sample.type",
                             "Experiment",
                             "Greenhouse",
                             "Sampling.date",
                             "Mesocosm..",
                             "Water.type",
                             "Plant.type" ,
                             "Greenhouse" ,     
                             "Temperature",      
                             "Material"))
ps_long$OTU <- as.factor(ps_long$OTU)

# Abundance of every ASVs in D0 and D62
ab_d0_D62<- ps_long %>% 
  group_by(OTU, Time) %>% 
  summarise(Sum_Abundance = sum(Abundance))%>%
  pivot_wider(names_from = Time, values_from = Sum_Abundance, values_fill = 0)

# Retrieve taxonomy
tax <- as.data.frame(tax_table(ps_sed_D0D62))
tax$OTU <- rownames(tax)

# Merge taxonomy and object
ab_tax_d0d62 <- merge(ab_d0_D62,tax,by="OTU")

# Merge previous object and significant diff. ab table 
sig_data <- merge(sigtab,ab_tax_d0d62,by="OTU")

# Change charcyers to factors
sig_data[sapply(sig_data, is.character)] <- lapply(sig_data[sapply(sig_data, is.character)], as.factor) # did it work? Check with str(meta)

# Subset D0 and D62
sig_data_plot <- subset(sig_data,sig_data$D0!=0 & sig_data$D62!=0 )

# Plot ASVs regrouped at order levels with  
# #significant differential abundance
# Colors is set for abundance in D62
# size of point is set for abundance in D0
ggplot(sig_data_plot,aes(y=Order, x=log2FoldChange,size=D0,color=D62))+
  geom_jitter(height=0.1)+
  theme_bw()+
  geom_vline(xintercept = 0, color ="red")+
  scale_color_viridis()+
  scale_size_continuous(breaks=c(1,1000,10000,40000))+
  guides(size = guide_legend(reverse=TRUE))


### sanity check 
# Check why Delftia is so abundant in D0
# Delftia <- subset(sig_data_plot$OTU, sig_data_plot$Genus=="Delftia")
# Delftia <- as.character(Delftia)
# 
# # Subset that summed abundance table keeping only Delftia ASV based on row.names (sequence)
# check <- subset(ab_d0_D62,ab_d0_D62$OTU==Delftia)
# # 40742 in D0 ; 516 in D62


#### Metacoder ####
ps_sed_D0D62 <- subset_samples(ps_sed,Time=="D0"|Time=="D62")
ps_sed_D0D62 <- prune_taxa(taxa_sums(ps_sed_D0D62)>0, ps_sed_D0D62)

ps.relab <- transform_sample_counts(ps_sed_D0D62,function(x) x/sum(x))

meta_obj <- parse_phyloseq(ps.relab) # transforms the phyloseq object in a metacoder object

#### Diff heat tree 

### Time ###
meta_obj$data$tax_abund <- calc_taxon_abund(meta_obj, "otu_table")

meta_obj$data$diff_table <- compare_groups(meta_obj, data = "tax_abund",
                                           cols = meta_obj$data$sample_data$sample_id,
                                           groups = meta_obj$data$sample_data$Time)# What columns of sample data to use
# Compares the taxonomic difference based on ASV abundance


meta_obj$data$diff_table$wilcox_FDR_p_value <- p.adjust(meta_obj$data$diff_table$wilcox_p_value,
                                                        method = "fdr")
range(meta_obj$data$diff_table$wilcox_FDR_p_value, finite = TRUE) 

meta_obj$data$diff_table$log2_median_ratio[meta_obj$data$diff_table$wilcox_FDR_p_value > 0.05] <- 0
# Sets the log2MedianRatio to 0 if the corrected p values are bigger than 0.05

res <- meta_obj$data$diff_table
res <- res[res$log2_median_ratio!=0,]

res_tax <- merge(tax,res,by="taxon_id",all=T)
res_tax_D0 <- res_tax[res_tax$log2_median_ratio>0,] # 90
res_tax_D62 <- res_tax[res_tax$log2_median_ratio<0,] # 452

heat_tree(meta_obj,
          node_size = n_obs, # n_obs is a function that calculates, in this case, the number of ASV per taxon
          node_label = taxon_names,
          node_color = log2_median_ratio, # A column from `obj$data$diff_table`
          node_color_range = diverging_palette(), # The built-in palette for diverging data                 
          node_color_interval = c(-10, 10), # The range of `log2_median_ratio` to display
          edge_color_interval = c(-10,10), # The range of `log2_median_ratio` to display
          node_size_axis_label = "Number of ASVs",
          node_color_axis_label = "Log2 ratio median proportions",
          repel_labels = TRUE,
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
          output_file = here("Results","Figures","Sed_D0_D62_DiffTree.pdf")) # Saves the plot as a pdf file



### Unlabel taxons with  log2_median_ratio = 0 
per_taxon_fold_changes <- obs(meta_obj, data = 'diff_table', value = 'log2_median_ratio')
per_taxon_max_change <- unlist(lapply(per_taxon_fold_changes, function(tax_changes) max(abs(tax_changes))))
meta_obj_simp <- filter_taxa(meta_obj, per_taxon_max_change !=0, supertaxa = TRUE, reassign_obs = c(diff_table = FALSE))

heat_tree(meta_obj_simp,
          node_size = n_obs, # n_obs is a function that calculates, in this case, the number of ASV per taxon
          node_label = ifelse(is_ambiguous(taxon_names), "", taxon_names),
          node_color = log2_median_ratio, # A column from `obj$data$diff_table`
          node_color_range = diverging_palette(), # The built-in palette for diverging data                 
          node_color_interval = c(-10, 10), # The range of `log2_median_ratio` to display
          edge_color_interval = c(-10,10), # The range of `log2_median_ratio` to display
          node_size_axis_label = "Number of ASVs",
          node_label_size_range = c(0,0.03),
          node_color_axis_label = "Log2 ratio median proportions",
          repel_labels = TRUE,
          initial_layout = "re", layout = "da",
          make_node_legend=T,
          #node_label_max = 150,
          #overlap_avoidance=5,
          output_file = here("Results","Figures","Sed_D4_D66_SimpLegDiffTree.pdf")) # Saves the plot as a pdf file

