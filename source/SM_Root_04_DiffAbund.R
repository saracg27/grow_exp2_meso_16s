### Packages#####
library(here)
source(here("source", "libraries.R"))
library(phyloseq)
library(viridis)
library(metacoder)
library(dplyr)
library(ggpp)
library(EnhancedVolcano)


#### Data ####

load(here("Rdata","ps_obj.RData"))

#### Separation by temperature #####

ps_root_Warm <- subset_samples(ps, Sample.type =="Roots"& Temperature =="20°C day 10°C night")
ps_root_Warm <- prune_taxa(taxa_sums(ps_root_Warm)>0, ps_root_Warm)
ps_root_Warm # 42 samples - 2050  taxa
# check to see if you kept only Roots samples in WARM condition
str(sample_data(ps_root_Warm)) 

ps_root_Cold <- subset_samples(ps, Sample.type =="Roots"& Temperature =="10°C day 5°C night")
ps_root_Cold <- prune_taxa(taxa_sums(ps_root_Cold)>0, ps_root_Cold)
ps_root_Cold # 39 samples - 1587 taxa
# check to see if you kept only Roots samples in COLD condition
str(sample_data(ps_root_Cold))  


##DESeq2 #####

# Taxonomic agglomeration

## Warm
#ps_glom_ord_Warm<- tax_glom(ps_root_Warm,taxrank="Order",NArm = F)
ps_glom_gen_Warm<- tax_glom(ps_root_Warm,taxrank="Genus",NArm = F)


## Cold
#ps_glom_ord_Cold <- tax_glom(ps_root_Cold,taxrank="Order",NArm = F)
ps_glom_gen_Cold <- tax_glom(ps_root_Cold,taxrank="Genus",NArm = F)


#### SAMPLE TYPE ####

#####  Warm #####

ps_obj <- ps_glom_gen_Warm


####|Triglochin VS Scirpus ####

ps_obj_Scirpus_Triglochin <-  subset_samples(ps_obj, Plant.type=="Scirpus" | Plant.type =="Triglochin")
ps_obj_Scirpus_Triglochin <- prune_taxa(taxa_sums(ps_obj_Scirpus_Triglochin)>0, ps_obj_Scirpus_Triglochin)
ps_obj_Scirpus_Triglochin <- prune_taxa(rowSums(otu_table(ps_obj_Scirpus_Triglochin) == 0) 
                                        < ncol(otu_table(ps_obj_Scirpus_Triglochin)) * 0.9, ps_obj_Scirpus_Triglochin)

diagdds = phyloseq_to_deseq2(ps_obj_Scirpus_Triglochin, ~ Plant.type)

diagdds = estimateSizeFactors(diagdds, type = "poscounts")

diagdds = DESeq(diagdds,
                test = "Wald",
                fitType="local")
plotDispEsts(diagdds) # local has a better fit

res = results(diagdds, alpha=0.05) # results of the test FDR accepted = 5% 
plotMA(res)


res2 = res[order(res$padj, na.last=NA), ] # ordered and removes the ASVs which have padj = NAs 
plotMA(res2)
summary(res2)


# Volcano plots
res2 <- as.data.frame(res2)
tax <- data.frame(tax_table(ps_obj))
res2_taxa <- merge(res2,tax,by='row.names')


res2_taxa$neg_log10_pvalue <- -log10(res2_taxa$padj)

# Filter signif_taxa baroot on cutoffs
filtered_signif_taxa <- res2_taxa[(abs(res2_taxa$log2FoldChange) > 1) & (res2_taxa$neg_log10_pvalue > 5), ]


result_diff <- filtered_signif_taxa[, c("Genus", "log2FoldChange", "neg_log10_pvalue")]
result_diff$Genus

volcano_root_Scirpus_Triglochin_Warm <- EnhancedVolcano(res2_taxa, 
                                                         lab = res2_taxa$Genus, 
                                                         selectLab = res2_taxa$Genus,
                                                         x = 'log2FoldChange', 
                                                         y = "padj",
                                                         pCutoff = 1e-05,
                                                         title = "Differential ab. of 16S reads in roots",
                                                         titleLabSize = 15,
                                                         subtitle = "Triglochin vs Scirpus - Genus level",
                                                         pointSize = 3,
                                                         labSize = 3.5,
                                                         labFace = "italic",
                                                         colAlpha = 0.5,
                                                         legendPosition = 'bottom',
                                                         legendLabSize = 10,
                                                         legendIconSize = 4.0,
                                                         drawConnectors = TRUE,
                                                         widthConnectors = 0.6,
                                                         max.overlaps = 20
)

volcano_root_Scirpus_Triglochin_Warm

#ggsave(filename = here("Results_W&C/Figures/", "Volcano_root_Warm_TriglochinVSScirpus.png"),height = 7.5, width = 10.5, dpi = 300)
#ggsave(filename = here("Results_W&C/Figures/", "Volcano_root_Warm_TriglochinVSScirpus.png"),height = 7.5, width = 10.5, dpi = 300)


#####  Cold #####

ps_obj <- ps_glom_gen_Cold


####|Triglochin VS Scirpus ####

ps_obj_Scirpus_Triglochin <-  subset_samples(ps_obj, Plant.type=="Scirpus" | Plant.type =="Triglochin")
ps_obj_Scirpus_Triglochin <- prune_taxa(taxa_sums(ps_obj_Scirpus_Triglochin)>0, ps_obj_Scirpus_Triglochin)
ps_obj_Scirpus_Triglochin <- prune_taxa(rowSums(otu_table(ps_obj_Scirpus_Triglochin) == 0) 
                                        < ncol(otu_table(ps_obj_Scirpus_Triglochin)) * 0.9, ps_obj_Scirpus_Triglochin)

diagdds = phyloseq_to_deseq2(ps_obj_Scirpus_Triglochin, ~ Plant.type)

diagdds = estimateSizeFactors(diagdds, type = "poscounts")

diagdds = DESeq(diagdds,
                test = "Wald",
                fitType="local")
plotDispEsts(diagdds) # local has a better fit

res = results(diagdds, alpha=0.05) # results of the test FDR accepted = 5% 
plotMA(res)


res2 = res[order(res$padj, na.last=NA), ] # ordered and removes the ASVs which have padj = NAs 
plotMA(res2)
summary(res2)


# Volcano plots
res2 <- as.data.frame(res2)
tax <- data.frame(tax_table(ps_obj))
res2_taxa <- merge(res2,tax,by='row.names')


res2_taxa$neg_log10_pvalue <- -log10(res2_taxa$padj)

# Filter signif_taxa baroot on cutoffs
filtered_signif_taxa <- res2_taxa[(abs(res2_taxa$log2FoldChange) > 1) & (res2_taxa$neg_log10_pvalue > 5), ]


result_diff <- filtered_signif_taxa[, c("Genus", "log2FoldChange", "neg_log10_pvalue")]
result_diff$Genus

volcano_root_Scirpus_Triglochin_Cold <- EnhancedVolcano(res2_taxa, 
                                                         lab = res2_taxa$Genus, 
                                                         selectLab = res2_taxa$Genus,
                                                         x = 'log2FoldChange', 
                                                         y = "padj",
                                                         pCutoff = 1e-05,
                                                         title = "Differential ab. of 16S reads in roots",
                                                         titleLabSize = 15,
                                                         subtitle = "Triglochin vs Scirpus - Genus level",
                                                         pointSize = 3,
                                                         labSize = 3.5,
                                                         labFace = "italic",
                                                         colAlpha = 0.5,
                                                         legendPosition = 'bottom',
                                                         legendLabSize = 10,
                                                         legendIconSize = 4.0,
                                                         drawConnectors = TRUE,
                                                         widthConnectors = 0.6,
                                                         max.overlaps = 20
)

volcano_root_Scirpus_Triglochin_Cold

#ggsave(filename = here("Results_W&C/Figures/", "Volcano_root_Cold_TriglochinVSScirpus.png"),height = 7.5, width = 10.5, dpi = 300)
#ggsave(filename = here("Results_W&C/Figures/", "Volcano_root_Cold_TriglochinVSScirpus.png"),height = 7.5, width = 10.5, dpi = 300)


#####Common plots #####


Triglo_Scirpus <- ggarrange(volcano_root_Scirpus_Triglochin_Warm,volcano_root_Scirpus_Triglochin_Cold,
                            ncol=2,
                            labels=c("Warm","Cold"))


ggsave(filename = here("Results_W&C/Figures/", "Volcano_Root_PlantTypes.pdf"),
       plot=Triglo_Scirpus, height = 7.5, width = 10.5, dpi = 300)

ggsave(filename = here("Results_W&C/Figures/", "Volcano_Root_PlantTypes.png"),
       plot=Triglo_Scirpus, height = 7.5, width = 10.5, dpi = 300)



####~~~~~~~~~~~~~~~~####

#### TIME ???  ####

####~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~####
########~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~####
########~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~####
########~~~~~~~~~~~~~~~~####
#### Not separated by temp ####
####~~~~~~~~~~~~~~~~####



### DIFFERENTIAL ABUNDANCE ####
#### Data ####

load(here("Rdata","ps_obj.RData"))

ps_root <-  subset_samples(ps, Sample.type =="Roots")
ps_root <- prune_taxa(taxa_sums(ps_root)>0, ps_root)
ps_root
# 81 samples 5954 taxa

#sample_names(ps_root) <- sub("Meso.*","",sample_names(ps_root))

## ______ #####

### ....Metacoder Object####
# Order agglomeration
ps_glom_ord <- tax_glom(ps_root,taxrank="Order")
#ps_glom_fam <- tax_glom(ps_root,taxrank="Family",NArm = F)

ps_obj <- ps_glom_ord


meta_obj <- parse_phyloseq(ps_obj) 

B <- data.frame(tax_table(ps_obj))
B$Order <- as.factor(B$Order)


# transforms the phyloseq object in a metacoder object


####******  Plant Type ####

meta_obj$data$otu_table <- calc_obs_props(meta_obj, data = "otu_table", 
                                          cols = meta_obj$data$sample_data$sample_id)


meta_obj$data$tax_table <- calc_taxon_abund(meta_obj, data = "otu_table",
                                            cols = meta_obj$data$sample_data$sample_id,
                                            group= rep("total_count", nrow(meta_obj$data$sample_data)))
                                                                                                  
meta_obj$data$tax_table <- subset(meta_obj$data$tax_table,
                                  meta_obj$data$tax_table$taxon_id %in%meta_obj$data$otu_table$taxon_id)
# Quantifies each taxon abundance -> sums the ASVs abundance with the same taxonomy across the whole dataset


meta_obj$data$diff_table <- compare_groups(meta_obj, data = "tax_table",
                                           cols = meta_obj$data$sample_data$sample_id,
                                           groups = meta_obj$data$sample_data$Plant.type)# What columns of sample data to use
# Compares the taxonomic difference based on ASV abundance

# Adjust the p_value with FDR correction
meta_obj$data$diff_table$wilcox_FDR_p_value <- p.adjust(meta_obj$data$diff_table$wilcox_p_value,
                                                        method = "fdr")
# See the range of FDR p_values
range(meta_obj$data$diff_table$wilcox_FDR_p_value, finite = TRUE) 

meta_obj$data$diff_table$log2_median_ratio[meta_obj$data$diff_table$wilcox_FDR_p_value > 0.05] <- 0
# Sets the log2MedianRatio to 0 if the corrected p values are bigger than 0.05

# 
# per_taxon_fold_changes <- obs(meta_obj, data = 'diff_table', value = 'log2_median_ratio')
# per_taxon_max_change <- unlist(lapply(per_taxon_fold_changes, function(tax_changes) max(abs(tax_changes))))
# meta_obj_simp <- filter_taxa(meta_obj, per_taxon_max_change !=0, supertaxa = TRUE, reassign_obs = c(diff_table = FALSE))


heat_tree(meta_obj,
          node_size = n_obs , # n_obs is a function that calculates, in this case, the number of ASV per taxon
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
          output_file = here("Results","Figures","Roots_Plants_DiffTree.pdf")) # Saves the plot as a pdf file





### Sanity check of metacoder results ###

Diff_table <- as.data.frame(meta_obj$data$diff_table)
tax <- meta_obj$data$tax_data
tax_order <- tax[,c("taxon_id","Order")] 
Diff_table_tax <- merge(tax_order,Diff_table,'taxon_id')


### According to the tree,
# Cyanobacterials should be more abundant in Triglochin than in Scirpus

ps_scirpus<- subset_samples(ps_root,Plant.type=="Scirpus")
ps_scirpus <- prune_taxa(taxa_sums(ps_scirpus)>0, ps_scirpus)

ps_cyano_sci <- subset_taxa(ps_scirpus,Order=="Cyanobacteriales")
ps_cyano_sci <- prune_taxa(taxa_sums(ps_cyano_sci)>0, ps_cyano_sci)

Ab_cyano_scirpus <- sum(rowSums(otu_table(ps_cyano_sci))) # 3398 reads of Cyanobacteriales in Scirpus


ps_triglo<- subset_samples(ps_root,Plant.type=="Triglochin")
ps_triglo <- prune_taxa(taxa_sums(ps_triglo)>0, ps_triglo)

ps_cyano_tri <- subset_taxa(ps_triglo,Order=="Cyanobacteriales")
ps_cyano_tri <- prune_taxa(taxa_sums(ps_cyano_tri)>0, ps_cyano_tri)

Ab_cyano_tri <- sum(rowSums(otu_table(ps_cyano_tri))) # 19416 reads of Cyanobacteriales in Triglochin


### Micavibrionales
ps_mica_sci <- subset_taxa(ps_scirpus,Order=="Micavibrionales")
ps_mica_sci <- prune_taxa(taxa_sums(ps_mica_sci)>0, ps_mica_sci)

Ab_mica_scirpus <- sum(rowSums(otu_table(ps_mica_sci))) # 3398 reads of micabacteriales in Scirpus

ps_mica_tri <- subset_taxa(ps_triglo,Order=="Micavibrionales")
ps_mica_tri <- prune_taxa(taxa_sums(ps_mica_tri)>0, ps_mica_tri)

Ab_mica_tri <- sum(rowSums(otu_table(ps_mica_tri))) # 19416 reads of micabacteriales in Triglochin



### Temperature  ####

ps_glom_ord <- tax_glom(ps_root,taxrank="Order",NArm = F)
ps.relab <- transform_sample_counts(ps_glom_ord,function(x) x/sum(x))

meta_obj <- parse_phyloseq(ps.relab) 

meta_obj$data$tax_abund <- calc_taxon_abund(meta_obj, "otu_table")


meta_obj$data$diff_table <- compare_groups(meta_obj, data = "tax_abund",
                                           cols = meta_obj$data$sample_data$sample_id,
                                           groups = meta_obj$data$sample_data$Temperature)# What columns of sample data to use
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


signif_taxa <- subset(meta_obj$data$diff_table,meta_obj$data$diff_table$wilcox_FDR_p_value <0.05)
signif_taxa_taxo <- subset(meta_obj$data$tax_data,meta_obj$data$tax_data$taxon_id%in%signif_taxa$taxon_id)



###____####
#####DESeq2 #####

#Order agglomeration
ps_glom_ord <- tax_glom(ps_root,taxrank="Order",NArm = F)

ps_obj <- ps_glom_ord

# Filter out features (OTUs or taxa) that have more 
#than 90% zeros across all samples
# This is done by checking the number of zeros in each row
#of the OTU table and comparing it to 90% of the total
#number of columns (samples)
ps_obj <- prune_taxa(rowSums(otu_table(ps_obj) == 0) 
                     < ncol(otu_table(ps_obj)) * 0.9, ps_obj)

# Removes 99 orders 


####**Plant Type ####
diagdds = phyloseq_to_deseq2(ps_obj, ~ Plant.type)

# calculate geometric means prior to estimate size factors
# gm_mean = function(x, na.rm=TRUE){
#   exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
# }
# 
# 
# geoMeans = apply(counts(diagdds), 1, gm_mean)
# geoMeansvis <- as.data.frame(geoMeans)
# 
# ?estimateSizeFactors
# diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

diagdds = estimateSizeFactors(diagdds, type = "poscounts")
# For each ASV -> calculates its geometric mean. Doesn't take into account the zeros
# exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
# https://support.bioconductor.org/p/62246/#62250


diagdds = DESeq(diagdds,
                test = "Wald",
                fitType="local")
plotDispEsts(diagdds)

# Which fit to use ?
# https://support.bioconductor.org/p/81094/

#  dispFit <-  diagdds@rowRanges@elementMetadata@listData$dispFit
# dispGeneEst <- diagdds@rowRanges@elementMetadata@listData$dispGeneEst
# test <- log(dispGeneEst)-log(dispFit)
#  med_abs_res <- median(abs(test),na.rm=TRUE)
# ## 
# ## 
# diagdds2 = DESeq(diagdds,
#                  test = "Wald",
#                  fitType="parametric")
# plotDispEsts(diagdds2)
# # 
#  dispFit2 <-  diagdds2@rowRanges@elementMetadata@listData$dispFit
#  dispGeneEst2 <- diagdds2@rowRanges@elementMetadata@listData$dispGeneEst
#  test2 <- log(dispGeneEst2)-log(dispFit2)
#  med_abs_res2 <- median(abs(test2),na.rm=TRUE)
# 
# Take the fitType minimizing the median absolute residuals?? 
# HERE THE LOCAL ! 

res = results(diagdds, alpha=0.05) # results of the test FDR accepted = 5% 
plotMA(res)

res2 = res[order(res$padj, na.last=NA), ] # ordered and removes the ASVs which have padj = NAs 
plotMA(res2)

summary(res2)
sigtab = res2[(res2$padj < 0.05), ] # Select significant padj
sigtab <- as.data.frame(sigtab)
sigtab <- cbind(rownames(sigtab),sigtab)
colnames(sigtab)[1] <- "OTU"

tax <- data.frame(tax_table(ps_obj))
tax_diffab <- subset(tax,rownames(tax)%in%rownames(sigtab))

signif_taxa <- merge(tax_diffab,sigtab,by='row.names')
res2

#### ¨¨¨¨¨¨Simple plot #####
signif_taxa$Differential_abundance <- "More abundant in Triglochin"
for (i in 1:nrow(signif_taxa)){
if(signif_taxa$log2FoldChange[i]<0){
  signif_taxa$Differential_abundance[i] <- "More abundant in Scirpus"
}
}
ggplot(signif_taxa,aes(x=log2FoldChange,y=Order,fill=Differential_abundance,shape=Differential_abundance))+
  theme_bw()+
  facet_grid(Phylum~.,scale='free',space='free',switch='y')+
  theme(strip.text.y.left = element_text(angle = 0))+
  geom_vline(xintercept = 0)+
  scale_shape_manual(values=c(22,23))+
  geom_point(size=3)


ggsave(filename = here("Results/Figures/", "deseq_diffab_roots_plant.png"),height = 7.5, width = 10.5, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_diffab_roots_plant.pdf"),height = 7.5, width = 10.5, dpi = 300)


#### ¨¨¨¨¨¨Volcano plot #####
# Convert p-values to -log10 scale
signif_taxa$neg_log10_pvalue <- -log10(signif_taxa$padj)

# Filter signif_taxa based on cutoffs
#filtered_signif_taxa <- signif_taxa[(abs(signif_taxa$log2FoldChange) > 1) & (signif_taxa$neg_log10_pvalue > 5), ]

# Create a new dataframe with the relevant columns
result_order_diff <- signif_taxa[, c("Order", "log2FoldChange", "neg_log10_pvalue")]
result_order_diff$Order

volcano_ro_plants <- EnhancedVolcano(signif_taxa, 
                                          lab = signif_taxa$Order, 
                                          selectLab = signif_taxa$Order,
                                          x = 'log2FoldChange', 
                                          y = "padj",
                                          title = "Differential abundance in plant roots",
                                          subtitle = "Triglochin vs Scirpus - Order",
                                          pointSize = 4.0,
                                          labSize = 3.0,
                                          colAlpha = 0.5,
                                          legendPosition = 'top',
                                          legendLabSize = 10,
                                          legendIconSize = 4.0,
                                          drawConnectors = TRUE,
                                          widthConnectors = 0.6,
                                          max.overlaps = 20
)
volcano_ro_plants

ggsave(filename = here("Results/Figures/", "deseq_volcano_roots_plant.png"),height = 7.5, width = 10.5, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_volcano_roots_plant.pdf"),height = 7.5, width = 10.5, dpi = 300)



####**Temperature  ####
diagdds = phyloseq_to_deseq2(ps_obj, ~ Temperature)

# calculate geometric means prior to estimate size factors
# gm_mean = function(x, na.rm=TRUE){
#   exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
# }
# 
# 
# geoMeans = apply(counts(diagdds), 1, gm_mean)
# geoMeansvis <- as.data.frame(geoMeans)
# 
# ?estimateSizeFactors
# diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

diagdds = estimateSizeFactors(diagdds, type = "poscounts")
# For each ASV -> calculates its geometric mean. Doesn't take into account the zeros
# exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
# https://support.bioconductor.org/p/62246/#62250


diagdds = DESeq(diagdds,
                test = "Wald",
                fitType="local")
plotDispEsts(diagdds)

# Which fit to use ?
# https://support.bioconductor.org/p/81094/

#  dispFit <-  diagdds@rowRanges@elementMetadata@listData$dispFit
# dispGeneEst <- diagdds@rowRanges@elementMetadata@listData$dispGeneEst
# test <- log(dispGeneEst)-log(dispFit)
#  med_abs_res <- median(abs(test),na.rm=TRUE)
# ## 
# ## 
# diagdds2 = DESeq(diagdds,
#                  test = "Wald",
#                  fitType="parametric")
# plotDispEsts(diagdds2)
# # 
#  dispFit2 <-  diagdds2@rowRanges@elementMetadata@listData$dispFit
#  dispGeneEst2 <- diagdds2@rowRanges@elementMetadata@listData$dispGeneEst
#  test2 <- log(dispGeneEst2)-log(dispFit2)
#  med_abs_res2 <- median(abs(test2),na.rm=TRUE)
# 
# Take the fitType minimizing the median absolute residuals?? 
# HERE THE LOCAL ! 
results(diagdds)
res = results(diagdds, alpha=0.05) # results of the test FDR accepted = 5% 
plotMA(res)

res2 = res[order(res$padj, na.last=NA), ] # ordered and removes the ASVs which have padj = NAs 
plotMA(res2)

summary(res2)
sigtab = res2[(res2$padj < 0.05), ] # Select significant padj
sigtab <- as.data.frame(sigtab)
sigtab <- cbind(rownames(sigtab),sigtab)
colnames(sigtab)[1] <- "OTU"

tax <- data.frame(tax_table(ps_obj))
tax_diffab <- subset(tax,rownames(tax)%in%rownames(sigtab))

signif_taxa <- merge(tax_diffab,sigtab,by='row.names')

#### ¨¨¨¨¨¨Simple plot #####


signif_taxa$Temp <- "Warm"
for (i in 1:nrow(signif_taxa)){
  if(signif_taxa$log2FoldChange[i]<0){
    signif_taxa$Temp[i] <- "Cold"
  }
}
ggplot(signif_taxa,aes(x=log2FoldChange,y=Order,color=Temp))+
  theme_bw()+
  facet_grid(Phylum~.,scale='free',space='free',switch='y')+
  theme(strip.text.y.left = element_text(angle = 0))+
  geom_vline(xintercept = 0)+
  scale_color_manual(values=c("#0000FE","#B02223"),name="Temperature")+
  geom_point(size=3)

ggsave(filename = here("Results/Figures/", "deseq_diffab_roots_temp.png"),height = 7.5, width = 10.5, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_diffab_roots_temp.pdf"),height = 7.5, width = 10.5, dpi = 300)

#### ¨¨¨¨¨¨Volcano plot #####
# Convert p-values to -log10 scale
signif_taxa$neg_log10_pvalue <- -log10(signif_taxa$padj)

# Filter signif_taxa based on cutoffs
#filtered_signif_taxa <- signif_taxa[(abs(signif_taxa$log2FoldChange) > 1) & (signif_taxa$neg_log10_pvalue > 5), ]

# Create a new dataframe with the relevant columns
result_order_diff <- signif_taxa[, c("Order", "log2FoldChange", "neg_log10_pvalue")]
result_order_diff$Order

volcano_ro_temp <- EnhancedVolcano(signif_taxa, 
                                     lab = signif_taxa$Order, 
                                     selectLab = signif_taxa$Order,
                                     x = 'log2FoldChange', 
                                     y = "padj",
                                     title = "Differential abundance in roots",
                                     subtitle = "Warm vs Cold - Order",
                                     pointSize = 4.0,
                                     labSize = 3.0,
                                     colAlpha = 0.5,
                                     legendPosition = 'top',
                                     legendLabSize = 10,
                                     legendIconSize = 4.0,
                                     drawConnectors = TRUE,
                                     widthConnectors = 0.6,
                                     max.overlaps = 20
)

volcano_ro_temp

ggsave(filename = here("Results/Figures/", "deseq_volcano_roots_temp.png"), plot = volcano_ro_temp, width = 6, height = 6, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_volcano_roots_temp.pdf"), plot = volcano_ro_temp, width = 6, height = 6, dpi = 300)

### sanity check #
ps_warm<- subset_samples(ps_root,Temperature=="20°C day 10°C night")
ps_warm <- prune_taxa(taxa_sums(ps_warm)>0, ps_warm)

ps_cold<- subset_samples(ps_root,Temperature=="10°C day 5°C night")
ps_cold <- prune_taxa(taxa_sums(ps_cold)>0, ps_cold)


ps_desulfo_cold <- subset_taxa(ps_cold,Order=="Desulfobulbales")
ps_desulfo_cold <- prune_taxa(taxa_sums(ps_desulfo_cold)>0, ps_desulfo_cold)

ps_desulfo_warm <- subset_taxa(ps_warm,Order=="Desulfobulbales")
ps_desulfo_warm <- prune_taxa(taxa_sums(ps_desulfo_warm)>0, ps_desulfo_warm)

Ab_desulfo_cold<- sum(rowSums(otu_table(ps_desulfo_cold))) # 854 reads of Desulfobulbales in cold T°
Ab_desulfo_warm <- sum(rowSums(otu_table(ps_desulfo_warm))) # 3500 reads of Desulfobulbales in warm T°

####**Time D0-D62  ####
####**Time D0-D62  ####

#Subset ps_rhizo to keep samples from D0 and D62 sampling dates
ps_obj_d0d62 <-  subset_samples(ps_root, Time =="D0" | Time =="D62")
ps_obj_d0d62 <- prune_taxa(taxa_sums(ps_obj_d0d62)>0, ps_obj_d0d62)

# agglomerate at order level
ps_obj_d0d62 <- tax_glom(ps_obj_d0d62,taxrank="Order",NArm = F)

# remove taxa with not present in 90% of the samples
ps_obj <- ps_obj_d0d62
ps_obj <- prune_taxa(rowSums(otu_table(ps_obj) == 0) 
                     < ncol(otu_table(ps_obj)) * 0.9, ps_obj)
# Results in the loss of 68 orders

diagdds = phyloseq_to_deseq2(ps_obj_d0d62, ~ Time)
# calculate geometric means prior to estimate size factors
# gm_mean = function(x, na.rm=TRUE){
#   exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
# }
# 
# 
# geoMeans = apply(counts(diagdds), 1, gm_mean)
# geoMeansvis <- as.data.frame(geoMeans)
# 
# ?estimateSizeFactors
# diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

diagdds = estimateSizeFactors(diagdds, type = "poscounts")
# For each ASV -> calculates its geometric mean. Doesn't take into account the zeros
# exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
# https://support.bioconductor.org/p/62246/#62250


diagdds = DESeq(diagdds,
                test = "Wald",
                fitType="local")
plotDispEsts(diagdds)

# Which fit to use ?
# https://support.bioconductor.org/p/81094/

#  dispFit <-  diagdds@rowRanges@elementMetadata@listData$dispFit
# dispGeneEst <- diagdds@rowRanges@elementMetadata@listData$dispGeneEst
# test <- log(dispGeneEst)-log(dispFit)
#  med_abs_res <- median(abs(test),na.rm=TRUE)
# ## 
# ## 
# diagdds2 = DESeq(diagdds,
#                  test = "Wald",
#                  fitType="parametric")
# plotDispEsts(diagdds2)
# # 
#  dispFit2 <-  diagdds2@rowRanges@elementMetadata@listData$dispFit
#  dispGeneEst2 <- diagdds2@rowRanges@elementMetadata@listData$dispGeneEst
#  test2 <- log(dispGeneEst2)-log(dispFit2)
#  med_abs_res2 <- median(abs(test2),na.rm=TRUE)
# 
# Take the fitType minimizing the median absolute residuals?? 
# HERE THE LOCAL ! 
results(diagdds)
res = results(diagdds, alpha=0.05) # results of the test FDR accepted = 5% 
plotMA(res)

res2 = res[order(res$padj, na.last=NA), ] # ordered and removes the ASVs which have padj = NAs 
plotMA(res2)

summary(res2)
sigtab = res2[(res2$padj < 0.05), ] # Select significant padj
sigtab <- as.data.frame(sigtab)
sigtab <- cbind(rownames(sigtab),sigtab)
colnames(sigtab)[1] <- "OTU"

tax <- data.frame(tax_table(ps_obj))
tax_diffab <- subset(tax,rownames(tax)%in%rownames(sigtab))

signif_taxa <- merge(tax_diffab,sigtab,by='row.names')

#### ¨¨¨¨¨¨Simple plot #####


signif_taxa$Time <- "More abundant in D62"
for (i in 1:nrow(signif_taxa)){
  if(signif_taxa$log2FoldChange[i]<0){
    signif_taxa$Time[i] <- "More abundant in D0"
  }
}
ggplot(signif_taxa,aes(x=log2FoldChange,y=Order,color=Time))+
  theme_bw()+
  facet_grid(Phylum~.,scale='free',space='free',switch='y')+
  theme(strip.text.y.left = element_text(angle = 0))+
  geom_vline(xintercept = 0)+
  scale_color_manual(values=c("black","yellow"),name="Differential abundance")+
  geom_point(size=3)

ggsave(filename = here("Results/Figures/", "deseq_diffab_roots_time.png"),height = 7.5, width = 10.5, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_diffab_roots_time.pdf"),height = 7.5, width = 10.5, dpi = 300)

#### ¨¨¨¨¨¨Volcano plot #####
# Convert p-values to -log10 scale
signif_taxa$neg_log10_pvalue <- -log10(signif_taxa$padj)

# Filter signif_taxa based on cutoffs
#filtered_signif_taxa <- signif_taxa[(abs(signif_taxa$log2FoldChange) > 1) & (signif_taxa$neg_log10_pvalue > 5), ]

# Create a new dataframe with the relevant columns
result_order_diff <- signif_taxa[, c("Order", "log2FoldChange", "neg_log10_pvalue")]
result_order_diff$Order

volcano_ro_time <- EnhancedVolcano(signif_taxa, 
                                   lab = signif_taxa$Order, 
                                   selectLab = signif_taxa$Order,
                                   x = 'log2FoldChange', 
                                   y = "padj",
                                   title = "Differential abundance in roots",
                                   subtitle = "D62 vs D0 - Order",
                                   pointSize = 4.0,
                                   labSize = 3.0,
                                   colAlpha = 0.5,
                                   legendPosition = 'top',
                                   legendLabSize = 10,
                                   legendIconSize = 4.0,
                                   drawConnectors = TRUE,
                                   widthConnectors = 0.6,
                                   max.overlaps = 20
)

volcano_ro_time

ggsave(filename = here("Results/Figures/", "deseq_volcano_roots_time.png"), plot = volcano_ro_time,height = 7.5, width = 10.5, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_volcano_roots_time.pdf"), plot = volcano_ro_time,height = 7.5, width = 10.5, dpi = 300)


#### BAS #####
meta_obj <- parse_phyloseq(ps_obj) 

B <- data.frame(tax_table(ps_obj))
B$Order <- as.factor(B$Order)

A <- data.frame(otu_table(ps_obj))

C <- data.frame(sample_data(ps_obj))

head(meta_obj$data$otu_table)

meta_obj$data$tax_abund <- calc_taxon_abund(meta_obj, "otu_table",
                                            cols = meta_obj$data$sample_data$sample_id,
                                            groups = rep("total_count", nrow(meta_obj$data$sample_data)))

meta_obj$data$tax_abund$relab <- meta_obj$data$tax_abund$total_count/1202370


meta_obj$data$diff_table <- compare_groups(meta_obj, data = "otu_table",
                                           cols = meta_obj$data$sample_data$sample_id,
                                           groups = meta_obj$data$sample_data$Plant.type)# What columns of sample data to use


A <- meta_obj$data$diff_table
tax <- meta_obj$data$tax_data
tax_abund <- meta_obj$data$tax_abund

# Compares the taxonomic difference based on ASV abundance

# Adjust the p_value with FDR correction
meta_obj$data$diff_table$wilcox_FDR_p_value <- p.adjust(meta_obj$data$diff_table$wilcox_p_value,
                                                        method = "fdr")
# See the range of FDR p_values
range(meta_obj$data$diff_table$wilcox_FDR_p_value, finite = TRUE) 

meta_obj$data$diff_table$log2_median_ratio[meta_obj$data$diff_table$wilcox_FDR_p_value > 0.05] <- 0
# Sets the log2MedianRatio to 0 if the corrected p values are bigger than 0.05

A <- meta_obj$data$diff_table
Otu <- data.frame(otu_table(ps_obj))


# 
per_taxon_fold_changes <- obs(meta_obj, data = 'diff_table', value = 'log2_median_ratio')
per_taxon_max_change <- unlist(lapply(per_taxon_fold_changes, function(tax_changes) max(abs(tax_changes))))
meta_obj_simp <- filter_taxa(meta_obj, per_taxon_max_change !=0, supertaxa = TRUE, reassign_obs = c(diff_table = FALSE))

options(digits = 2)    
heat_tree(meta_obj,
          node_size = relab, # n_obs is a function that calculates, in this case, the number of ASV per taxon
          node_label = taxon_names,
          node_color = log2_median_ratio, # A column from `obj$data$diff_table`
          node_color_range = diverging_palette(), # The built-in palette for diverging data                 
          #node_color_interval = c(-10, 10), # The range of `log2_median_ratio` to display
          #node_edge_interval = c(-10, 10), # The range of `log2_median_ratio` to display
          #node_size_axis_label = "Relative abundance in dataset",
          #node_size_digits=2,
          #node_size_interval=c(0,1),
          node_color_axis_label = "Log2 ratio median proportions",
          repel_labels = TRUE,
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
          output_file = here("Results","Figures","Roots_Plants_DiffTree.pdf")) # Saves the plot as a pdf file

?heat_tree


otu <- as.data.frame(otu_table(ps_obj))
sum(rowSums(otu))
