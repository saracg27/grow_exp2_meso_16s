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

load(here("Rdata","ps_water_obj.RData"))


### DIFFERENTIAL ABUNDANCE ####
###____####

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
##DESeq2 #####

#Order agglomeration
#ps_glom_ord <- tax_glom(ps,taxrank="Order",NArm = F)
# 96 orders in 140 samples

ps_obj <- ps

# Filter out features (OTUs or taxa) that have more 
#than 90% zeros across all samples
# This is done by checking the number of zeros in each row
#of the OTU table and comparing it to 90% of the total
#number of columns (samples)
ps_obj <- prune_taxa(rowSums(otu_table(ps_obj) == 0) 
                     < ncol(otu_table(ps_obj)) * 0.9, ps_obj)

# Removes 55 orders if applied on order agglom
# Removes 707 Asvs, 83% of the dataset, keeps 143 taxa..

######Sample Type ####

####|Scirpus VS Unplanted ####
ps_obj_NP_Scirpus <-  subset_samples(ps_obj, Sample_type=="No_plant" | Sample_type =="Scirpus")
ps_obj_NP_Scirpus <- prune_taxa(taxa_sums(ps_obj_NP_Scirpus)>0, ps_obj_NP_Scirpus)


diagdds = phyloseq_to_deseq2(ps_obj_NP_Scirpus, ~ Sample_type)

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
# Take the fitType minimizing the median absolute residuals
# here the local fit! 

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


signif_taxa$Differential_abundance <- "More abundant in Scirpus"
for (i in 1:nrow(signif_taxa)){
if(signif_taxa$log2FoldChange[i]<0){
  signif_taxa$Differential_abundance[i] <- "More abundant in Unplanted"
}
}
ggplot(signif_taxa,aes(x=log2FoldChange,y=Family,fill=Differential_abundance,shape=Differential_abundance))+
  theme_bw()+
  facet_grid(Phylum~.,scale='free',space='free',switch='y')+
  theme(strip.text.y.left = element_text(angle = 0))+
  geom_vline(xintercept = 0)+
  scale_shape_manual(values=c(22,23))+
  geom_point(size=3)


ggsave(filename = here("Results/Figures/", "deseq_diffab_water_ScirpusVSnoplant.png"),height = 7.5, width = 10.5, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_diffab_water_ScirpusVSnoplant.pdf"),height = 7.5, width = 10.5, dpi = 300)


# Convert p-values to -log10 scale
signif_taxa$neg_log10_pvalue <- -log10(signif_taxa$padj)

# Filter signif_taxa based on cutoffs
#filtered_signif_taxa <- signif_taxa[(abs(signif_taxa$log2FoldChange) > 1) & (signif_taxa$neg_log10_pvalue > 5), ]

# Create a new dataframe with the relevant columns
result_order_diff <- signif_taxa[, c("Order", "log2FoldChange", "neg_log10_pvalue")]
result_order_diff$Order

volcano_water_ScVSNp <- EnhancedVolcano(signif_taxa, 
                                          lab = signif_taxa$Family, 
                                          selectLab = signif_taxa$Family,
                                          x = 'log2FoldChange', 
                                          y = "padj",
                                          title = "Differential abundance in water",
                                          subtitle = "Scirpus vs Unplanted - Family",
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
volcano_water_ScVSNp

ggsave(filename = here("Results/Figures/", "deseq_volcano_water_ScirpusVSnoplant.png"),height = 7.5, width = 10.5, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_volcano_water_ScirpusVSnoplant.pdf"),height = 7.5, width = 10.5, dpi = 300)

######|Scirpus VS Triglochin ####
ps_obj_Triglo_Scirpus <-  subset_samples(ps_obj, Sample_type=="Triglochin" | Sample_type =="Scirpus")
ps_obj_Triglo_Scirpus <- prune_taxa(taxa_sums(ps_obj_Triglo_Scirpus)>0, ps_obj_Triglo_Scirpus)
ps_obj_Triglo_Scirpus <- prune_taxa(rowSums(otu_table(ps_obj_Triglo_Scirpus) == 0) 
           < ncol(otu_table(ps_obj_Triglo_Scirpus)) * 0.9, ps_obj_Triglo_Scirpus)

a <- data.frame(rowSums(otu_table(ps_obj_Triglo_Scirpus) == 0))
      
      
      
      
diagdds = phyloseq_to_deseq2(ps_obj_Triglo_Scirpus, ~ Sample_type)

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
# Take the fitType minimizing the median absolute residuals
# here the local fit! 

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


signif_taxa$Differential_abundance <- "More abundant in Triglochin"
for (i in 1:nrow(signif_taxa)){
  if(signif_taxa$log2FoldChange[i]<0){
    signif_taxa$Differential_abundance[i] <- "More abundant in Scirpus"
  }
}
ggplot(signif_taxa,aes(x=log2FoldChange,y=Family,fill=Differential_abundance,shape=Differential_abundance))+
  theme_bw()+
  facet_grid(Phylum~.,scale='free',space='free',switch='y')+
  theme(strip.text.y.left = element_text(angle = 0))+
  geom_vline(xintercept = 0)+
  scale_shape_manual(values=c(22,23))+
  geom_point(size=3)


ggsave(filename = here("Results/Figures/", "deseq_diffab_water_TrigloVSScirpus.png"),height = 7.5, width = 10.5, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_diffab_water_TrigloVSScirpus.pdf"),height = 7.5, width = 10.5, dpi = 300)



# Convert p-values to -log10 scale
signif_taxa$neg_log10_pvalue <- -log10(signif_taxa$padj)

# Filter signif_taxa based on cutoffs
#filtered_signif_taxa <- signif_taxa[(abs(signif_taxa$log2FoldChange) > 1) & (signif_taxa$neg_log10_pvalue > 5), ]

# Create a new dataframe with the relevant columns
result_order_diff <- signif_taxa[, c("Order", "log2FoldChange", "neg_log10_pvalue")]
result_order_diff$Order

volcano_water_TgVSSc <- EnhancedVolcano(signif_taxa, 
                                        lab = signif_taxa$Family, 
                                        selectLab = signif_taxa$Family,
                                        x = 'log2FoldChange', 
                                        y = "padj",
                                        title = "Differential abundance in water",
                                        subtitle = "Triglochin vs Scirpus - Family",
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
volcano_water_TgVSSc

# 
ggsave(filename = here("Results/Figures/", "deseq_volcano_water_TrigloVSScirpus.png"),height = 7.5, width = 10.5, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_volcano_water_TrigloVSScirpus.pdf"),height = 7.5, width = 10.5, dpi = 300)


# # Check LKM11
# LKM11 <- subset_taxa(ps_obj_Triglo_Scirpus,Family=="LKM11")
# 
# melt_LKM11 <- psmelt(LKM11)
# 
# melt_LKM11%>%
#   group_by(Sample_type)%>%
#   summarise(mean=mean(Abundance),
#             max=max(Abundance),
#             med=median(Abundance),
#             presence=sum(Abundance!=0))
# The LKM11 outlier has 15000 reads in a single sample 

# Other NA outlier is just present in one Scirpus sample


######|Scirpus VS Unplanted ####
ps_obj_NP_Triglochin <-  subset_samples(ps_obj, Sample_type=="No_plant" | Sample_type =="Triglochin")
ps_obj_NP_Triglochin <- prune_taxa(taxa_sums(ps_obj_NP_Triglochin)>0, ps_obj_NP_Triglochin)
ps_obj_NP_Triglochin <- prune_taxa(rowSums(otu_table(ps_obj_NP_Triglochin) == 0) 
                                    < ncol(otu_table(ps_obj_NP_Triglochin)) * 0.9, ps_obj_NP_Triglochin)


diagdds = phyloseq_to_deseq2(ps_obj_NP_Triglochin, ~ Sample_type)

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
# Take the fitType minimizing the median absolute residuals
# here the local fit! 

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


signif_taxa$Differential_abundance <- "More abundant in Triglochin"
for (i in 1:nrow(signif_taxa)){
  if(signif_taxa$log2FoldChange[i]<0){
    signif_taxa$Differential_abundance[i] <- "More abundant in Unplanted"
  }
}
ggplot(signif_taxa,aes(x=log2FoldChange,y=Family,fill=Differential_abundance,shape=Differential_abundance))+
  theme_bw()+
  facet_grid(Phylum~.,scale='free',space='free',switch='y')+
  theme(strip.text.y.left = element_text(angle = 0))+
  geom_vline(xintercept = 0)+
  scale_shape_manual(values=c(22,23))+
  geom_point(size=3)


ggsave(filename = here("Results/Figures/", "deseq_diffab_water_TriglochinVSnoplant.png"),height = 7.5, width = 10.5, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_diffab_water_TriglochinVSnoplant.pdf"),height = 7.5, width = 10.5, dpi = 300)



# Convert p-values to -log10 scale
signif_taxa$neg_log10_pvalue <- -log10(signif_taxa$padj)

# Filter signif_taxa based on cutoffs
#filtered_signif_taxa <- signif_taxa[(abs(signif_taxa$log2FoldChange) > 1) & (signif_taxa$neg_log10_pvalue > 5), ]

# Create a new dataframe with the relevant columns
result_order_diff <- signif_taxa[, c("Order", "log2FoldChange", "neg_log10_pvalue")]
result_order_diff$Order

volcano_water_TgVSNp <- EnhancedVolcano(signif_taxa, 
                                        lab = signif_taxa$Family, 
                                        selectLab = signif_taxa$Family,
                                        x = 'log2FoldChange', 
                                        y = "padj",
                                        title = "Differential abundance in water",
                                        subtitle = "Triglochin vs Unplanted - Family",
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
volcano_water_TgVSNp

ggsave(filename = here("Results/Figures/", "deseq_volcano_water_TriglochinVSnoplant.png"),height = 7.5, width = 10.5, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_volcano_water_TriglochinVSnoplant.pdf"),height = 7.5, width = 10.5, dpi = 300)








####Temperature  ####
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
# 
# 
#  diagdds2 = DESeq(diagdds,
#                   test = "Wald",
#                   fitType="parametric")
#  plotDispEsts(diagdds2)
#   
#   dispFit2 <-  diagdds2@rowRanges@elementMetadata@listData$dispFit
#   dispGeneEst2 <- diagdds2@rowRanges@elementMetadata@listData$dispGeneEst
#   test2 <- log(dispGeneEst2)-log(dispFit2)
# med_abs_res2 <- median(abs(test2),na.rm=TRUE)
# 
# Take the fitType minimizing the median absolute residuals?? 
# Here they are similar
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

res2

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

ggsave(filename = here("Results/Figures/", "deseq_diffab_water_temp.png"),height = 7.5, width = 10.5, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_diffab_water_temp.pdf"),height = 7.5, width = 10.5, dpi = 300)


# Convert p-values to -log10 scale
signif_taxa$neg_log10_pvalue <- -log10(signif_taxa$padj)

# Filter signif_taxa based on cutoffs
#filtered_signif_taxa <- signif_taxa[(abs(signif_taxa$log2FoldChange) > 1) & (signif_taxa$neg_log10_pvalue > 5), ]

# Create a new dataframe with the relevant columns
result_order_diff <- signif_taxa[, c("Order", "log2FoldChange", "neg_log10_pvalue")]
result_order_diff$Order

volcano_water_temp <- EnhancedVolcano(signif_taxa, 
                                     lab = signif_taxa$Family, 
                                     selectLab = signif_taxa$Family,
                                     x = 'log2FoldChange', 
                                     y = "padj",
                                     title = "Differential abundance in water",
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

volcano_water_temp

ggsave(filename = here("Results/Figures/", "deseq_volcano_water_temp.png"), plot = volcano_water_temp, width = 6, height = 6, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_volcano_roots_water.pdf"), plot = volcano_water_temp, width = 6, height = 6, dpi = 300)

### sanity check #
ps_warm<- subset_samples(ps,Temperature=="20C day - 10C night")
ps_warm <- prune_taxa(taxa_sums(ps_warm)>0, ps_warm)
sample_data(ps_warm)

ps_cold<- subset_samples(ps,Temperature=="10C day - 5C night")
ps_cold <- prune_taxa(taxa_sums(ps_cold)>0, ps_cold)


ps_Conthreep_cold <- subset_taxa(ps_cold,Order=="Conthreep")
ps_Conthreep_cold <- prune_taxa(taxa_sums(ps_Conthreep_cold)>0, ps_Conthreep_cold)

ps_Conthreep_warm <- subset_taxa(ps_warm,Order=="Conthreep")
ps_Conthreep_warm <- prune_taxa(taxa_sums(ps_Conthreep_warm)>0, ps_Conthreep_warm)

Ab_Conthreep_cold<- sum(rowSums(otu_table(ps_Conthreep_cold))) # 14525 reads of Conthreepbulbales in cold T°
Ab_Conthreep_warm <- sum(rowSums(otu_table(ps_Conthreep_warm))) # 6450 reads of Conthreepbulbales in warm T°




####Time  ####
#####|D0-D62 ####
#Subset ps_rhizo to keep samples from D0 and D62 sampling dates
ps_obj_d0d62 <-  subset_samples(ps, Time =="D0" | Time =="D62")
ps_obj_d0d62 <- prune_taxa(taxa_sums(ps_obj_d0d62)>0, ps_obj_d0d62)



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



