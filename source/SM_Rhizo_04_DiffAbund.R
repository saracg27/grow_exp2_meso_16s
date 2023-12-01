### Packages#####
library(here)
source(here("source", "libraries.R"))
library(phyloseq)
library(viridis)
library(metacoder)
library(dplyr)
library(ggpp)
library(DESeq2)
library(EnhancedVolcano)

#### Data ####

load(here("Rdata","ps_obj.RData"))

#### Separation by temperature #####

ps_rhizo_Warm <- subset_samples(ps, Sample.type =="Rhizosphere"& Temperature =="20°C day 10°C night")
ps_rhizo_Warm <- prune_taxa(taxa_sums(ps_rhizo_Warm)>0, ps_rhizo_Warm)
ps_rhizo_Warm # 43 samples - 4630  taxa
# check to see if you kept only RHIZOSPHERE samples in WARM condition
str(sample_data(ps_rhizo_Warm)) 

ps_rhizo_Cold <- subset_samples(ps, Sample.type =="Rhizosphere"& Temperature =="10°C day 5°C night")
ps_rhizo_Cold <- prune_taxa(taxa_sums(ps_rhizo_Cold)>0, ps_rhizo_Cold)
ps_rhizo_Cold # 38 samples - 4522 taxa
# check to see if you kept only RHIZOSPHERE samples in COLD condition
str(sample_data(ps_rhizo_Cold)) 




##DESeq2 #####

# Taxonomic agglomeration

## Warm
#ps_glom_ord_Warm<- tax_glom(ps_rhizo_Warm,taxrank="Order",NArm = F)
ps_glom_gen_Warm<- tax_glom(ps_rhizo_Warm,taxrank="Genus",NArm = F)


## Cold
#ps_glom_ord_Cold <- tax_glom(ps_rhizo_Cold,taxrank="Order",NArm = F)
ps_glom_gen_Cold <- tax_glom(ps_rhizo_Cold,taxrank="Genus",NArm = F)


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

# Filter signif_taxa barhizo on cutoffs
filtered_signif_taxa <- res2_taxa[(abs(res2_taxa$log2FoldChange) > 1) & (res2_taxa$neg_log10_pvalue > 5), ]


result_diff <- filtered_signif_taxa[, c("Genus", "log2FoldChange", "neg_log10_pvalue")]
result_diff$Genus

volcano_rhizo_Scirpus_Triglochin_Warm <- EnhancedVolcano(res2_taxa, 
                                                       lab = res2_taxa$Genus, 
                                                       selectLab = res2_taxa$Genus,
                                                       x = 'log2FoldChange', 
                                                       y = "padj",
                                                       pCutoff = 1e-05,
                                                       title = "Differential ab. of 16S reads in rhizosphere",
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

volcano_rhizo_Scirpus_Triglochin_Warm

#ggsave(filename = here("Results_W&C/Figures/", "Volcano_rhizo_Warm_TriglochinVSScirpus.png"),height = 7.5, width = 10.5, dpi = 300)
#ggsave(filename = here("Results_W&C/Figures/", "Volcano_rhizo_Warm_TriglochinVSScirpus.png"),height = 7.5, width = 10.5, dpi = 300)


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

# Filter signif_taxa barhizo on cutoffs
filtered_signif_taxa <- res2_taxa[(abs(res2_taxa$log2FoldChange) > 1) & (res2_taxa$neg_log10_pvalue > 5), ]


result_diff <- filtered_signif_taxa[, c("Genus", "log2FoldChange", "neg_log10_pvalue")]
result_diff$Genus

volcano_rhizo_Scirpus_Triglochin_Cold <- EnhancedVolcano(res2_taxa, 
                                                       lab = res2_taxa$Genus, 
                                                       selectLab = res2_taxa$Genus,
                                                       x = 'log2FoldChange', 
                                                       y = "padj",
                                                       pCutoff = 1e-05,
                                                       title = "Differential ab. of 16S reads in rhizosphere",
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

volcano_rhizo_Scirpus_Triglochin_Cold

#ggsave(filename = here("Results_W&C/Figures/", "Volcano_rhizo_Cold_TriglochinVSScirpus.png"),height = 7.5, width = 10.5, dpi = 300)
#ggsave(filename = here("Results_W&C/Figures/", "Volcano_rhizo_Cold_TriglochinVSScirpus.png"),height = 7.5, width = 10.5, dpi = 300)


#####Common plots #####


Triglo_Scirpus <- ggarrange(volcano_rhizo_Scirpus_Triglochin_Warm,volcano_rhizo_Scirpus_Triglochin_Cold,
                            ncol=2,
                            labels=c("Warm","Cold"))


ggsave(filename = here("Results_W&C/Figures/", "Volcano_rhizo_PlantTypes.pdf"),
       plot=Triglo_Scirpus, height = 7.5, width = 10.5, dpi = 300)

ggsave(filename = here("Results_W&C/Figures/", "Volcano_rhizo_PlantTypes.png"),
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

ps_rhizo <-  subset_samples(ps, Sample.type =="Rhizosphere")
ps_rhizo <- prune_taxa(taxa_sums(ps_rhizo)>0, ps_rhizo)
ps_rhizo
# 81 samples 5954 taxa

###____####
####Metacoder on all sed samples####

### ....Metacoder Object####
# Order agglomeration
ps_glom_ord <- tax_glom(ps_rhizo,taxrank="Order",NArm = F)
ps_glom_fam <- tax_glom(ps_rhizo,taxrank="Family",NArm = F)

ps_obj <- ps_glom_ord
# Relative abundance transformation
ps.relab <- transform_sample_counts(ps_rhizo,function(x) x/sum(x))

meta_obj <- parse_phyloseq(ps.relab) 
# transforms the phyloseq object in a metacoder object


### Plant Type ###
meta_obj$data$tax_abund <- calc_taxon_abund(meta_obj, "otu_table")


meta_obj$data$diff_table <- compare_groups(meta_obj, data = "tax_abund",
                                           cols = meta_obj$data$sample_data$sample_id,
                                           groups = meta_obj$data$sample_data$Plant.type)# What columns of sample data to use
# Compares the taxonomic difference based on ASV abundance

# Adjust the p_value with FDR correction
meta_obj$data$diff_table$wilcox_FDR_p_value <- p.adjust(meta_obj$data$diff_table$wilcox_p_value,
                                                        method = "fdr")

# See the range of FDR p_values
range(meta_obj$data$diff_table$wilcox_FDR_p_value, finite = TRUE) 
# No significance at the ASV level
# No significance at the family level
# No significance at the order level




### Temperature  ####

ps_glom_ord <- tax_glom(ps_rhizo,taxrank="Order",NArm = F)
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
ps_glom_ord <- tax_glom(ps_rhizo,taxrank="Order",NArm = F)
#ps_glom_fam <- tax_glom(ps_rhizo,taxrank="Family",NArm = F)

ps_obj <- ps_glom_ord

# Filter out features (OTUs or taxa) that have more 
#than 90% zeros across all samples
# This is done by checking the number of zeros in each row
#of the OTU table and comparing it to 90% of the total
#number of columns (samples)
ps_obj <- prune_taxa(rowSums(otu_table(ps_obj) == 0) 
                              < ncol(otu_table(ps_obj)) * 0.9, ps_obj)

# Removes 73 orders 


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


ggsave(filename = here("Results/Figures/", "deseq_diffab_rhizo_plant.png"),height = 7.5, width = 10.5, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_diffab_rhizo_plant.pdf"),height = 7.5, width = 10.5, dpi = 300)


#### ¨¨¨¨¨¨Volcano plot #####
# Convert p-values to -log10 scale
signif_taxa$neg_log10_pvalue <- -log10(signif_taxa$padj)

# Filter signif_taxa based on cutoffs
#filtered_signif_taxa <- signif_taxa[(abs(signif_taxa$log2FoldChange) > 1) & (signif_taxa$neg_log10_pvalue > 5), ]

# Create a new dataframe with the relevant columns
result_order_diff <- signif_taxa[, c("Order", "log2FoldChange", "neg_log10_pvalue")]
result_order_diff$Order

volcano_rh_plants <- EnhancedVolcano(signif_taxa, 
                                     lab = signif_taxa$Order, 
                                     selectLab = signif_taxa$Order,
                                     x = 'log2FoldChange', 
                                     y = "padj",
                                     title = "Differential abundance in plant rhizosphere",
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
volcano_rh_plants

ggsave(filename = here("Results/Figures/", "deseq_volcano_rhizo_plant.png"),height = 7.5, width = 10.5, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_volcano_rhizo_plant.pdf"),height = 7.5, width = 10.5, dpi = 300)



#check

ps_scirpus <- subset_samples(ps_obj,Plant.type =="Scirpus")
ps_triglo  <- subset_samples(ps_obj,Plant.type =="Triglochin")

ps_scirpus_Acidithiobacillales <- subset_taxa(ps_scirpus,Order =="Acidithiobacillales")
ps_triglo_Acidithiobacillales <- subset_taxa(ps_triglo,Order =="Acidithiobacillales")

Ab_scirpus_acidi<- sum(rowSums(otu_table(ps_scirpus_Acidithiobacillales))) # 854 reads of Desulfobulbales in cold T°
Ab_triglo_acidi<- sum(rowSums(otu_table(ps_triglo_Acidithiobacillales))) # 854 reads of Desulfobulbales in cold T°


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


signif_taxa$Temp <- "More abundant in warm conditions"
for (i in 1:nrow(signif_taxa)){
  if(signif_taxa$log2FoldChange[i]<0){
    signif_taxa$Temp[i] <- "More abundant in cold conditions"
  }
}
ggplot(signif_taxa,aes(x=log2FoldChange,y=Order,color=Temp))+
  theme_bw()+
  facet_grid(Phylum~.,scale='free',space='free',switch='y')+
  theme(strip.text.y.left = element_text(angle = 0))+
  geom_vline(xintercept = 0)+
  scale_color_manual(values=c("#0000FE","#B02223"),name="Differential abundance")+
  geom_point(size=3)

ggsave(filename = here("Results/Figures/", "deseq_diffab_rhizo_temp.png"),height = 7.5, width = 10.5, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_diffab_rhizo_temp.pdf"),height = 7.5, width = 10.5, dpi = 300)

#### ¨¨¨¨¨¨Volcano plot #####
# Convert p-values to -log10 scale
signif_taxa$neg_log10_pvalue <- -log10(signif_taxa$padj)

# Filter signif_taxa based on cutoffs
#filtered_signif_taxa <- signif_taxa[(abs(signif_taxa$log2FoldChange) > 1) & (signif_taxa$neg_log10_pvalue > 5), ]

# Create a new dataframe with the relevant columns
result_order_diff <- signif_taxa[, c("Order", "log2FoldChange", "neg_log10_pvalue")]
result_order_diff$Order

volcano_rh_temp <- EnhancedVolcano(signif_taxa, 
                                   lab = signif_taxa$Order, 
                                   selectLab = signif_taxa$Order,
                                   x = 'log2FoldChange', 
                                   y = "padj",
                                   title = "Differential abundance in rhizosphere",
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

volcano_rh_temp

ggsave(filename = here("Results/Figures/", "deseq_volcano_rhizo_temp.png"),height = 7.5, width = 10.5, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_volcano_rhizo_temp.pdf"),height = 7.5, width = 10.5, dpi = 300)

### sanity check #
ps_warm<- subset_samples(ps_rhizo,Temperature=="20°C day 10°C night")
ps_cold<- subset_samples(ps_rhizo,Temperature=="10°C day 5°C night")


ps_xantho_cold <- subset_taxa(ps_cold,Order=="Xanthomonadales")
ps_xantho_warm <- subset_taxa(ps_warm,Order=="Xanthomonadales")

Ab_desulfo_cold<- sum(rowSums(otu_table(ps_xantho_cold))) # 854 reads of Desulfobulbales in cold T°
Ab_desulfo_warm <- sum(rowSums(otu_table(ps_xantho_warm))) # 3500 reads of Desulfobulbales in warm T°

####**Time D0-D62  ####

#Subset ps_rhizo to keep samples from D0 and D62 sampling dates
ps_obj_d0d62 <-  subset_samples(ps_rhizo, Time =="D0" | Time =="D62")
ps_obj_d0d62 <- prune_taxa(taxa_sums(ps_obj_d0d62)>0, ps_obj_d0d62)

# agglomerate t order level
ps_obj_d0d62 <- tax_glom(ps_obj_d0d62,taxrank="Order",NArm = F)

# remove taxa with not present in 90% of the samples
ps_obj <- ps_obj_d0d62
ps_obj <- prune_taxa(rowSums(otu_table(ps_obj) == 0) 
                     < ncol(otu_table(ps_obj)) * 0.9, ps_obj)
# Results in the loss of 47 orders

diagdds = phyloseq_to_deseq2(ps_obj, ~ Time)
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

ggsave(filename = here("Results/Figures/", "deseq_diffab_rhizo_d0-62.png"),height = 7.5, width = 10.5, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_diffab_rhizo_d0-62.pdf"),height = 7.5, width = 10.5, dpi = 300)

#### ¨¨¨¨¨¨Volcano plot #####
# Convert p-values to -log10 scale
signif_taxa$neg_log10_pvalue <- -log10(signif_taxa$padj)

# Filter signif_taxa based on cutoffs
#filtered_signif_taxa <- signif_taxa[(abs(signif_taxa$log2FoldChange) > 1) & (signif_taxa$neg_log10_pvalue > 5), ]

# Create a new dataframe with the relevant columns
result_order_diff <- signif_taxa[, c("Order", "log2FoldChange", "neg_log10_pvalue")]
result_order_diff$Order

volcano_rh_time <- EnhancedVolcano(signif_taxa, 
                                   lab = signif_taxa$Order, 
                                   selectLab = signif_taxa$Order,
                                   x = 'log2FoldChange', 
                                   y = "padj",
                                   title = "Differential abundance in rhizosphere",
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

volcano_rh_time

ggsave(filename = here("Results/Figures/", "deseq_volcano_rhizo_d0-62.png"),height = 7.5, width = 10.5, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_volcano_rhizo_d0-62.pdf"),height = 7.5, width = 10.5, dpi = 300)


####**Time D0-D35  ####
#Subset ps_rhizo to keep samples from D0 and D35 sampling dates
ps_obj_d0d35 <-  subset_samples(ps_rhizo, Time =="D0" | Time =="D35")
ps_obj_d0d35 <- prune_taxa(taxa_sums(ps_obj_d0d35)>0, ps_obj_d0d35)

# agglomerate t order level
ps_obj_d0d35 <- tax_glom(ps_obj_d0d35,taxrank="Order",NArm = F)

# remove taxa with not present in 90% of the samples
ps_obj <- ps_obj_d0d35
ps_obj <- prune_taxa(rowSums(otu_table(ps_obj) == 0) 
                     < ncol(otu_table(ps_obj)) * 0.9, ps_obj)
# Results in the loss of 56 orders

diagdds = phyloseq_to_deseq2(ps_obj_d0d35, ~ Time)
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
# https://support.bioconductor.org/p/35246/#35250


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
    signif_taxa$Time[i] <- "More abundant in D35"
  }
}
ggplot(signif_taxa,aes(x=log2FoldChange,y=Order,color=Time))+
  theme_bw()+
  facet_grid(Phylum~.,scale='free',space='free',switch='y')+
  theme(strip.text.y.left = element_text(angle = 0))+
  geom_vline(xintercept = 0)+
  scale_color_manual(values=c("black","orange"),name="Differential abundance")+
  geom_point(size=3)

ggsave(filename = here("Results/Figures/", "deseq_diffab_rhizo_d0-35.png"),height = 7.5, width = 10.5, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_diffab_rhizo_d0-35.pdf"),height = 7.5, width = 10.5, dpi = 300)

#### ¨¨¨¨¨¨Volcano plot #####
# Convert p-values to -log10 scale
signif_taxa$neg_log10_pvalue <- -log10(signif_taxa$padj)

# Filter signif_taxa based on cutoffs
#filtered_signif_taxa <- signif_taxa[(abs(signif_taxa$log2FoldChange) > 1) & (signif_taxa$neg_log10_pvalue > 5), ]

# Create a new dataframe with the relevant columns
result_order_diff <- signif_taxa[, c("Order", "log2FoldChange", "neg_log10_pvalue")]
result_order_diff$Order

volcano_rh_time <- EnhancedVolcano(signif_taxa, 
                                   lab = signif_taxa$Order, 
                                   selectLab = signif_taxa$Order,
                                   x = 'log2FoldChange', 
                                   y = "padj",
                                   title = "Differential abundance in rhizosphere",
                                   subtitle = "D35 vs D0 - Order",
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

volcano_rh_time

ggsave(filename = here("Results/Figures/", "deseq_volcano_rhizo_d0-35.png"),height = 7.5, width = 10.5, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_volcano_rhizo_d0-35.pdf"),height = 7.5, width = 10.5, dpi = 300)

#check

ps_d0<- subset_samples(ps_rhizo,Time =="D0")
ps_d35  <- subset_samples(ps_rhizo,Time =="D35")

ps_d0_Desulfo <- subset_taxa(ps_d0,Order =="Desulfobulbales")
ps_d35_Desulfo <- subset_taxa(ps_d35,Order =="Desulfobulbales")

Ab_d0_Desulfo<- sum(rowSums(otu_table(ps_d0_Desulfo))) # 743 reads of Desulfobulbales
#mostly present in one sample

Ab_d35_Desulfo<- sum(rowSums(otu_table(ps_d35_Desulfo))) # 1711 reads of Desulfobulbales 
#mostly present in 3 samples


sample_data(ps_obj)



####**Time D35-D62  ####

#Subset ps_rhizo to keep samples from D35 and D62 sampling dates
ps_obj_d35d62 <-  subset_samples(ps_rhizo, Time =="D35" | Time =="D62")
ps_obj_d35d62 <- prune_taxa(taxa_sums(ps_obj_d35d62)>0, ps_obj_d35d62)

# agglomerate t order level
ps_obj_d35d62 <- tax_glom(ps_obj_d35d62,taxrank="Order",NArm = F)

# remove taxa with not present in 90% of the samples
ps_obj <- ps_obj_d35d62
ps_obj <- prune_taxa(rowSums(otu_table(ps_obj) == 0) 
                     < ncol(otu_table(ps_obj)) * 0.9, ps_obj)
# Results in the loss of 53 order


diagdds = phyloseq_to_deseq2(ps_obj_d35d62, ~ Time)
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
# https://support.bioconductor.org/p/35246/#35250


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

summary(res2) # No significant difference in abudance



###### test ######

# ab_data <- data.frame(otu_table(ps))
# nrow(ab_data)*ncol(ab_data)
# zero_persample <- colSums(ab_data==0)
# ((sum(zero_persample))/(nrow(ab_data)*ncol(ab_data)))*100
# 
# ps_rhizo <-  subset_samples(ps, Sample.type =="Rhizosphere")
# ps_rhizo <- prune_taxa(taxa_sums(ps_rhizo)>0, ps_rhizo)
# ab_data <- data.frame(otu_table(ps_rhizo))
# nrow(ab_data)*ncol(ab_data)
# zero_persample <- colSums(ab_data==0)
# ((sum(zero_persample))/(nrow(ab_data)*ncol(ab_data)))*100

