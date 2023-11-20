### Packages#####
library(here)
source(here("source", "libraries.R"))
library(phyloseq)
library(viridis)
library(dplyr)
library(ggpp)
library(EnhancedVolcano)
library(metacoder)
#### Data ####

load(here("Rdata","ps_16S_water_obj.RData"))

###____####
##DESeq2 #####

#Order agglomeration
ps_glom_ord <- tax_glom(ps,taxrank="Order",NArm = F)
# 150  orders in 139 samples
ps_glom_gen <- tax_glom(ps,taxrank="Genus",NArm = F)
#451 generea in 139 samples 
#A <- data.frame(tax_table(ps_glom_gen))

# Run to select which object you want to run the analysis on
#ps_obj <- ps
ps_obj <- ps_glom_ord
ps_obj <- ps_glom_gen
# Filter out features (OTUs or taxa) that have more 
#than 90% zeros across all samples
# This is done by checking the number of zeros in each row
#of the OTU table and comparing it to 90% of the total
#number of columns (samples)
ps_obj <- prune_taxa(rowSums(otu_table(ps_obj) == 0) 
                     < ncol(otu_table(ps_obj)) * 0.9, ps_obj)

# Removes 49 orders if applied on order agglomerated data

######Sample Type ####

####|Scirpus VS Unplanted ####
ps_obj_NP_Scirpus <-  subset_samples(ps_obj, Sample_type=="No_plant" | Sample_type =="Scirpus")
ps_obj_NP_Scirpus <- prune_taxa(taxa_sums(ps_obj_NP_Scirpus)>0, ps_obj_NP_Scirpus)
ps_obj_NP_Scirpus <- prune_taxa(rowSums(otu_table(ps_obj_NP_Scirpus) == 0) 
                     < ncol(otu_table(ps_obj_NP_Scirpus)) * 0.9, ps_obj_NP_Scirpus)

diagdds = phyloseq_to_deseq2(ps_obj_NP_Scirpus, ~ Sample_type)

diagdds = estimateSizeFactors(diagdds, type = "poscounts")

diagdds = DESeq(diagdds,
                test = "Wald",
                fitType="parametric")
plotDispEsts(diagdds)

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

# Simple graph
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


#ggsave(filename = here("Results/Figures/", "deseq_diffab_water_ScirpusVSnoplant.png"),height = 7.5, width = 10.5, dpi = 300)
#ggsave(filename = here("Results/Figures/", "deseq_diffab_water_ScirpusVSnoplant.pdf"),height = 7.5, width = 10.5, dpi = 300)

# Volcano plots

# Convert p-values to -log10 scale
signif_taxa$neg_log10_pvalue <- -log10(signif_taxa$padj)

# Filter signif_taxa based on cutoffs
filtered_signif_taxa <- signif_taxa[(abs(signif_taxa$log2FoldChange) > 1) & (signif_taxa$neg_log10_pvalue > 2), ]
# -log10 p = 2 means there's a FDR of 1% 


# Create a new dataframe with the relevant columns
result_order_diff <- signif_taxa[, c("Order", "log2FoldChange", "neg_log10_pvalue")]
result_order_diff$Order

volcano_water_ScVSNp <- EnhancedVolcano(signif_taxa, 
                                          lab =signif_taxa$Order, 
                                          selectLab =signif_taxa$Order,
                                          x = 'log2FoldChange', 
                                          y = "padj",
                                          pCutoff = 1e-02,
                                          title = "Differential abundance in water",
                                          subtitle = "Scirpus vs Unplanted - Order",
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

ggsave(filename = here("Results/Figures/", "deseq_volcano_16S_water_ScirpusVSnoplant.png"),height = 7.5, width = 10.5, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_volcano_16S_water_ScirpusVSnoplant.pdf"),height = 7.5, width = 10.5, dpi = 300)


# Check outliers 
#ps_long <- psmelt(ps)
Obscuri <- subset_taxa(ps,Order=="Obscuribacterales")

p <- plot_bar(Obscuri,"Sample",fill="Genus")
p + theme_bw()+
    facet_grid(~Sample_type,drop=T,scales = "free")+
    geom_bar(aes(fill=Genus), stat="identity", position="stack")





######|Triglochin VS Unplanted ####
ps_obj_NP_Triglochin <-  subset_samples(ps_obj, Sample_type=="No_plant" | Sample_type =="Triglochin")
ps_obj_NP_Triglochin <- prune_taxa(taxa_sums(ps_obj_NP_Triglochin)>0, ps_obj_NP_Triglochin)
ps_obj_NP_Triglochin <- prune_taxa(rowSums(otu_table(ps_obj_NP_Triglochin) == 0) 
                                   < ncol(otu_table(ps_obj_NP_Triglochin)) * 0.9, ps_obj_NP_Triglochin)


diagdds = phyloseq_to_deseq2(ps_obj_NP_Triglochin, ~ Sample_type)

diagdds = estimateSizeFactors(diagdds, type = "poscounts")

diagdds = DESeq(diagdds,
                test = "Wald",
                fitType="parametric")
plotDispEsts(diagdds)

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

# Simple plot FDR=0.05
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


#ggsave(filename = here("Results/Figures/", "deseq_diffab_water_TriglochinVSnoplant.png"),height = 7.5, width = 10.5, dpi = 300)
#ggsave(filename = here("Results/Figures/", "deseq_diffab_water_TriglochinVSnoplant.pdf"),height = 7.5, width = 10.5, dpi = 300)


# Volcano plot
# Convert p-values to -log10 scale
signif_taxa$neg_log10_pvalue <- -log10(signif_taxa$padj)

# Filter signif_taxa based on cutoffs
filtered_signif_taxa <- signif_taxa[(abs(signif_taxa$log2FoldChange) > 1) & (signif_taxa$neg_log10_pvalue > 2), ]
# -log10 p = 2 means there's a FDR of 1% 

# Create a new dataframe with the relevant columns
result_order_diff <- signif_taxa[, c("Order", "log2FoldChange", "neg_log10_pvalue")]
result_order_diff$Order

volcano_water_TgVSNp <- EnhancedVolcano(signif_taxa, 
                                        lab = signif_taxa$Order, 
                                        selectLab = signif_taxa$Order,
                                        #xlim=c(-5,5),
                                        #ylim=c(0,20),
                                        x = 'log2FoldChange', 
                                        y = "padj",
                                        pCutoff = 1e-02,
                                        title = "Differential abundance in water",
                                        subtitle = "Triglochin vs Unplanted - Order",
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


ggsave(filename = here("Results/Figures/", "deseq_volcano_16S_water_TriglochinVSnoplant.png"),height = 7.5, width = 10.5, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_volcano_16S_water_TriglochinVSnoplant.pdf"),height = 7.5, width = 10.5, dpi = 300)

# Check outliers
Cyano <- subset_taxa(ps,Order=="Cyanobacteriales")

p <- plot_bar(Cyano,"Sample",fill="Genus")
p + theme_bw()+
    facet_grid(~Sample_type,drop=T,scales = "free")+
    geom_bar(aes(fill=Genus), stat="identity", position="stack")

Pseudo <- subset_taxa(ps,Order=="Pseudomonadales")

p <- plot_bar(Pseudo,"Sample",fill="Genus")
p + theme_bw()+
    facet_grid(~Sample_type,drop=T,scales = "free")+
    geom_bar(aes(fill=Genus), stat="identity", position="stack")


Pseuda <- subset_taxa(ps,Order=="Pseudanabaenales")
p <- plot_bar(Pseuda,"Sample",fill="Genus")
p + theme_bw()+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    facet_grid(~Sample_type,drop=T,scales = "free")+
    geom_bar(aes(fill=Genus), stat="identity", position="stack")


######|Scirpus VS Triglochin ####
ps_obj_Triglo_Scirpus <-  subset_samples(ps_obj, Sample_type=="Triglochin" | Sample_type =="Scirpus")
ps_obj_Triglo_Scirpus <- prune_taxa(taxa_sums(ps_obj_Triglo_Scirpus)>0, ps_obj_Triglo_Scirpus)
ps_obj_Triglo_Scirpus <- prune_taxa(rowSums(otu_table(ps_obj_Triglo_Scirpus) == 0) 
           < ncol(otu_table(ps_obj_Triglo_Scirpus)) * 0.9, ps_obj_Triglo_Scirpus)
      
      
diagdds = phyloseq_to_deseq2(ps_obj_Triglo_Scirpus, ~ Sample_type)

diagdds = estimateSizeFactors(diagdds, type = "poscounts")


diagdds = DESeq(diagdds,
                test = "Wald",
                fitType="parametric")
plotDispEsts(diagdds)


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


#ggsave(filename = here("Results/Figures/", "deseq_diffab_water_TrigloVSScirpus.png"),height = 7.5, width = 10.5, dpi = 300)
#ggsave(filename = here("Results/Figures/", "deseq_diffab_water_TrigloVSScirpus.pdf"),height = 7.5, width = 10.5, dpi = 300)


# Volcano plot

# Convert p-values to -log10 scale
signif_taxa$neg_log10_pvalue <- -log10(signif_taxa$padj)

# Filter signif_taxa based on cutoffs
filtered_signif_taxa <- signif_taxa[(abs(signif_taxa$log2FoldChange) > 1) & (signif_taxa$neg_log10_pvalue > 2), ]
# -log10 p = 2 means there's a FDR of 1% 

# Create a new dataframe with the relevant columns
result_order_diff <- signif_taxa[, c("Order", "log2FoldChange", "neg_log10_pvalue")]
result_order_diff$Order

volcano_water_TgVSSc <- EnhancedVolcano(signif_taxa, 
                                        lab = signif_taxa$Order, 
                                        selectLab = signif_taxa$Order,
                                        x = 'log2FoldChange', 
                                        y = "padj",
                                        pCutoff = 1e-02,
                                        title = "Differential abundance in water",
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
volcano_water_TgVSSc

# 
ggsave(filename = here("Results/Figures/", "deseq_volcano_16S_water_TrigloVSScirpus.png"),height = 7.5, width = 10.5, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_volcano_16S_water_TrigloVSScirpus.pdf"),height = 7.5, width = 10.5, dpi = 300)



# Check outliers
Pseudo <- subset_taxa(ps,Order=="Pseudomonadales")

p <- plot_bar(Pseudo,"Sample",fill="Genus")
p + theme_bw()+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    facet_grid(~Sample_type,drop=T,scales = "free")+
    geom_bar(aes(fill=Genus), stat="identity", position="stack")

Gammapro <- subset_taxa(ps,Order=="Gammaproteobacteria")

p <- plot_bar(Gammapro,"Sample",fill="Genus")
p + theme_bw()+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    facet_grid(~Sample_type,drop=T,scales = "free")+
    geom_bar(aes(fill=Genus), stat="identity", position="stack")

Elsterales <- subset_taxa(ps,Order=="Elsterales")

p <- plot_bar(Elsterales,"Sample",fill="Genus")
p + theme_bw()+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    facet_grid(~Sample_type,drop=T,scales = "free")+
    geom_bar(aes(fill=Genus), stat="identity", position="stack")




####Temperature  ####
diagdds = phyloseq_to_deseq2(ps_obj, ~ Temperature)

diagdds = estimateSizeFactors(diagdds, type = "poscounts")


diagdds = DESeq(diagdds,
                test = "Wald",
                fitType="parametric")
plotDispEsts(diagdds)

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

#### Simple plot ##
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

#ggsave(filename = here("Results/Figures/", "deseq_diffab_water_temp.png"),height = 7.5, width = 10.5, dpi = 300)
#ggsave(filename = here("Results/Figures/", "deseq_diffab_water_temp.pdf"),height = 7.5, width = 10.5, dpi = 300)

#### Volcano plot
# Convert p-values to -log10 scale
signif_taxa$neg_log10_pvalue <- -log10(signif_taxa$padj)

# Filter signif_taxa based on cutoffs
filtered_signif_taxa <- signif_taxa[(abs(signif_taxa$log2FoldChange) > 1) & (signif_taxa$neg_log10_pvalue > 2), ]

# Create a new dataframe with the relevant columns
result_order_diff <- signif_taxa[, c("Order", "log2FoldChange", "neg_log10_pvalue")]
result_order_diff$Order

volcano_water_temp <- EnhancedVolcano(signif_taxa, 
                                     lab = signif_taxa$Order, 
                                     selectLab = signif_taxa$Order,
                                     x = 'log2FoldChange', 
                                     y = "padj",
                                     pCutoff = 1e-02,
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

ggsave(filename = here("Results/Figures/", "deseq_volcano_16S_water_temp.png"),height = 7.5, width = 10.5,dpi=300)
ggsave(filename = here("Results/Figures/", "deseq_volcano_16S_water_temp.pdf"),height = 7.5, width = 10.5,dpi=300)


Synechococcales <- subset_taxa(ps,Order=="Synechococcales")
p <- plot_bar(Synechococcales,"Sample",fill="Genus")
p + theme_bw()+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    facet_grid(~Temperature,drop=T,scales = "free")+
    geom_bar(aes(fill=Genus), stat="identity", position="stack")


Rhodospirillales <- subset_taxa(ps,Order=="Rhodospirillales")
p <- plot_bar(Araeolaimida,"Sample",fill="Genus")
p + theme_bw()+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    facet_grid(~Temperature,drop=T,scales = "free")+
    geom_bar(aes(fill=Genus), stat="identity", position="stack")




####Time  ####
#####|D0-D62 ####
#Subset ps_rhizo to keep samples from D0 and D62 sampling dates
ps_obj_d0d62 <-  subset_samples(ps_obj, Time =="D0" | Time =="D62")
ps_obj_d0d62 <- prune_taxa(taxa_sums(ps_obj_d0d62)>0, ps_obj_d0d62)
ps_obj_d0d62 <- prune_taxa(rowSums(otu_table(ps_obj_d0d62) == 0) 
                                    < ncol(otu_table(ps_obj_d0d62)) * 0.9, ps_obj_d0d62)


diagdds = phyloseq_to_deseq2(ps_obj_d0d62, ~ Time)

diagdds = estimateSizeFactors(diagdds, type = "poscounts")

diagdds = DESeq(diagdds,
                test = "Wald",
                fitType="parametric")
plotDispEsts(diagdds)


results(diagdds)
res = results(diagdds, alpha=0.05) # results of the test FDR accepted = 5% 
plotMA(res)

res2 = res[order(res$padj, na.last=NA), ] # ordered and removes the ASVs which have padj = NAs 
plotMA(res2)

summary(res2)


# Simple plot
#sigtab = res2[(res2$padj < 0.05), ] # Select significant padj
#sigtab <- as.data.frame(sigtab)
#sigtab <- cbind(rownames(sigtab),sigtab)
#colnames(sigtab)[1] <- "OTU"
#
#tax <- data.frame(tax_table(ps_obj))
#tax_diffab <- subset(tax,rownames(tax)%in%rownames(sigtab))
#
#signif_taxa <- merge(tax_diffab,sigtab,by='row.names')

# 
# signif_taxa$Time <- "More abundant in D62"
# for (i in 1:nrow(signif_taxa)){
#   if(signif_taxa$log2FoldChange[i]<0){
#     signif_taxa$Time[i] <- "More abundant in D0"
#   }
# }
# ggplot(signif_taxa,aes(x=log2FoldChange,y=Order,color=Time))+
#   theme_bw()+
#   facet_grid(Phylum~.,scale='free',space='free',switch='y')+
#   theme(strip.text.y.left = element_text(angle = 0))+
#   geom_vline(xintercept = 0)+
#   scale_color_manual(values=c("black","yellow"),name="Differential abundance")+
#   geom_point(size=3)
# 
# ggsave(filename = here("Results/Figures/", "deseq_diffab_roots_time.png"),height = 7.5, width = 10.5, dpi = 300)
# ggsave(filename = here("Results/Figures/", "deseq_diffab_roots_time.pdf"),height = 7.5, width = 10.5, dpi = 300)
# 
# 
# Volcano plot 
# Convert p-values to -log10 scale
res2 <- as.data.frame(res2)
tax <- data.frame(tax_table(ps_obj))
res2_taxa <- merge(res2,tax,by='row.names')


res2_taxa$neg_log10_pvalue <- -log10(res2_taxa$padj)

# Filter signif_taxa based on cutoffs
filtered_signif_taxa <- res2_taxa[(abs(res2_taxa$log2FoldChange) > 1) & (res2_taxa$neg_log10_pvalue > 5), ]

# Create a new dataframe with the relevant columns
result_diff <- filtered_signif_taxa[, c("Genus", "log2FoldChange", "neg_log10_pvalue")]
result_diff$Genus

volcano_water_d0d62<- EnhancedVolcano(res2_taxa, 
                                      lab = res2_taxa$Genus, 
                                      selectLab = res2_taxa$Genus,
                                      x = 'log2FoldChange', 
                                      y = "padj",
                                      pCutoff = 1e-05,
                                      title = "Differential abundance of 16S reads in water",
                                      subtitle = "D62 vs D0 - Genus level",
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

volcano_water_d0d62


ggsave(filename = here("Results/Figures/", "deseq_volcano_16S_water_Genusd0d62.png"), plot = volcano_water_d0d62,height = 7.5, width = 10.5, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_volcano_16S_water_Genusd0d62.pdf"), plot = volcano_water_d0d62,height = 7.5, width = 10.5, dpi = 300)

#####|D0-D77 ####
#Subset ps_rhizo to keep samples from D0 and D62 sampling dates
ps_obj_d0d77 <-  subset_samples(ps_obj, Time =="D0" | Time =="D77")
ps_obj_d0d77 <- prune_taxa(taxa_sums(ps_obj_d0d77)>0, ps_obj_d0d77)
ps_obj_d0d77 <- prune_taxa(rowSums(otu_table(ps_obj_d0d77) == 0) 
                           < ncol(otu_table(ps_obj_d0d77)) * 0.9, ps_obj_d0d77)
#244 genera in 45 samples
#
diagdds = phyloseq_to_deseq2(ps_obj_d0d77, ~ Time)
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
                fitType="parametric")
plotDispEsts(diagdds)


results(diagdds)
res = results(diagdds, alpha=0.05) # results of the test FDR accepted = 5% 
plotMA(res)

res2 = res[order(res$padj, na.last=NA), ] # ordered and removes the ASVs which have padj = NAs 
plotMA(res2)

summary(res2)


# Simple plot
 
# sigtab = res2[(res2$padj < 0.05), ] # Select significant padj
# sigtab <- as.data.frame(sigtab)
# sigtab <- cbind(rownames(sigtab),sigtab)
# colnames(sigtab)[1] <- "OTU"
# 
# tax <- data.frame(tax_table(ps_obj))
# tax_diffab <- subset(tax,rownames(tax)%in%rownames(sigtab))
# 
# signif_taxa <- merge(tax_diffab,sigtab,by='row.names')
# 
# signif_taxa$Time <- "More abundant in D77"
# for (i in 1:nrow(signif_taxa)){
#     if(signif_taxa$log2FoldChange[i]<0){
#         signif_taxa$Time[i] <- "More abundant in D0"
#     }
# }
# ggplot(signif_taxa,aes(x=log2FoldChange,y=Order,color=Time))+
#     theme_bw()+
#     facet_grid(Phylum~.,scale='free',space='free',switch='y')+
#     theme(strip.text.y.left = element_text(angle = 0))+
#     geom_vline(xintercept = 0)+
#     scale_color_manual(values=c("black","yellow"),name="Differential abundance")+
#     geom_point(size=3)
# 
#ggsave(filename = here("Results/Figures/", "deseq_diffab_roots_time.png"),height = 7.5, width = 10.5, dpi = 300)
#ggsave(filename = here("Results/Figures/", "deseq_diffab_roots_time.pdf"),height = 7.5, width = 10.5, dpi = 300)


# Volcano plot 
# Convert p-values to -log10 scale

res2 <- as.data.frame(res2)
tax <- data.frame(tax_table(ps_obj))
res2_taxa <- merge(res2,tax,by='row.names')


res2_taxa$neg_log10_pvalue <- -log10(res2_taxa$padj)

# Filter signif_taxa based on cutoffs
filtered_signif_taxa <- res2_taxa[(abs(res2_taxa$log2FoldChange) > 1) & (res2_taxa$neg_log10_pvalue > 5), ]

# Create a new dataframe with the relevant columns
result_diff <- filtered_signif_taxa[, c("Genus", "log2FoldChange", "neg_log10_pvalue")]
result_diff$Genus

volcano_water_d0d77<- EnhancedVolcano(res2_taxa, 
                                      lab = res2_taxa$Genus, 
                                      selectLab = res2_taxa$Genus,
                                      x = 'log2FoldChange', 
                                      y = "padj",
                                      pCutoff = 1e-05,
                                      title = "Differential abundance of 16S reads in water",
                                      subtitle = "D77 vs D0 - Genus level",
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

volcano_water_d0d77

ggsave(filename = here("Results/Figures/", "deseq_volcano_16S_water_Genusd0d77.png"), plot = volcano_water_d0d77,height = 7.5, width = 10.5, dpi = 300)
ggsave(filename = here("Results/Figures/", "deseq_volcano_16S_water_Genusd0d77.pdf"), plot = volcano_water_d0d77,height = 7.5, width = 10.5, dpi = 300)

####Metacoder####
load(here("Rdata","ps_16S_water_obj.RData"))

### ....Metacoder Object####
# Taxonomic agglomeration
#ps_glom_order <- tax_glom(ps,taxrank="Order")
ps_glom_class <- tax_glom(ps,taxrank="Class")

ps_obj <- ps_glom_class

# Relative abundance transformation
ps.relab <- transform_sample_counts(ps_obj,function(x) x/sum(x))

meta_obj <- parse_phyloseq(ps.relab) 
# transforms the phyloseq object in a metacoder object

### ....Heat tree #####

### Sample Type ###
meta_obj$data$tax_abund <- calc_taxon_abund(meta_obj, "otu_table")


meta_obj$data$sample_data$Time <- factor(as.character(meta_obj$data$sample_data$Time), 
                                         levels = c("D0","D10","D16","D35","D62","D77"),
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
                  key_size = 0.8,
                  label_small_trees = F,
                  layout = "davidson-harel", # The primary layout algorithm
                  initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                  output_file = here("Results","Figures","Watwr16S_Time_DiffTree.pdf")) # Saves the plot as a pdf file

