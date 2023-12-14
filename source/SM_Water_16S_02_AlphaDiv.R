### Packages#####
library(here)
source(here("source", "libraries.R"))
library(phyloseq)
library(viridis)
library(metacoder)
library(dplyr)
library(ggpp)

#### Data ####

load(here("Rdata","ps_16S_water_obj.RData"))

#### Separation by temperature #####

ps_water_16S_Warm <- subset_samples(ps, Temperature =="20C day - 10C night")
ps_water_16S_Warm <- prune_taxa(taxa_sums(ps_water_16S_Warm)>0, ps_water_16S_Warm)
ps_water_16S_Warm # 68 samples - 2208 taxa
# check to see if you kept only water_16S samples in WARM condition
str(sample_data(ps_water_16S_Warm)) 

ps_water_16S_Cold <- subset_samples(ps, Temperature =="10C day - 5C night")
ps_water_16S_Cold <- prune_taxa(taxa_sums(ps_water_16S_Cold)>0, ps_water_16S_Cold)
ps_water_16S_Cold # 71 samples - 2019 taxa
# check to see if you kept only water_16S samples in COLD condition
str(sample_data(ps_water_16S_Cold)) 

####~~~~~~~~~~~~~~~~####
### ALPHA DIV ####
#### Raw abundance ####

#####Plant type ####
######| Warm ####
alpha.div.warm<- plot_richness(ps_water_16S_Warm, x="Sample_type", measures=c("Observed","Simpson", "Shannon"))

# Simplified data for plotting
plot <- select(alpha.div.warm$data,c("samples","Time","Sample_type","variable","value"))

Warm_PlantType <- ggplot(plot,aes(x=variable,y=value,shape=Sample_type))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(aes(x=variable,y=value,shape=Sample_type),position=position_dodge(width=0.75),size=2.5,fill='black')+
    theme_bw()+
    theme(title = element_text(size=20),
          axis.text.x =element_blank(),
          axis.text.y =element_text(size=12),
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size=15),
          legend.text.align = 0,
          axis.ticks.x = element_blank(),
          legend.text=element_text(size=20),
          legend.position = "bottom",
          strip.text = element_text(size=15),
          strip.background = element_rect(fill='white'))+
    stat_pwc(method="wilcox_test",hide.ns=T,
             p.adjust.method="fdr",
             label="p.adj.signif",
             tip.length = 0, 
             step.increase = 0.05,
             vjust=0)+
    # Computes pairwise wilcoxon rank-sum test in each facet
    scale_shape_manual(values=c(19,22,23),name="Plant type",
                       labels=c("No plant",expression(paste(italic("Scirpus"))),expression(italic("Triglochin"))))+
    facet_wrap(~variable,scales="free")+
    ylab("Alpha diversity score")

# Save generated figure
# ggsave(here("Results_W&C","Figures","water_16S_Warm_SType_AlphaDiv.pdf"), plot = Warm_PlantType ,device='pdf',height = 7.5, width = 10.5)
# ggsave(here("Results_W&C","Figures","water_16S_Warm_SType_AlphaDiv.png"), plot = Warm_PlantType, device='png',height = 7.5, width = 10.5)

######| Cold ####
alpha.div.Cold<- plot_richness(ps_water_16S_Cold, x="Sample_type", measures=c("Observed","Simpson", "Shannon"))

# Simplified data for plotting
plot <- select(alpha.div.Cold$data,c("samples","Time","Sample_type","variable","value"))

Cold_PlantType <- ggplot(plot,aes(x=variable,y=value,shape=Sample_type))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(aes(x=variable,y=value,shape=Sample_type),position=position_dodge(width=0.75),size=2.5,fill='black')+
    theme_bw()+
    theme(title = element_text(size=20),
          axis.text.x =element_blank(),
          axis.text.y =element_text(size=12),
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size=15),
          legend.text.align = 0,
          axis.ticks.x = element_blank(),
          legend.text=element_text(size=20),
          legend.position = "bottom",
          strip.text = element_text(size=15),
          strip.background = element_rect(fill='white'))+
    stat_pwc(method="wilcox_test",hide.ns=T,
             p.adjust.method="fdr",
             label="p.adj.signif",
             tip.length = 0, 
             step.increase = 0.05,
             vjust=0)+
    # Computes pairwise wilcoxon rank-sum test in each facet
    scale_shape_manual(values=c(19,22,23),name="Plant type",
                       labels=c("No plant",expression(paste(italic("Scirpus"))),expression(italic("Triglochin"))))+
    facet_wrap(~variable,scales="free")+
    ylab("Alpha diversity score")

# Save generated figure
# ggsave(here("Results_W&C","Figures","water_16S_Cold_SType_AlphaDiv.pdf"), plot = Cold_PlantType ,device='pdf',height = 7.5, width = 10.5)
# ggsave(here("Results_W&C","Figures","water_16S_Cold_SType_AlphaDiv.png"), plot = Cold_PlantType, device='png',height = 7.5, width = 10.5)

Plant_type_alpha <- ggarrange(Warm_PlantType,Cold_PlantType,
                              nrow = 2,common.legend = T,
                              legend = "bottom",
                              labels = c("Warm","Cold"),
                              hjust = -0.1)

ggsave(here("Results_W&C","Figures","Water_16S_SType_AlphaDiv.pdf"), plot = Plant_type_alpha ,device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results_W&C","Figures","Water_16S_SType_AlphaDiv.png"), plot = Plant_type_alpha, device='png',height = 7.5, width = 10.5)


#####Time ####
######| Warm ####
alpha.div.warm<- plot_richness(ps_water_16S_Warm, x="Time", measures=c("Observed","Simpson", "Shannon"))

# Simplified data for plotting
plot <- select(alpha.div.warm$data,c("samples","Time","Sample_type","variable","value"))

Warm_Time <- ggplot(plot,aes(x=variable,y=value,color=Time))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(aes(x=variable,y=value),position=position_dodge(width=0.75),size=2.5,fill='black')+
    theme_bw()+
    theme(title = element_text(size=20),
          axis.text.x =element_blank(),
          axis.text.y =element_text(size=12),
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size=15),
          legend.text.align = 0,
          axis.ticks.x = element_blank(),
          legend.text=element_text(size=20),
          legend.position = "bottom",
          strip.text = element_text(size=15),
          strip.background = element_rect(fill='white'))+
    stat_pwc(method="wilcox_test",
             hide.ns=T,
             p.adjust.method="fdr",
             label="p.adj.signif",
             tip.length = 0, 
             step.increase = 0.1,
             vjust=0.5,
             size = 0.5,
             label.size=4.2)+
    # Computes pairwise wilcoxon rank-sum test in each facet
    scale_colour_viridis(option="magma",discrete=T,name="Time")+
    facet_wrap(~variable,scales="free")+
    ylab("Alpha diversity score")

# Save generated figure
# ggsave(here("Results_W&C","Figures","water_16S_Warm_Time_AlphaDiv.pdf"), plot = Warm_PlantType ,device='pdf',height = 7.5, width = 10.5)
# ggsave(here("Results_W&C","Figures","water_16S_Warm_Time_AlphaDiv.png"), plot = Warm_PlantType, device='png',height = 7.5, width = 10.5)

######| Cold ####
alpha.div.Cold<- plot_richness(ps_water_16S_Cold, x="Time", measures=c("Observed","Simpson", "Shannon"))

# Simplified data for plotting
plot <- select(alpha.div.Cold$data,c("samples","Time","Sample_type","variable","value"))

Cold_Time <- ggplot(plot,aes(x=variable,y=value,color=Time))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(aes(x=variable,y=value),position=position_dodge(width=0.75),size=2.5,fill='black')+
    theme_bw()+
    theme(title = element_text(size=20),
          axis.text.x =element_blank(),
          axis.text.y =element_text(size=12),
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size=15),
          legend.text.align = 0,
          axis.ticks.x = element_blank(),
          legend.text=element_text(size=20),
          legend.position = "bottom",
          strip.text = element_text(size=15),
          strip.background = element_rect(fill='white'))+
    stat_pwc(method="wilcox_test",
             hide.ns=T,
             p.adjust.method="fdr",
             label="p.adj.signif",
             tip.length = 0, 
             step.increase = 0.1,
             vjust=0.5,
             size = 0.5,
             label.size=4.2)+
    # Computes pairwise wilcoxon rank-sum test in each facet
    scale_colour_viridis(option="magma",discrete=T,name="Time")+
    facet_wrap(~variable,scales="free")+
    ylab("Alpha diversity score")

# Save generated figure
# ggsave(here("Results_W&C","Figures","water_16S_Cold_Time_AlphaDiv.pdf"), plot = Cold_PlantType ,device='pdf',height = 7.5, width = 10.5)
# ggsave(here("Results_W&C","Figures","water_16S_Cold_Time_AlphaDiv.png"), plot = Cold_PlantType, device='png',height = 7.5, width = 10.5)

Time_alpha <- ggarrange(Warm_Time,Cold_Time,
                        nrow = 2,common.legend = T,
                        legend = "bottom",
                        labels = c("Warm","Cold"),
                        hjust = -0.1)

ggsave(here("Results_W&C","Figures","Water_16S_Time_AlphaDiv.pdf"), 
       plot = Time_alpha ,
       device='pdf',height = 7.5, width = 10.5)

ggsave(here("Results_W&C","Figures","water_16S_Time_AlphaDiv.png"), 
       plot = Time_alpha, 
       device='png',height = 7.5, width = 10.5)


##### Time & Plant type####
######| Warm ####
alpha.div.warm<- plot_richness(ps_water_16S_Warm, x="Time", measures=c("Observed","Simpson", "Shannon"))

# Simplified data for plotting
plot <- select(alpha.div.warm$data,c("samples","Time","Sample_type","variable","value"))

# Rename plant labels
Plant_labs <- c("No Plant","Scirpus", "Triglochin")
names(Plant_labs) <- c("No_plant","Scirpus", "Triglochin")


Warm_PT_Time <- ggplot(plot,aes(x=variable,y=value,color=Time))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(aes(x=variable,y=value),position=position_dodge(width=0.75),size=2.5,fill='black')+
    theme_bw()+
    theme(title = element_text(size=20),
          axis.text.x =element_blank(),
          axis.text.y =element_text(size=12),
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size=15),
          legend.text.align = 0,
          axis.ticks.x = element_blank(),
          legend.text=element_text(size=20),
          legend.position = "bottom",
          strip.text = element_text(size=15),
          strip.text.x = element_text(face="italic"),
          strip.background = element_rect(fill='white'))+
    stat_pwc(method="wilcox_test",
             hide.ns=T,
             p.adjust.method="fdr",
             label="p.adj.signif",
             tip.length = 0, 
             step.increase = 0.1,
             vjust=0.5,
             size = 0.5,
             label.size=4.2)+
    # Computes pairwise wilcoxon rank-sum test in each facet
    scale_colour_viridis(option="magma",discrete=T,name="Time")+
    facet_nested(variable~Sample_type,scales="free",
                 drop=T,independent = "all",shrink = F,
                 labeller = labeller(Sample_type=Plant_labs))+
    ylab("Alpha diversity score")

# Save generated figure
# ggsave(here("Results_W&C","Figures","water_16S_Warm_PT_Time_AlphaDiv.pdf"), plot = Warm_PlantType ,device='pdf',height = 7.5, width = 10.5)
# ggsave(here("Results_W&C","Figures","water_16S_Warm_PT_Time_AlphaDiv.png"), plot = Warm_PlantType, device='png',height = 7.5, width = 10.5)

######| Cold ####
alpha.div.Cold<- plot_richness(ps_water_16S_Cold, x="Time", measures=c("Observed","Simpson", "Shannon"))

# Simplified data for plotting
plot <- select(alpha.div.Cold$data,c("samples","Time","Sample_type","variable","value"))

Cold_PT_Time <- ggplot(plot,aes(x=variable,y=value,color=Time))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(aes(x=variable,y=value),position=position_dodge(width=0.75),size=2.5,fill='black')+
    theme_bw()+
    theme(title = element_text(size=20),
          axis.text.x =element_blank(),
          axis.text.y =element_text(size=12),
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size=15),
          legend.text.align = 0,
          axis.ticks.x = element_blank(),
          legend.text=element_text(size=20),
          legend.position = "bottom",
          strip.text = element_text(size=15),
          strip.text.x = element_text(face="italic"),
          strip.background = element_rect(fill='white'))+
    stat_pwc(method="wilcox_test",
             hide.ns=T,
             p.adjust.method="fdr",
             label="p.adj.signif",
             tip.length = 0, 
             step.increase = 0.1,
             vjust=0.5,
             size = 0.5,
             label.size=4.2)+
    # Computes pairwise wilcoxon rank-sum test in each facet
    scale_colour_viridis(option="magma",discrete=T,name="Time")+
    facet_nested(variable~Sample_type,scales="free",
                 drop=T,independent = "all",shrink = F,
                 labeller = labeller(Sample_type=Plant_labs))+
    ylab("Alpha diversity score")

# Save generated figure
# ggsave(here("Results_W&C","Figures","water_16S_Cold_PT_Time_AlphaDiv.pdf"), plot = Cold_PlantType ,device='pdf',height = 7.5, width = 10.5)
# ggsave(here("Results_W&C","Figures","water_16S_Cold_PT_Time_AlphaDiv.png"), plot = Cold_PlantType, device='png',height = 7.5, width = 10.5)

Time_PT_alpha <- ggarrange(Warm_PT_Time,Cold_PT_Time,
                           nrow = 2,common.legend = T,
                           legend = "bottom",
                           labels = c("Warm","Cold"),
                           hjust = -0.1)

ggsave(here("Results_W&C","Figures","Water_16S_PT_Time_AlphaDiv.pdf"), 
       plot = Time_PT_alpha ,
       device='pdf',height = 7.5, width = 10.5)

ggsave(here("Results_W&C","Figures","Water_16S_PT_Time_AlphaDiv.png"), 
       plot = Time_PT_alpha, 
       device='png',height = 7.5, width = 10.5)


#### Supp Figure ######

######| Warm ####
alpha.div.warm<- plot_richness(ps_water_16S_Warm, x="Time", measures= "Shannon")

# Simplified data for plotting
plot <- select(alpha.div.warm$data,c("samples","Time","Sample_type","variable","value"))

Warm <- ggplot(plot,aes(x=Time,y=value,fill=Time))+
    geom_point(aes(x=Time,y=value,shape=Sample_type,fill=Time),size=4,position=position_dodge(width=0.5))+
    geom_boxplot(outlier.shape = NA,alpha=0.6)+
    theme_bw()+
    theme(title = element_text(size=20),
          axis.text.x =element_blank(),
          axis.text.y =element_text(size=12),
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size=15),
          legend.text.align = 0,
          axis.ticks.x = element_blank(),
          legend.text=element_text(size=18),
          legend.title = element_text(size=20),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.margin=margin())+
    stat_pwc(method="wilcox_test",
             hide.ns=T,
             p.adjust.method="fdr",
             label="p.adj.signif",
             tip.length = 0, 
             step.increase = 0.05,
             vjust=0.5,
             size = 0.5,
             label.size=4.2)+
    # Computes pairwise wilcoxon rank-sum test in each facet
    scale_colour_viridis(option="magma",discrete=T,name="Time")+
    scale_fill_viridis(option="magma",discrete=T,name="Time")+
    scale_shape_manual(name="Plant type :",
                       labels=c("No plant",expression(italic(Scirpus)),expression(italic(Triglochin))),
                       values = c(21,22,23))+
    ylab("Alpha diversity score")
Warm



######| Cold ####
alpha.div.cold<- plot_richness(ps_water_16S_Cold, x="Time", measures= "Shannon")

# Simplified data for plotting
plot <- select(alpha.div.cold$data,c("samples","Time","Sample_type","variable","value"))

Cold <- ggplot(plot,aes(x=Time,y=value,fill=Time))+
    geom_point(aes(x=Time,y=value,shape=Sample_type,fill=Time),size=4,position=position_dodge(width=0.5))+
    geom_boxplot(outlier.shape = NA,alpha=0.6)+
    theme_bw()+
    theme(title = element_text(size=20),
          axis.text.x =element_blank(),
          axis.text.y =element_text(size=12),
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size=15),
          legend.text.align = 0,
          axis.ticks.x = element_blank(),
          legend.text=element_text(size=18),
          legend.title = element_text(size=20),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.margin=margin())+
    stat_pwc(method="wilcox_test",
             hide.ns=T,
             p.adjust.method="fdr",
             label="p.adj.signif",
             tip.length = 0, 
             step.increase = 0.05,
             vjust=0.5,
             size = 0.5,
             label.size=4.2)+
    # Computes pairwise wilcoxon rank-sum test in each facet
    scale_colour_viridis(option="magma",discrete=T,name="Time")+
    scale_fill_viridis(option="magma",discrete=T,name="Time")+
    scale_shape_manual(name="Plant type :",
                       labels=c("No plant",expression(italic(Scirpus)),expression(italic(Triglochin))),
                       values = c(21,22,23))+
    ylab("Alpha diversity score")
Cold

Alpha <- ggarrange(Warm,Cold,
                   ncol = 2,common.legend = T,
                   legend = "bottom", labels=c("A","B"))


ggsave(here("Results_W&C","Figures","SupFig_Water16SAlphaDiv.pdf"), 
       plot = Alpha ,
       device='pdf',height = 7.5, width = 10.5)

ggsave(here("Results_W&C","Figures","SupFig_Water16SAlphaDiv.png"), 
       plot = Alpha, 
       device='png',height = 7.5, width = 10.5)





####~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~####

#### Not separated by temperature #####
####~~~~~~~~~~~~~~~~####

#### Data ####

load(here("Rdata","ps_16S_water_obj.RData"))

###______#####
### ALPHA DIV ####
###Raw abundance ####

#####Sample type ####
alpha.div<- plot_richness(ps, x="Sample_type", color="Sample_type", measures=c("Observed","Simpson", "Shannon"))

plot <- select(alpha.div$data,c("samples","Time","Sample_type","Temperature","variable","value"))


ggplot(plot,aes(x=variable,y=value,shape=Sample_type))+
  geom_boxplot(outlier.shape = NA)+
  #geom_boxplot(outlier.shape = NA,position=position_dodge(width=1.5))+
  geom_point(aes(x=variable,y=value,shape=Sample_type),position=position_dodge(width=0.75),size=2.5,fill='black')+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank())+
  stat_pwc(method="wilcox_test",hide.ns=T,
           p.adjust.method="fdr",label="p.adj.signif",
           tip.length = 0, 
           step.increase = 0.05,
           vjust=0)+
  scale_shape_manual(values=c(19,22,23),name="Sample type")+
  facet_wrap(~variable,scales="free")+
  ylab("Alpha diversity score")

ggsave(here("Results","Figures","16S_Water_SType_AlphaDiv.pdf"),device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results","Figures","16S_Water_SType_AlphaDiv.png"),device='png',height = 7.5, width = 10.5)


#####Temperature ####
alpha.div<- plot_richness(ps, x="Temperature", color="Temperature", measures=c("Observed","Simpson", "Shannon"))

plot <- select(alpha.div$data,c("samples","Time","Sample_type","Temperature","variable","value"))



ggplot(plot,aes(x=variable,y=value,color=Temperature))+
  geom_boxplot(outlier.shape = NA)+
  #geom_boxplot(outlier.shape = NA,position=position_dodge(width=1.5))+
  geom_point(aes(x=variable,y=value),position=position_dodge(width=0.75),size=2.5,fill='black')+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank())+
  stat_pwc(method="wilcox_test",
           hide.ns=T,
           p.adjust.method="fdr",
           label="p.adj.signif",
           tip.length = 0, 
           step.increase = 0.5,
           vjust=0)+
  scale_color_manual(values=c("#0000FE","#B02223"),name="Temperature")+
  facet_wrap(~variable,scales="free")+
  ylab("Alpha diversity score")

ggsave(here("Results","Figures","16S_Water_Temp_AlphaDiv.pdf"),device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results","Figures","16S_Water_Temp_AlphaDiv.png"),device='png',height = 7.5, width = 10.5)

#####Time #####
library(ggh4x)
alpha.div <- plot_richness(ps, x="Time", measures=c("Observed","Simpson", "Shannon"))

plot <- select(alpha.div$data,c("samples","Time","Sample_type","Temperature","variable","value"))


ggplot(plot,aes(x=variable,y=value,shape=Sample_type,color=Time))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(x=variable,y=value,shape=Sample_type),position=position_dodge(width=0.75),size=2,fill='grey60')+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank())+
  stat_pwc(method="wilcox_test",
           p.adjust.method="fdr",
           hide.ns=T,
           label="p.adj.signif",
           tip.length = 0, 
           vjust=0.7,
           size = 0.5,
           label.size=4.2)+
  scale_colour_viridis(option="magma",discrete=T,name="Time")+
  scale_shape_manual(values=c(21,22,23),name="Sample type")+
  facet_nested(variable~Sample_type,scales="free",drop=T,independent = "all")+
  ylab("Alpha diversity score")


ggsave(here("Results","Figures","16S_Water_TimeSType_AlphaDiv.pdf"),device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results","Figures","16S_Water_TimeSType_AlphaDiv.png"),device='png',height = 7.5, width = 10.5)



###Rarefied abundance ####

# Rarefaction curves
rar_obj<- as.data.frame(t(otu_table(ps)))
rarecurve(rar_obj, step=100, cex=0.5)


### Number of reads per sample
min(sample_sums(ps)) #4780
max(sample_sums(ps)) #75408

readsumsdf = data.frame(nreads = sample_sums(ps), 
                        type = "Samples")

sdat <- data.frame(sample_data(ps))

reads <- cbind(readsumsdf,sdat)

ggplot(reads, aes(x = rownames(reads), y = nreads,fill=Time))+
  theme_bw()+
  geom_bar(stat = "identity",color='black',size=0.2)+
  facet_grid(~Time,drop=T,scales="free")+
  theme(axis.text.x =element_text(angle=45,size=5))+
  scale_y_continuous(expand = expansion(c(0, 0.05)))+
  scale_fill_viridis(option="magma",discrete=T,name="Time")

sdat <- data.frame(sample_data(ps))
sum_sdat <- summary(sdat)
data.frame(sum_sdat[1:5,6])
# Not the same number of samples per sampling date


### Rarefaction 
set.seed(12)
ps.rar <- rarefy_even_depth(ps) #
min(sample_sums(ps.rar))
max(sample_sums(ps.rar))
# Rarefied to 4780 reads, 19 OTUs removed


#####Sample type ####
alpha.div<- plot_richness(ps.rar, x="Sample_type", color="Sample_type", measures=c("Observed","Simpson", "Shannon"))

plot <- select(alpha.div$data,c("samples","Time","Sample_type","Temperature","variable","value"))


ggplot(plot,aes(x=variable,y=value,shape=Sample_type))+
  geom_boxplot(outlier.shape = NA)+
  #geom_boxplot(outlier.shape = NA,position=position_dodge(width=1.5))+
  geom_point(aes(x=variable,y=value,shape=Sample_type),position=position_dodge(width=0.75),size=2.5,fill='black')+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank())+
  stat_pwc(method="wilcox_test",hide.ns=T,
           p.adjust.method="fdr",label="p.adj.signif",
           tip.length = 0, 
           step.increase = 0.05,
           vjust=0)+
  scale_shape_manual(values=c(19,22,23),name="Sample type")+
  facet_wrap(~variable,scales="free")+
  ylab("Alpha diversity score")

# No significant change


#####Temperature ####
alpha.div<- plot_richness(ps.rar, x="Temperature", color="Temperature", measures=c("Observed","Simpson", "Shannon"))

plot <- select(alpha.div$data,c("samples","Time","Sample_type","Temperature","variable","value"))



ggplot(plot,aes(x=variable,y=value,color=Temperature))+
  geom_boxplot(outlier.shape = NA)+
  #geom_boxplot(outlier.shape = NA,position=position_dodge(width=1.5))+
  geom_point(aes(x=variable,y=value),position=position_dodge(width=0.75),size=2.5,fill='black')+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank())+
  stat_pwc(method="wilcox_test",
           hide.ns=T,
           p.adjust.method="fdr",
           label="p.adj.signif",
           tip.length = 0, 
           step.increase = 0.5,
           vjust=0)+
  scale_color_manual(values=c("#0000FE","#B02223"),name="Temperature")+
  facet_wrap(~variable,scales="free")+
  ylab("Alpha diversity score")

# No significant change

#####Time #####
library(ggh4x)
alpha.div <- plot_richness(ps.rar, x="Time", measures=c("Observed","Simpson", "Shannon"))

plot <- select(alpha.div$data,c("samples","Time","Sample_type","Temperature","variable","value"))


ggplot(plot,aes(x=variable,y=value,shape=Sample_type,color=Time))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(x=variable,y=value,shape=Sample_type),position=position_dodge(width=0.75),size=2,fill='grey60')+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank())+
  stat_pwc(method="wilcox_test",
           p.adjust.method="fdr",
           hide.ns=T,
           label="p.adj.signif",
           tip.length = 0, 
           vjust=0.7,
           size = 0.5,
           label.size=4.2)+
  scale_colour_viridis(option="magma",discrete=T,name="Time")+
  scale_shape_manual(values=c(21,22,23),name="Sample type")+
  facet_nested(variable~Sample_type,scales="free",drop=T,independent = "all")+
  ylab("Alpha diversity score")

# Some changes occur 
# (new significant differences or changes in pvalue significance)



