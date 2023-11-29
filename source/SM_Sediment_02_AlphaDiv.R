### Packages#####
library(here)
source(here("source", "libraries.R"))
library(phyloseq)
library(viridis)
library(metacoder)
library(dplyr)
library(ggpp)
library(ggh4x)
#### Data ####

load(here("Rdata","ps_obj.RData"))

#### Separation by temperature #####

ps_sed_Warm <- subset_samples(ps, Sample.type =="Sediment"& Temperature =="20°C day 10°C night")
ps_sed_Warm <- prune_taxa(taxa_sums(ps_sed_Warm)>0, ps_sed_Warm)
ps_sed_Warm # 61 samples - 11912 taxa
# check to see if you kept only SEDIMENT samples in WARM condition
str(sample_data(ps_sed_Warm)) 

ps_sed_Cold <- subset_samples(ps, Sample.type =="Sediment"& Temperature =="10°C day 5°C night")
ps_sed_Cold <- prune_taxa(taxa_sums(ps_sed_Cold)>0, ps_sed_Cold)
ps_sed_Cold # 59 samples - 12973 taxa
# check to see if you kept only SEDIMENT samples in COLD condition
str(sample_data(ps_sed_Cold)) 

####~~~~~~~~~~~~~~~~####
### ALPHA DIV ####
#### Raw abundance ####

#####Plant type ####
######| Warm ####
alpha.div.warm<- plot_richness(ps_sed_Warm, x="Plant.type", measures=c("Observed","Simpson", "Shannon"))

# Simplified data for plotting
plot <- select(alpha.div.warm$data,c("samples","Time","Plant.type","variable","value"))

Warm_PlantType <- ggplot(plot,aes(x=variable,y=value,shape=Plant.type))+
                geom_boxplot(outlier.shape = NA)+
                geom_point(aes(x=variable,y=value,shape=Plant.type),position=position_dodge(width=0.75),size=2.5,fill='black')+
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
# ggsave(here("Results_W&C","Figures","Sed_Warm_SType_AlphaDiv.pdf"), plot = Warm_PlantType ,device='pdf',height = 7.5, width = 10.5)
# ggsave(here("Results_W&C","Figures","Sed_Warm_SType_AlphaDiv.png"), plot = Warm_PlantType, device='png',height = 7.5, width = 10.5)

######| Cold ####
alpha.div.Cold<- plot_richness(ps_sed_Cold, x="Plant.type", measures=c("Observed","Simpson", "Shannon"))

# Simplified data for plotting
plot <- select(alpha.div.Cold$data,c("samples","Time","Plant.type","variable","value"))

Cold_PlantType <- ggplot(plot,aes(x=variable,y=value,shape=Plant.type))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(aes(x=variable,y=value,shape=Plant.type),position=position_dodge(width=0.75),size=2.5,fill='black')+
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
# ggsave(here("Results_W&C","Figures","Sed_Cold_SType_AlphaDiv.pdf"), plot = Cold_PlantType ,device='pdf',height = 7.5, width = 10.5)
# ggsave(here("Results_W&C","Figures","Sed_Cold_SType_AlphaDiv.png"), plot = Cold_PlantType, device='png',height = 7.5, width = 10.5)

Plant_type_alpha <- ggarrange(Warm_PlantType,Cold_PlantType,
                              nrow = 2,common.legend = T,
                              legend = "bottom",
                              labels = c("Warm","Cold"),
                              hjust = -0.1)

ggsave(here("Results_W&C","Figures","Sed_SType_AlphaDiv.pdf"), plot = Plant_type_alpha ,device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results_W&C","Figures","Sed_SType_AlphaDiv.png"), plot = Plant_type_alpha, device='png',height = 7.5, width = 10.5)


#####Time ####
######| Warm ####
alpha.div.warm<- plot_richness(ps_sed_Warm, x="Time", measures=c("Observed","Simpson", "Shannon"))

# Simplified data for plotting
plot <- select(alpha.div.warm$data,c("samples","Time","Plant.type","variable","value"))

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
# ggsave(here("Results_W&C","Figures","Sed_Warm_Time_AlphaDiv.pdf"), plot = Warm_PlantType ,device='pdf',height = 7.5, width = 10.5)
# ggsave(here("Results_W&C","Figures","Sed_Warm_Time_AlphaDiv.png"), plot = Warm_PlantType, device='png',height = 7.5, width = 10.5)

######| Cold ####
alpha.div.Cold<- plot_richness(ps_sed_Cold, x="Time", measures=c("Observed","Simpson", "Shannon"))

# Simplified data for plotting
plot <- select(alpha.div.Cold$data,c("samples","Time","Plant.type","variable","value"))

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
# ggsave(here("Results_W&C","Figures","Sed_Cold_Time_AlphaDiv.pdf"), plot = Cold_PlantType ,device='pdf',height = 7.5, width = 10.5)
# ggsave(here("Results_W&C","Figures","Sed_Cold_Time_AlphaDiv.png"), plot = Cold_PlantType, device='png',height = 7.5, width = 10.5)

Time_alpha <- ggarrange(Warm_Time,Cold_Time,
                              nrow = 2,common.legend = T,
                              legend = "bottom",
                              labels = c("Warm","Cold"),
                              hjust = -0.1)

ggsave(here("Results_W&C","Figures","Sed_Time_AlphaDiv.pdf"), 
       plot = Time_alpha ,
       device='pdf',height = 7.5, width = 10.5)

ggsave(here("Results_W&C","Figures","Sed_Time_AlphaDiv.png"), 
       plot = Time_alpha, 
       device='png',height = 7.5, width = 10.5)


##### Time & Plant type####
######| Warm ####
alpha.div.warm<- plot_richness(ps_sed_Warm, x="Time", measures=c("Observed","Simpson", "Shannon"))

# Simplified data for plotting
plot <- select(alpha.div.warm$data,c("samples","Time","Plant.type","variable","value"))

# Rename plant labels
Plant_labs <- c("No Plant","Scirpus", "Triglochin")
names(Plant_labs) <- c("no plant","Scirpus", "Triglochin")


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
    facet_nested(variable~Plant.type,scales="free",
                 drop=T,independent = "all",shrink = F,
                 labeller = labeller(Plant.type=Plant_labs))+
    ylab("Alpha diversity score")

# Save generated figure
# ggsave(here("Results_W&C","Figures","Sed_Warm_PT_Time_AlphaDiv.pdf"), plot = Warm_PlantType ,device='pdf',height = 7.5, width = 10.5)
# ggsave(here("Results_W&C","Figures","Sed_Warm_PT_Time_AlphaDiv.png"), plot = Warm_PlantType, device='png',height = 7.5, width = 10.5)

######| Cold ####
alpha.div.Cold<- plot_richness(ps_sed_Cold, x="Time", measures=c("Observed","Simpson", "Shannon"))

# Simplified data for plotting
plot <- select(alpha.div.Cold$data,c("samples","Time","Plant.type","variable","value"))

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
    facet_nested(variable~Plant.type,scales="free",
                 drop=T,independent = "all",shrink = F,
                 labeller = labeller(Plant.type=Plant_labs))+
    ylab("Alpha diversity score")

# Save generated figure
# ggsave(here("Results_W&C","Figures","Sed_Cold_PT_Time_AlphaDiv.pdf"), plot = Cold_PlantType ,device='pdf',height = 7.5, width = 10.5)
# ggsave(here("Results_W&C","Figures","Sed_Cold_PT_Time_AlphaDiv.png"), plot = Cold_PlantType, device='png',height = 7.5, width = 10.5)

Time_PT_alpha <- ggarrange(Warm_PT_Time,Cold_PT_Time,
                        nrow = 2,common.legend = T,
                        legend = "bottom",
                        labels = c("Warm","Cold"),
                        hjust = -0.1)

ggsave(here("Results_W&C","Figures","Sed_PT_Time_AlphaDiv.pdf"), 
       plot = Time_PT_alpha ,
       device='pdf',height = 7.5, width = 10.5)

ggsave(here("Results_W&C","Figures","Sed_PT_Time_AlphaDiv.png"), 
       plot = Time_PT_alpha, 
       device='png',height = 7.5, width = 10.5)








### Rarefied Abundance ####

## Warm ##
# Abundance extremes
min(sample_sums(ps_sed_Warm))
max(sample_sums(ps_sed_Warm))

# Object with read number per sample
readsumsdf = data.frame(nreads = sample_sums(ps_sed_Warm), 
                        type = "Samples")

# Add sample data to the object
sdat <- as.data.frame(sample_data(ps_sed_Warm))
reads <- cbind(readsumsdf,sdat)

#Plot read number per sample grouped by sampling date
ggplot(reads, aes(x = rownames(reads), y = nreads,fill=Time))+
    theme_bw()+
    geom_bar(stat = "identity",color='black',size=0.2)+
    facet_grid(. ~ Time,scales="free_x",drop=T)+
    theme(axis.text.x =element_text(angle=45,size=5))+
    scale_y_continuous(expand = expansion(c(0, 0.05)))+
    scale_fill_viridis(option="magma",discrete=T,name="Time")


### Rarefaction 
ps_sed_Warm_sub <- prune_samples(sample_sums(ps_sed_Warm)>10000,ps_sed_Warm)
ps_sed_Warm_sub <- prune_taxa(taxa_sums(ps_sed_Warm_sub)>0, ps_sed_Warm_sub)
# Loss of 4 D0 samples and 27 ASVs were only present in those samples

set.seed(12)
ps.Warm.rar <- rarefy_even_depth(ps_sed_Warm_sub) # 88 ASVs lost
min(sample_sums(ps.Warm.rar)) #20816
max(sample_sums(ps.Warm.rar)) #20816


## Cold ##
# Abundance extremes
min(sample_sums(ps_sed_Cold))
max(sample_sums(ps_sed_Cold))

# Object with read number per sample
readsumsdf = data.frame(nreads = sample_sums(ps_sed_Cold), 
                        type = "Samples")

# Add sample data to the object
sdat <- as.data.frame(sample_data(ps_sed_Cold))
reads <- cbind(readsumsdf,sdat)

#Plot read number per sample grouped by sampling date
ggplot(reads, aes(x = rownames(reads), y = nreads,fill=Time))+
    theme_bw()+
    geom_bar(stat = "identity",color='black',size=0.2)+
    facet_grid(. ~ Time,scales="free_x",drop=T)+
    theme(axis.text.x =element_text(angle=45,size=5))+
    scale_y_continuous(expand = expansion(c(0, 0.05)))+
    scale_fill_viridis(option="magma",discrete=T,name="Time")


### Rarefaction 
ps_sed_Cold_sub <- prune_samples(sample_sums(ps_sed_Cold)>10000,ps_sed_Cold)
ps_sed_Cold_sub <- prune_taxa(taxa_sums(ps_sed_Cold_sub)>0, ps_sed_Cold_sub)
# Loss of 1 D0 sample and 0 ASV only present in that sample

set.seed(12)
ps.Cold.rar <- rarefy_even_depth(ps_sed_Cold_sub) # 54 ASVs lost
min(sample_sums(ps.Cold.rar)) #22709
max(sample_sums(ps.Cold.rar)) #22709


#####Plant type ####
######| Warm ####

alpha.div.Warm.rar<- plot_richness(ps.Warm.rar, x="Plant.type", measures=c("Observed","Simpson", "Shannon"))

# Simplified data for plotting
plot <- select(alpha.div.Warm.rar$data,c("samples","Time","Plant.type","variable","value"))

Warm_PlantType_rar <- ggplot(plot,aes(x=variable,y=value,shape=Plant.type))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(aes(x=variable,y=value,shape=Plant.type),position=position_dodge(width=0.75),size=2.5,fill='black')+
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


######| Cold ####
alpha.div.Cold.rar<- plot_richness(ps.Cold.rar, x="Plant.type", measures=c("Observed","Simpson", "Shannon"))

# Simplified data for plotting
plot <- select(alpha.div.Cold.rar$data,c("samples","Time","Plant.type","variable","value"))

Cold_PlantType_rar <- ggplot(plot,aes(x=variable,y=value,shape=Plant.type))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(aes(x=variable,y=value,shape=Plant.type),position=position_dodge(width=0.75),size=2.5,fill='black')+
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


Plant_type_alpha_rar <- ggarrange(Warm_PlantType_rar,Cold_PlantType_rar,
                              nrow = 2,common.legend = T,
                              legend = "bottom",
                              labels = c("Warm","Cold"),
                              hjust = -0.1)


#  Change in Cold Shannon 
#  Triglo no longer signif diff from Scirpus and unplanted 
#  compared to unrarefied

#####Time ####
######| Warm ####
alpha.div.Warm.rar<- plot_richness(ps.Warm.rar, x="Time", measures=c("Observed","Simpson", "Shannon"))

# Simplified data for plotting
plot <- select(alpha.div.Warm.rar$data,c("samples","Time","Plant.type","variable","value"))

Warm_Time_rar <- ggplot(plot,aes(x=variable,y=value,color=Time))+
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

######| Cold ####
alpha.div.Cold.rar<- plot_richness(ps.Cold.rar, x="Time", measures=c("Observed","Simpson", "Shannon"))

# Simplified data for plotting
plot <- select(alpha.div.Cold.rar$data,c("samples","Time","Plant.type","variable","value"))

Cold_Time_rar <- ggplot(plot,aes(x=variable,y=value,color=Time))+
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

Time_alpha_rar <- ggarrange(Warm_Time_rar,Cold_Time_rar,
                        nrow = 2,common.legend = T,
                        legend = "bottom",
                        labels = c("Warm","Cold"),
                        hjust = -0.1)

# Some changes in Warm Observed  
# but D62 still signif. diff from the other dates

####~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~####

#### Not separated by temperature #####
####~~~~~~~~~~~~~~~~####

# Keep only sediment samples and remove ASVs with 0 reads
# ps_sed <-  subset_samples(ps, Sample.type =="Sediment")
# ps_sed <- prune_taxa(taxa_sums(ps_sed)>0, ps_sed)
# ps_sed
# 120 samples 18851 taxa

# Simplify sample names 
#sample_names(ps_sed) <- sub("Meso.*","",sample_names(ps_sed))



### ALPHA DIV ####
##### Raw abundance ####

######Plant type ####
alpha.div<- plot_richness(ps_sed, x="Plant.type", measures=c("Observed","Simpson", "Shannon"))

# Simplified data for plotting
plot <- select(alpha.div$data,c("samples","Time","Plant.type","Temperature","variable","value"))

ggplot(plot,aes(x=variable,y=value,shape=Plant.type))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(aes(x=variable,y=value,shape=Plant.type),position=position_dodge(width=0.75),size=2.5,fill='black')+
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
                       labels=c("No plant",expression(paste("ab",italic("Scirpus"))),expression(italic("Triglochin"))))+
    facet_wrap(~variable,scales="free")+
    ylab("Alpha diversity score")

# Save generated figure
ggsave(here("Results","Figures","Sed_SType_AlphaDiv.pdf"),device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results","Figures","Sed_SType_AlphaDiv.png"),device='png',height = 7.5, width = 10.5)


######Temperature ####
alpha.div<- plot_richness(ps, x="Temperature", measures=c("Observed","Simpson", "Shannon"))


# Simplified data for plotting
plot <- select(alpha.div$data,c("samples","Time","Plant.type","Temperature","variable","value"))

ggplot(plot,aes(x=variable,y=value,color=Temperature))+
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
             step.increase = 0.5,
             vjust=0)+
    # Computes pairwise wilcoxon rank-sum test in each facet
    scale_color_manual(values=c("#0000FE","#B02223"),name="Temperature")+
    facet_wrap(~variable,scales="free")+
    ylab("Alpha diversity score")

# Save generated figure
ggsave(here("Results","Figures","Sed_Temp_AlphaDiv.pdf"),device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results","Figures","Sed_Temp_AlphaDiv.png"),device='png',height = 7.5, width = 10.5)


###### Time ####

alpha.div<- plot_richness(ps_sed, x="Time", measures=c("Observed","Simpson", "Shannon"))

# Simplied object for plotting
plot <- select(alpha.div$data,c("samples","Time","Plant.type","Temperature","variable","value"))

# Rename plant labels
Plant_labs <- c("No Plant","Scirpus", "Triglochin")
names(Plant_labs) <- c("no plant","Scirpus", "Triglochin")

ggplot(plot,aes(x=variable,y=value,color=Time))+
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
                  facet_nested(variable~Plant.type,scales="free",
                               drop=T,independent = "all",shrink = F,
                               labeller = labeller(Plant.type=Plant_labs))+
                  ylab("Alpha diversity score")

# Save generated figure
ggsave(here("Results","Figures","Sed_Time_AlphaDiv.pdf"),device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results","Figures","Sed_Time_AlphaDiv.png"),device='png',height = 7.5, width = 10.5)


### *** Rarefied Abundance ####
ps_sed<- subset_samples(ps, Sample.type =="Sediment")

# Rarefaction curves
rar_obj<- as.data.frame(t(otu_table(ps_sed)))
rarecurve(rar_obj, step=100, cex=0.5)

# Abundance extremes
min(sample_sums(ps_sed))
max(sample_sums(ps_sed))

readsumsdf = data.frame(nreads = sample_sums(ps_sed), 
                        type = "Samples")

sdat <- as.data.frame(sample_data(ps_sed))

reads <- cbind(readsumsdf,sdat)

# Read number per sample grouped by sampling date
ggplot(reads, aes(x = rownames(reads), y = nreads,fill=Time))+
    theme_bw()+
    geom_bar(stat = "identity",color='black',size=0.2)+
    facet_wrap(Time ~ Temperature,scales="free_x",drop=T)+
    theme(axis.text.x =element_text(angle=45,size=5))+
    scale_y_continuous(expand = expansion(c(0, 0.05)))+
    scale_fill_viridis(option="magma",discrete=T,name="Time")

# ggsave(here("Results","Figures","Sed_Read#.pdf"),device='pdf',height = 7.5, width = 10.5)


### Rarefaction 
ps_sed_sub <- prune_samples(sample_sums(ps_sed)>10000,ps_sed)
ps_sed_sub <- prune_taxa(taxa_sums(ps_sed_sub)>0, ps_sed_sub)
# Loss of 5 D4 samples (8,15,27,30,31) and 13 ASVs were only present in those samples

set.seed(12)
ps.rar <- rarefy_even_depth(ps_sed_sub)
min(sample_sums(ps.rar))
max(sample_sums(ps.rar))
# Rarefied to 20816 reads,  36 ASVs removed 


######Plant type ####
alpha.div.rar<- plot_richness(ps.rar, x="Plant.type", measures=c("Observed","Simpson", "Shannon"))

# Simplified object for plotting purposes 
plot <- select(alpha.div.rar$data,c("samples","Time","Plant.type","Temperature","variable","value"))


ggplot(plot,aes(x=variable,y=value,shape=Plant.type))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(aes(x=variable,y=value,shape=Plant.type),position=position_dodge(width=0.75),size=2.5,fill='black')+
    theme_bw()+
    theme(title = element_text(size=20),
          axis.text.x =element_blank(),
          axis.text.y =element_text(size=12),
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size=15),
          legend.text.align = 0,
          axis.ticks.x = element_blank(),
          legend.text=element_text(size=20),
          strip.text = element_text(size=15),
          strip.background = element_rect(fill='white'))+
    stat_pwc(method="wilcox_test",hide.ns=T,
             p.adjust.method="fdr",label="p.adj.signif",
             tip.length = 0, 
             step.increase = 0.05,
             vjust=0)+
    scale_shape_manual(values=c(19,22,23),name="Plant type",
                       labels=c("No plant",expression(italic("Scirpus")),expression(italic("Triglochin"))))+
    facet_wrap(~variable,scales="free")+
    ylab("Alpha diversity score")

######Temperature ####
alpha.div.rar <- plot_richness(ps.rar, x="Temperature", measures=c("Observed","Simpson", "Shannon"))

plot <- select(alpha.div.rar$data,c("samples","Time","Plant.type","Temperature","variable","value"))


ggplot(plot,aes(x=variable,y=value,color=Temperature))+
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
             step.increase = 0.5,
             vjust=0)+
    scale_color_manual(values=c("#0000FE","#B02223"),name="Temperature")+
    facet_wrap(~variable,scales="free")+
    ylab("Alpha diversity score")

###### Time ####

alpha.div.rar<- plot_richness(ps_sed, x="Time", color="Time", measures=c("Observed","Simpson", "Shannon"))

plot <- select(alpha.div.rar$data,c("samples","Time","Plant.type","Temperature","variable","value"))

ggplot(plot,aes(x=variable,y=value,color=Time))+
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
             step.increase = 0.05,
             vjust=0)+
    # Computes pairwise wilcoxon rank-sum test in each facet
    scale_colour_viridis(option="magma",discrete=T,name="Time")+
    facet_wrap(~variable,scales="free")+
    ylab("Alpha diversity score")




