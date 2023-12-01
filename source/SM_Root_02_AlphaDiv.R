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
sample_data(ps)

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


####~~~~~~~~~~~~~~~~####
### ALPHA DIV ####
#### Raw abundance ####

#####Plant type ####
######| Warm ####
alpha.div.warm<- plot_richness(ps_root_Warm, x="Plant.type", measures=c("Observed","Simpson", "Shannon"))

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
    scale_shape_manual(values=c(22,23),name="Plant type :",
                       labels=c(expression(paste(italic("Scirpus"))),expression(italic("Triglochin"))))+
    facet_wrap(~variable,scales="free")+
    ylab("Alpha diversity score")

# Save generated figure
# ggsave(here("Results_W&C","Figures","root_Warm_SType_AlphaDiv.pdf"), plot = Warm_PlantType ,device='pdf',height = 7.5, width = 10.5)
# ggsave(here("Results_W&C","Figures","root_Warm_SType_AlphaDiv.png"), plot = Warm_PlantType, device='png',height = 7.5, width = 10.5)

######| Cold ####
alpha.div.Cold<- plot_richness(ps_root_Cold, x="Plant.type", measures=c("Observed","Simpson", "Shannon"))

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
    scale_shape_manual(values=c(22,23),name="Plant type :",
                       labels=c(expression(paste(italic("Scirpus"))),expression(italic("Triglochin"))))+
    facet_wrap(~variable,scales="free")+
    ylab("Alpha diversity score")

# Save generated figure
# ggsave(here("Results_W&C","Figures","root_Cold_SType_AlphaDiv.pdf"), plot = Cold_PlantType ,device='pdf',height = 7.5, width = 10.5)
# ggsave(here("Results_W&C","Figures","root_Cold_SType_AlphaDiv.png"), plot = Cold_PlantType, device='png',height = 7.5, width = 10.5)

Plant_type_alpha <- ggarrange(Warm_PlantType,Cold_PlantType,
                              nrow = 2,common.legend = T,
                              legend = "bottom",
                              labels = c("Warm","Cold"),
                              hjust = -0.1)

ggsave(here("Results_W&C","Figures","Root_SType_AlphaDiv.pdf"), plot = Plant_type_alpha ,device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results_W&C","Figures","Root_SType_AlphaDiv.png"), plot = Plant_type_alpha, device='png',height = 7.5, width = 10.5)


#####Time ####
######| Warm ####
alpha.div.warm<- plot_richness(ps_root_Warm, x="Time", measures=c("Observed","Simpson", "Shannon"))

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
# ggsave(here("Results_W&C","Figures","root_Warm_Time_AlphaDiv.pdf"), plot = Warm_PlantType ,device='pdf',height = 7.5, width = 10.5)
# ggsave(here("Results_W&C","Figures","root_Warm_Time_AlphaDiv.png"), plot = Warm_PlantType, device='png',height = 7.5, width = 10.5)

######| Cold ####
alpha.div.Cold<- plot_richness(ps_root_Cold, x="Time", measures=c("Observed","Simpson", "Shannon"))

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
# ggsave(here("Results_W&C","Figures","root_Cold_Time_AlphaDiv.pdf"), plot = Cold_PlantType ,device='pdf',height = 7.5, width = 10.5)
# ggsave(here("Results_W&C","Figures","root_Cold_Time_AlphaDiv.png"), plot = Cold_PlantType, device='png',height = 7.5, width = 10.5)

Time_alpha <- ggarrange(Warm_Time,Cold_Time,
                        nrow = 2,common.legend = T,
                        legend = "bottom",
                        labels = c("Warm","Cold"),
                        hjust = -0.1)

ggsave(here("Results_W&C","Figures","Root_Time_AlphaDiv.pdf"), 
       plot = Time_alpha ,
       device='pdf',height = 7.5, width = 10.5)

ggsave(here("Results_W&C","Figures","Root_Time_AlphaDiv.png"), 
       plot = Time_alpha, 
       device='png',height = 7.5, width = 10.5)


##### Time & Plant type####
######| Warm ####
alpha.div.warm<- plot_richness(ps_root_Warm, x="Time", measures=c("Observed","Simpson", "Shannon"))

# Simplified data for plotting
plot <- select(alpha.div.warm$data,c("samples","Time","Plant.type","variable","value"))


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
    facet_nested(variable~Plant.type,scales="free",
                 drop=T,independent = "all",shrink = F)+
    ylab("Alpha diversity score")

# Save generated figure
# ggsave(here("Results_W&C","Figures","root_Warm_PT_Time_AlphaDiv.pdf"), plot = Warm_PlantType ,device='pdf',height = 7.5, width = 10.5)
# ggsave(here("Results_W&C","Figures","root_Warm_PT_Time_AlphaDiv.png"), plot = Warm_PlantType, device='png',height = 7.5, width = 10.5)

######| Cold ####
alpha.div.Cold<- plot_richness(ps_root_Cold, x="Time", measures=c("Observed","Simpson", "Shannon"))

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
                 drop=T,independent = "all",shrink = F)+
    ylab("Alpha diversity score")

# Save generated figure
# ggsave(here("Results_W&C","Figures","root_Cold_PT_Time_AlphaDiv.pdf"), plot = Cold_PlantType ,device='pdf',height = 7.5, width = 10.5)
# ggsave(here("Results_W&C","Figures","root_Cold_PT_Time_AlphaDiv.png"), plot = Cold_PlantType, device='png',height = 7.5, width = 10.5)

Time_PT_alpha <- ggarrange(Warm_PT_Time,Cold_PT_Time,
                           nrow = 2,common.legend = T,
                           legend = "bottom",
                           labels = c("Warm","Cold"),
                           hjust = -0.1)

ggsave(here("Results_W&C","Figures","Root_PT_Time_AlphaDiv.pdf"), 
       plot = Time_PT_alpha ,
       device='pdf',height = 7.5, width = 10.5)

ggsave(here("Results_W&C","Figures","Root_PT_Time_AlphaDiv.png"), 
       plot = Time_PT_alpha, 
       device='png',height = 7.5, width = 10.5)































####~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~####

#### Not separated by temperature #####
####~~~~~~~~~~~~~~~~####
#### Data ####

load(here("Rdata","ps_obj.RData"))

ps_root <-  subset_samples(ps, Sample.type =="Roots")
ps_root <- prune_taxa(taxa_sums(ps_root)>0, ps_root)
ps_root
# 81 samples 5954 taxa


###______#####
### ALPHA DIV ####
### *** Raw abundance ####
alpha.div<- plot_richness(ps_root, x="Plant.type", color="Plant.type", measures=c("Observed","Simpson", "Shannon"))

plot <- select(alpha.div$data,c("samples","Time","Plant.type","Temperature","variable","value"))
plot$Time <- as.factor(plot$Time)


ggplot(plot,aes(x=variable,y=value,shape=Plant.type))+
  geom_boxplot(outlier.shape = NA)+
  #geom_boxplot(outlier.shape = NA,position=position_dodge(width=1.5))+
  geom_point(aes(x=variable,y=value,shape=Plant.type),position=position_dodge(width=0.75),size=2.5,fill='black')+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank())+
  stat_pwc(method="t_test",hide.ns=F,
           p.adjust.method="BH",label="p.adj.signif",
           tip.length = 0, 
           step.increase = 0.5,
           vjust=0)+
  scale_shape_manual(values=c(22,23),name="Plant")+
  facet_wrap(~variable,scales="free")+
  ylab("Alpha diversity score")

ggsave(here("Results","Figures","Rh_AlphaDiv.pdf"),device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results","Figures","Rh_AlphaDiv.png"),device='png',height = 7.5, width = 10.5)


