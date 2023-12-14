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


####~~~~~~~~~~~~~~~~####
### ALPHA DIV ####
#### Raw abundance ####

#####Plant type ####
######| Warm ####
alpha.div.warm<- plot_richness(ps_rhizo_Warm, x="Plant.type", measures=c("Observed","Simpson", "Shannon"))

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
# ggsave(here("Results_W&C","Figures","Rhizo_Warm_SType_AlphaDiv.pdf"), plot = Warm_PlantType ,device='pdf',height = 7.5, width = 10.5)
# ggsave(here("Results_W&C","Figures","Rhizo_Warm_SType_AlphaDiv.png"), plot = Warm_PlantType, device='png',height = 7.5, width = 10.5)

######| Cold ####
alpha.div.Cold<- plot_richness(ps_rhizo_Cold, x="Plant.type", measures=c("Observed","Simpson", "Shannon"))

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
# ggsave(here("Results_W&C","Figures","Rhizo_Cold_SType_AlphaDiv.pdf"), plot = Cold_PlantType ,device='pdf',height = 7.5, width = 10.5)
# ggsave(here("Results_W&C","Figures","Rhizo_Cold_SType_AlphaDiv.png"), plot = Cold_PlantType, device='png',height = 7.5, width = 10.5)



### Check significant resutls with rarefied data ##

# sample_sums(ps_rhizo_Cold)
# # One sample (69.70MesoRh.11G.1.E2.16S.D13) 
# # has a low count number compared to the rest 
# # remove it and its unique ASVs before rarefying 
# 
# ps_rhizo_Cold_sub<- subset_samples(ps_rhizo_Cold, sample_sums(ps_rhizo_Cold)>5000)
# ps_rhizo_Cold_sub <- prune_taxa(taxa_sums(ps_rhizo_Cold_sub)>0, ps_rhizo_Cold_sub)
# 
# ps_rhizo_Cold_rar <- rarefy_even_depth(ps_rhizo_Cold_sub)
# 
# alpha.div.Cold.rar<- plot_richness(ps_rhizo_Cold_rar, x="Plant.type", measures=c("Observed","Simpson", "Shannon"))
# 
# # Simplified data for plotting
# plot <- select(alpha.div.Cold.rar$data,c("samples","Time","Plant.type","variable","value"))
# 
# Cold_PlantType_rar <- ggplot(plot,aes(x=variable,y=value,shape=Plant.type))+
#     geom_boxplot(outlier.shape = NA)+
#     geom_point(aes(x=variable,y=value,shape=Plant.type),position=position_dodge(width=0.75),size=2.5,fill='black')+
#     theme_bw()+
#     theme(title = element_text(size=20),
#           axis.text.x =element_blank(),
#           axis.text.y =element_text(size=12),
#           axis.title.x = element_blank(), 
#           axis.title.y = element_text(size=15),
#           legend.text.align = 0,
#           axis.ticks.x = element_blank(),
#           legend.text=element_text(size=20),
#           legend.position = "bottom",
#           strip.text = element_text(size=15),
#           strip.background = element_rect(fill='white'))+
#     stat_pwc(method="wilcox_test",hide.ns=T,
#              p.adjust.method="fdr",
#              label="p.adj.signif",
#              tip.length = 0, 
#              step.increase = 0.05,
#              vjust=0)+
#     # Computes pairwise wilcoxon rank-sum test in each facet
#     scale_shape_manual(values=c(22,23),name="Plant type :",
#                        labels=c(expression(paste(italic("Scirpus"))),expression(italic("Triglochin"))))+
#     facet_wrap(~variable,scales="free")+
#     ylab("Alpha diversity score")


### Common plot ###
Plant_type_alpha <- ggarrange(Warm_PlantType,Cold_PlantType,
                              nrow = 2,common.legend = T,
                              legend = "bottom",
                              labels = c("Warm","Cold"),
                              hjust = -0.1)

ggsave(here("Results_W&C","Figures","Rhizo_SType_AlphaDiv.pdf"), plot = Plant_type_alpha ,device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results_W&C","Figures","Rhizo_SType_AlphaDiv.png"), plot = Plant_type_alpha, device='png',height = 7.5, width = 10.5)


#####Time ####
######| Warm ####
alpha.div.warm<- plot_richness(ps_rhizo_Warm, x="Time", measures=c("Observed","Simpson", "Shannon"))

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
# ggsave(here("Results_W&C","Figures","Rhizo_Warm_Time_AlphaDiv.pdf"), plot = Warm_PlantType ,device='pdf',height = 7.5, width = 10.5)
# ggsave(here("Results_W&C","Figures","Rhizo_Warm_Time_AlphaDiv.png"), plot = Warm_PlantType, device='png',height = 7.5, width = 10.5)

######| Cold ####
alpha.div.Cold<- plot_richness(ps_rhizo_Cold, x="Time", measures=c("Observed","Simpson", "Shannon"))

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
# ggsave(here("Results_W&C","Figures","Rhizo_Cold_Time_AlphaDiv.pdf"), plot = Cold_PlantType ,device='pdf',height = 7.5, width = 10.5)
# ggsave(here("Results_W&C","Figures","Rhizo_Cold_Time_AlphaDiv.png"), plot = Cold_PlantType, device='png',height = 7.5, width = 10.5)

Time_alpha <- ggarrange(Warm_Time,Cold_Time,
                        nrow = 2,common.legend = T,
                        legend = "bottom",
                        labels = c("Warm","Cold"),
                        hjust = -0.1)

ggsave(here("Results_W&C","Figures","Rhizo_Time_AlphaDiv.pdf"), 
       plot = Time_alpha ,
       device='pdf',height = 7.5, width = 10.5)

ggsave(here("Results_W&C","Figures","Rhizo_Time_AlphaDiv.png"), 
       plot = Time_alpha, 
       device='png',height = 7.5, width = 10.5)


##### Time & Plant type####
######| Warm ####
alpha.div.warm<- plot_richness(ps_rhizo_Warm, x="Time", measures=c("Observed","Simpson", "Shannon"))

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
# ggsave(here("Results_W&C","Figures","Rhizo_Warm_PT_Time_AlphaDiv.pdf"), plot = Warm_PlantType ,device='pdf',height = 7.5, width = 10.5)
# ggsave(here("Results_W&C","Figures","Rhizo_Warm_PT_Time_AlphaDiv.png"), plot = Warm_PlantType, device='png',height = 7.5, width = 10.5)

######| Cold ####
alpha.div.Cold<- plot_richness(ps_rhizo_Cold, x="Time", measures=c("Observed","Simpson", "Shannon"))

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
                 drop=T,independent = "all",shrink = F,
                 labeller = labeller(Plant.type=Plant_labs))+
    ylab("Alpha diversity score")

# Save generated figure
# ggsave(here("Results_W&C","Figures","rhizo_Cold_PT_Time_AlphaDiv.pdf"), plot = Cold_PlantType ,device='pdf',height = 7.5, width = 10.5)
# ggsave(here("Results_W&C","Figures","rhizo_Cold_PT_Time_AlphaDiv.png"), plot = Cold_PlantType, device='png',height = 7.5, width = 10.5)

Time_PT_alpha <- ggarrange(Warm_PT_Time,Cold_PT_Time,
                           nrow = 2,common.legend = T,
                           legend = "bottom",
                           labels = c("Warm","Cold"),
                           hjust = -0.1)

ggsave(here("Results_W&C","Figures","Rhizo_PT_Time_AlphaDiv.pdf"), 
       plot = Time_PT_alpha ,
       device='pdf',height = 7.5, width = 10.5)

ggsave(here("Results_W&C","Figures","Rhizo_PT_Time_AlphaDiv.png"), 
       plot = Time_PT_alpha, 
       device='png',height = 7.5, width = 10.5)


#### Supp Figure ######

######| Warm ####
alpha.div.warm<- plot_richness(ps_rhizo_Warm, x="Time", measures= "Shannon")

# Simplified data for plotting
plot <- select(alpha.div.warm$data,c("samples","Time","Plant.type","variable","value"))

Warm <- ggplot(plot,aes(x=Time,y=value,fill=Time))+
    geom_point(aes(x=Time,y=value,shape=Plant.type,fill=Time),size=4,position=position_dodge(width=0.5))+
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
    scale_shape_manual(values=c(22,23),name="Plant type :",
                       labels=c(expression(paste(italic("Scirpus"))),expression(italic("Triglochin"))))+
    ylab("Alpha diversity score")
Warm



######| Cold ####
alpha.div.cold<- plot_richness(ps_rhizo_Cold, x="Time", measures= "Shannon")

# Simplified data for plotting
plot <- select(alpha.div.cold$data,c("samples","Time","Plant.type","variable","value"))

Cold <- ggplot(plot,aes(x=Time,y=value,fill=Time))+
    geom_point(aes(x=Time,y=value,shape=Plant.type,fill=Time),size=4,position=position_dodge(width=0.5))+
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
    scale_shape_manual(values=c(22,23),name="Plant type :",
                       labels=c(expression(paste(italic("Scirpus"))),expression(italic("Triglochin"))))+
    ylab("Alpha diversity score")
Cold

Alpha <- ggarrange(Warm,Cold,
                   ncol = 2,common.legend = T,
                   legend = "bottom", labels=c("A","B"))


ggsave(here("Results_W&C","Figures","SupFig_RhizoAlphaDiv.pdf"), 
       plot = Alpha ,
       device='pdf',height = 7.5, width = 10.5)

ggsave(here("Results_W&C","Figures","SupFig_RhizoAlphaDiv.png"), 
       plot = Alpha, 
       device='png',height = 7.5, width = 10.5)






























####~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~####

#### Not separated by temperature #####
####~~~~~~~~~~~~~~~~####
#### Data ####

load(here("Rdata","ps_obj.RData"))

ps_rhizo <-  subset_samples(ps, Sample.type =="Rhizosphere")
ps_rhizo <- prune_taxa(taxa_sums(ps_rhizo)>0, ps_rhizo)
ps_rhizo
# 81 samples 5954 taxa


###______#####
### ALPHA DIV ####
### *** Raw abundance ####
alpha.div<- plot_richness(ps_rhizo, x="Plant.type", color="Plant.type", measures=c("Observed","Simpson", "Shannon"))

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


