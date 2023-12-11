### Packages#####
library(here)
source(here("source", "libraries.R"))
library(phyloseq)
library(viridis)
library(metacoder)
library(dplyr)
library(ggpp)
library(pairwiseAdonis)

##### Data ####

load(here("Rdata","ps_obj.RData"))

##SEDIMENTS #####


##### Separation by temperature #####

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


#### Beta Diversity ####

#### | Warm ####
ps_Warm_rclr <- microbiome::transform(ps_sed_Warm, "rclr") # Robust Aitchison transformation

# Permanova
rclr_dist_matrix <- phyloseq::distance(ps_Warm_rclr, method = "euclidean") # Robust Aitchison distance
metadata <- as(sample_data(ps_Warm_rclr), "data.frame") # Metadata

Permanova.Warm.rclr<-adonis2(rclr_dist_matrix~Plant.type*Time,
                             data = metadata,
                             permutations = 999)

# PCA plots #
ord_Warm_rclr <- phyloseq::ordinate(ps_Warm_rclr, "RDA", distance = "euclidean")

Sed_warm_fig<- phyloseq::plot_ordination(ps_Warm_rclr, ord_Warm_rclr, type="samples",shape='Plant.type',color='Time') + 
    theme_bw()+
    theme(legend.text.align = 0,
          text = element_text(size=17))+
    ggtitle(label = "Sediments - Warm temperatures")+
    geom_point(size=4,aes(fill=sample_data(ps_Warm_rclr)$Time))+
    geom_point(size=4.1,color='black')+ # Workaround to get black outline on previous geom_points 
    scale_shape_manual(values=c(21,22,23),name="Plant type",
                       labels=c("No plant",expression(italic(Scirpus)),expression(italic(Triglochin))))+
    scale_fill_viridis(option="magma",discrete=T,name="Time")+ # Add both color and fill to get the desired graphical output
    scale_colour_viridis(option="magma",discrete=T,name="Time")
# If only fill -> legend is black..  plot_ordination does not have a fill option
Sed_warm_fig$layers <- Sed_warm_fig$layers[-1] 
Sed_warm_fig

#### | Cold ####
ps_Cold_rclr <- microbiome::transform(ps_sed_Cold, "rclr") # Robust Aitchison transformation

# Permanova
rclr_dist_matrix <- phyloseq::distance(ps_Cold_rclr, method = "euclidean") # Robust Aitchison distance
metadata <- as(sample_data(ps_Cold_rclr), "data.frame") # Metadata

Permanova.Cold.rclr<-adonis2(rclr_dist_matrix~Plant.type*Time,
                             data = metadata,
                             permutations = 999)

# PCA plots #
ord_Cold_rclr <- phyloseq::ordinate(ps_Cold_rclr, "RDA", distance = "euclidean")

Sed_cold_fig<- phyloseq::plot_ordination(ps_Cold_rclr, ord_Cold_rclr, type="samples",shape='Plant.type',color='Time')+
    theme_bw()+
    theme(legend.text.align = 0,
          text = element_text(size=17))+
    ggtitle(label = "Sediments - Cold temperatures")+
    geom_point(size=4,aes(fill=sample_data(ps_Cold_rclr)$Time))+
    geom_point(size=4.1,color='black')+ # Workaround to get black outline on previous geom_points 
    scale_shape_manual(values=c(21,22,23),name="Plant type",
                       labels=c("No plant",expression(italic(Scirpus)),expression(italic(Triglochin))))+
    scale_fill_viridis(option="magma",discrete=T,name="Time")+ # Add both color and fill to get the desired graphical output
    scale_colour_viridis(option="magma",discrete=T,name="Time")
# If only fill -> legend is black..  plot_ordination does not have a fill option
Sed_cold_fig$layers <- Sed_cold_fig$layers[-1] 
Sed_cold_fig

Sed_beta_div <- ggarrange(Sed_warm_fig,Sed_cold_fig,ncol=2,common.legend = T,legend = "bottom")

## RHIZOSPHERE #####


##### Separation by temperature #####

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


#### Beta Diversity ####

#### | Warm ####
ps_Warm_rclr <- microbiome::transform(ps_rhizo_Warm, "rclr") # Robust Aitchison transformation

# Permanova
rclr_dist_matrix <- phyloseq::distance(ps_Warm_rclr, method = "euclidean") # Robust Aitchison distance
metadata <- as(sample_data(ps_Warm_rclr), "data.frame") # Metadata

Permanova.Warm.rclr<-adonis2(rclr_dist_matrix~Plant.type*Time,
                             data = metadata,
                             permutations = 999)

# PCA plots #
ord_Warm_rclr <- phyloseq::ordinate(ps_Warm_rclr, "RDA", distance = "euclidean")

rhizo_warm_fig<- phyloseq::plot_ordination(ps_Warm_rclr, ord_Warm_rclr, type="samples",shape='Plant.type',color='Time') + 
    theme_bw()+
    theme(legend.text.align = 0,
          text = element_text(size=17))+
    ggtitle(label = "Rhizosphere - Warm temperatures")+
    geom_point(size=4,aes(fill=sample_data(ps_Warm_rclr)$Time))+
    geom_point(size=4.1,color='black')+ # Workaround to get black outline on previous geom_points 
    scale_shape_manual(values=c(22,23),name="Plant type",
                       labels=c(expression(italic(Scirpus)),expression(italic(Triglochin))))+
    scale_fill_viridis(option="magma",discrete=T,name="Time")+ # Add both color and fill to get the desired graphical output
    scale_colour_viridis(option="magma",discrete=T,name="Time")
# If only fill -> legend is black..  plot_ordination does not have a fill option
rhizo_warm_fig$layers <- rhizo_warm_fig$layers[-1] 
rhizo_warm_fig

#### | Cold ####
ps_Cold_rclr <- microbiome::transform(ps_rhizo_Cold, "rclr") # Robust Aitchison transformation

# Permanova
rclr_dist_matrix <- phyloseq::distance(ps_Cold_rclr, method = "euclidean") # Robust Aitchison distance
metadata <- as(sample_data(ps_Cold_rclr), "data.frame") # Metadata

Permanova.Cold.rclr<-adonis2(rclr_dist_matrix~Plant.type*Time,
                             data = metadata,
                             permutations = 999)

# PCA plots #
ord_Cold_rclr <- phyloseq::ordinate(ps_Cold_rclr, "RDA", distance = "euclidean")

rhizo_cold_fig<- phyloseq::plot_ordination(ps_Cold_rclr, ord_Cold_rclr, type="samples",shape='Plant.type',color='Time')+
    theme_bw()+
    theme(legend.text.align = 0,
          text = element_text(size=17))+
    ggtitle(label = "Rhizosphere - Cold temperatures")+
    geom_point(size=4,aes(fill=sample_data(ps_Cold_rclr)$Time))+
    geom_point(size=4.1,color='black')+ # Workaround to get black outline on previous geom_points 
    scale_shape_manual(values=c(22,23),name="Plant type",
                       labels=c(expression(italic(Scirpus)),expression(italic(Triglochin))))+
    scale_fill_viridis(option="magma",discrete=T,name="Time")+ # Add both color and fill to get the desired graphical output
    scale_colour_viridis(option="magma",discrete=T,name="Time")
# If only fill -> legend is black..  plot_ordination does not have a fill option
rhizo_cold_fig$layers <- rhizo_cold_fig$layers[-1] 
rhizo_cold_fig


rhizo_beta_div <- ggarrange(rhizo_warm_fig,rhizo_cold_fig,ncol=2,common.legend = T,legend = "bottom")



## ROOTS #####


##### Separation by temperature #####

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

#### Beta Diversity ####

#### | Warm ####
ps_Warm_rclr <- microbiome::transform(ps_root_Warm, "rclr") # Robust Aitchison transformation

# Permanova
rclr_dist_matrix <- phyloseq::distance(ps_Warm_rclr, method = "euclidean") # Robust Aitchison distance
metadata <- as(sample_data(ps_Warm_rclr), "data.frame") # Metadata

Permanova.Warm.rclr<-adonis2(rclr_dist_matrix~Plant.type*Time,
                             data = metadata,
                             permutations = 999)

# PCA plots #
ord_Warm_rclr <- phyloseq::ordinate(ps_Warm_rclr, "RDA", distance = "euclidean")

root_warm_fig<- phyloseq::plot_ordination(ps_Warm_rclr, ord_Warm_rclr, type="samples",shape='Plant.type',color='Time') + 
    theme_bw()+
    theme(legend.text.align = 0,
          text = element_text(size=17))+
    ggtitle(label = "Roots - Warm temperatures")+
    geom_point(size=4,aes(fill=sample_data(ps_Warm_rclr)$Time))+
    geom_point(size=4.1,color='black')+ # Workaround to get black outline on previous geom_points 
    scale_shape_manual(values=c(22,23),name="Plant type",
                       labels=c(expression(italic(Scirpus)),expression(italic(Triglochin))))+
    scale_fill_viridis(option="magma",discrete=T,name="Time")+ # Add both color and fill to get the desired graphical output
    scale_colour_viridis(option="magma",discrete=T,name="Time")
# If only fill -> legend is black..  plot_ordination does not have a fill option
root_warm_fig$layers <- root_warm_fig$layers[-1] 
root_warm_fig

#### | Cold ####
ps_Cold_rclr <- microbiome::transform(ps_root_Cold, "rclr") # Robust Aitchison transformation

# Permanova
rclr_dist_matrix <- phyloseq::distance(ps_Cold_rclr, method = "euclidean") # Robust Aitchison distance
metadata <- as(sample_data(ps_Cold_rclr), "data.frame") # Metadata

Permanova.Cold.rclr<-adonis2(rclr_dist_matrix~Plant.type*Time,
                             data = metadata,
                             permutations = 999)

# PCA plots #
ord_Cold_rclr <- phyloseq::ordinate(ps_Cold_rclr, "RDA", distance = "euclidean")

root_cold_fig<- phyloseq::plot_ordination(ps_Cold_rclr, ord_Cold_rclr, type="samples",shape='Plant.type',color='Time')+
    theme_bw()+
    theme(legend.text.align = 0,
          text = element_text(size=17))+
    ggtitle(label = "Roots - Cold temperatures")+
    geom_point(size=4,aes(fill=sample_data(ps_Cold_rclr)$Time))+
    geom_point(size=4.1,color='black')+ # Workaround to get black outline on previous geom_points 
    scale_shape_manual(values=c(22,23),name="Plant type",
                       labels=c(expression(italic(Scirpus)),expression(italic(Triglochin))))+
    scale_fill_viridis(option="magma",discrete=T,name="Time")+ # Add both color and fill to get the desired graphical output
    scale_colour_viridis(option="magma",discrete=T,name="Time")
# If only fill -> legend is black..  plot_ordination does not have a fill option
root_cold_fig$layers <- root_cold_fig$layers[-1] 
root_cold_fig

root_beta_div <- ggarrange(root_warm_fig,root_cold_fig,ncol=2,common.legend = T,legend = "bottom")


## Figure 1#####

Beta_div <- ggarrange(Sed_beta_div,rhizo_beta_div,root_beta_div,nrow=3)


ggsave(filename = here("Results_W&C/Figures/", "Figure_1.pdf"),
       plot=Beta_div, height = 15, width = 12, dpi = 300)

ggsave(filename = here("Results_W&C/Figures/", "Figure_1.png"),
       plot=Beta_div, height = 15, width = 12, dpi = 300)

