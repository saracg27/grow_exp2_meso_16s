### Packages#####
library(here)
source(here("source", "libraries.R"))
library(phyloseq)
library(viridis)
library(metacoder)
library(dplyr)
library(ggpp)
library(pairwiseAdonis)
library(gridExtra)



##16S Water  #####
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


#### Beta Diversity ####
#### #### | Warm ####
ps_Warm_rclr <- microbiome::transform(ps_water_16S_Warm, "rclr") # Robust Aitchison transformation

# Permanova 
rclr_dist_matrix <- phyloseq::distance(ps_Warm_rclr, method = "euclidean") # Robust Aitchison distance
metadata <- as(sample_data(ps_Warm_rclr), "data.frame") # Metadata

Permanova.Warm.rclr<-adonis2(rclr_dist_matrix~Sample_type*Time,
                             data = metadata,
                             permutations = 999)

# PCA plots #
ord_Warm_rclr <- phyloseq::ordinate(ps_Warm_rclr, "RDA", distance = "euclidean")

Water16S_warm_fig<- phyloseq::plot_ordination(ps_Warm_rclr, ord_Warm_rclr, type="samples",shape='Sample_type',color='Time') + 
    theme_bw()+
    theme(legend.position = "bottom",
          legend.text.align = 0,
          text = element_text(size=17))+
    ggtitle(label = "Water 16S - Warm temperatures")+
    geom_point(size=4,aes(fill=sample_data(ps_Warm_rclr)$Time))+
    geom_point(size=4.1,color='black')+ # Workaround to get black outline on previous geom_points 
    scale_shape_manual(values=c(21,22,23),name="Plant type",
                       labels=c("No plant",expression(italic(Scirpus)),expression(italic(Triglochin))))+
    scale_fill_viridis(option="magma",discrete=T,name="Time")+ # Add both color and fill to get the desired graphical output
    scale_colour_viridis(option="magma",discrete=T,name="Time")+
    guides(fill=guide_legend(nrow=1,byrow=TRUE))

# If only fill -> legend is black..  plot_ordination does not have a fill option
Water16S_warm_fig$layers <- Water16S_warm_fig$layers[-1] 
Water16S_warm_fig

g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}

mylegend<-g_legend(Water16S_warm_fig)

#### | Cold ####
ps_cold_rclr <- microbiome::transform(ps_water_16S_Cold, "rclr") # Robust Aitchison transformation

# Permanova #
rclr_dist_matrix <- phyloseq::distance(ps_cold_rclr, method = "euclidean") # Robust Aitchison distance
metadata <- as(sample_data(ps_cold_rclr), "data.frame") # Metadata

Permanova.cold.rclr<-adonis2(rclr_dist_matrix~Sample_type*Time,
                             data = metadata,
                             permutations = 999)

# PCA plots #
ord_cold_rclr <- phyloseq::ordinate(ps_cold_rclr, "RDA", distance = "euclidean")

Water16S_cold_fig<- phyloseq::plot_ordination(ps_cold_rclr, ord_cold_rclr, type="samples",shape='Sample_type',color='Time') + 
    theme_bw()+
    theme(legend.text.align = 0,
          text = element_text(size=17))+
    ggtitle(label = "Water 16S - Cold temperatures")+
    geom_point(size=4,aes(fill=sample_data(ps_cold_rclr)$Time))+
    geom_point(size=4.1,color='black')+ # Workaround to get black outline on previous geom_points 
    scale_shape_manual(values=c(21,22,23),name="Plant type",
                       labels=c("No plant",expression(italic(Scirpus)),expression(italic(Triglochin))))+
    scale_fill_viridis(option="magma",discrete=T,name="Time")+ # Add both color and fill to get the desired graphical output
    scale_colour_viridis(option="magma",discrete=T,name="Time")
# If only fill -> legend is black..  plot_ordination does not have a fill option
Water16S_cold_fig$layers <- Water16S_cold_fig$layers[-1] 
Water16S_cold_fig

Water16S_beta_div <- ggarrange(Water16S_warm_fig,Water16S_cold_fig,ncol=2,common.legend = T,legend = "none")

###________#####
### 18S Water #####
#### Data ####

load(here("Rdata","ps_18S_water_obj.RData"))

#### Separation by temperature #####

ps_water_18S_Warm <- subset_samples(ps, Temperature =="20C day - 10C night")
ps_water_18S_Warm <- prune_taxa(taxa_sums(ps_water_18S_Warm)>0, ps_water_18S_Warm)
ps_water_18S_Warm # 70 samples - 637 taxa
# check to see if you kept only water_18S samples in WARM condition
str(sample_data(ps_water_18S_Warm)) 

ps_water_18S_Cold <- subset_samples(ps, Temperature =="10C day - 5C night")
ps_water_18S_Cold <- prune_taxa(taxa_sums(ps_water_18S_Cold)>0, ps_water_18S_Cold)
ps_water_18S_Cold # 70 samples - 627 taxa
# check to see if you kept only water_18S samples in COLD condition
str(sample_data(ps_water_18S_Cold)) 

#### Beta Diversity ####
#### #### | Warm ####
ps_Warm_rclr <- microbiome::transform(ps_water_18S_Warm, "rclr") # Robust Aitchison transformation

# Permanova 
rclr_dist_matrix <- phyloseq::distance(ps_Warm_rclr, method = "euclidean") # Robust Aitchison distance
metadata <- as(sample_data(ps_Warm_rclr), "data.frame") # Metadata

Permanova.Warm.rclr<-adonis2(rclr_dist_matrix~Sample_type*Time,
                             data = metadata,
                             permutations = 999)

# PCA plots #
ord_Warm_rclr <- phyloseq::ordinate(ps_Warm_rclr, "RDA", distance = "euclidean")

Water18S_warm_fig<- phyloseq::plot_ordination(ps_Warm_rclr, ord_Warm_rclr, type="samples",shape='Sample_type',color='Time') + 
    theme_bw()+
    theme(legend.text.align = 0,
          text = element_text(size=17))+
    ggtitle(label = "Water 18S - Warm temperatures")+
    geom_point(size=4,aes(fill=sample_data(ps_Warm_rclr)$Time))+
    geom_point(size=4.1,color='black')+ # Workaround to get black outline on previous geom_points 
    scale_shape_manual(values=c(21,22,23),name="Plant type",
                       labels=c("No plant",expression(italic(Scirpus)),expression(italic(Triglochin))))+
    scale_fill_viridis(option="magma",discrete=T,name="Time")+ # Add both color and fill to get the desired graphical output
    scale_colour_viridis(option="magma",discrete=T,name="Time")
# If only fill -> legend is black..  plot_ordination does not have a fill option
Water18S_warm_fig$layers <- Water18S_warm_fig$layers[-1] 
Water18S_warm_fig

#### | Cold ####
ps_cold_rclr <- microbiome::transform(ps_water_18S_Cold, "rclr") # Robust Aitchison transformation

# Permanova #
rclr_dist_matrix <- phyloseq::distance(ps_cold_rclr, method = "euclidean") # Robust Aitchison distance
metadata <- as(sample_data(ps_cold_rclr), "data.frame") # Metadata

Permanova.cold.rclr<-adonis2(rclr_dist_matrix~Sample_type*Time,
                             data = metadata,
                             permutations = 999)

# PCA plots #
ord_cold_rclr <- phyloseq::ordinate(ps_cold_rclr, "RDA", distance = "euclidean")

Water18S_cold_fig<- phyloseq::plot_ordination(ps_cold_rclr, ord_cold_rclr, type="samples",shape='Sample_type',color='Time') + 
    theme_bw()+
    theme(legend.text.align = 0,
          text = element_text(size=17))+
    ggtitle(label = "Water 18S - Cold temperatures")+
    geom_point(size=4,aes(fill=sample_data(ps_cold_rclr)$Time))+
    geom_point(size=4.1,color='black')+ # Workaround to get black outline on previous geom_points 
    scale_shape_manual(values=c(21,22,23),name="Plant type",
                       labels=c("No plant",expression(italic(Scirpus)),expression(italic(Triglochin))))+
    scale_fill_viridis(option="magma",discrete=T,name="Time")+ # Add both color and fill to get the desired graphical output
    scale_colour_viridis(option="magma",discrete=T,name="Time")
# If only fill -> legend is black..  plot_ordination does not have a fill option
Water18S_cold_fig$layers <- Water18S_cold_fig$layers[-1] 
Water18S_cold_fig

Water18S_beta_div <- ggarrange(Water18S_warm_fig,Water18S_cold_fig,ncol=2,common.legend = T,legend = "none")


#### FIGURE 2#####
Beta_div <- grid.arrange(arrangeGrob(Water16S_beta_div,
                               Water18S_beta_div,
                               nrow=2),
                   mylegend, nrow=2,heights=c(10, 1))


ggsave(filename = here("Results_W&C/Figures/", "Figure_2.pdf"),
       plot=Beta_div, height = 10, width = 12, dpi = 300)

ggsave(filename = here("Results_W&C/Figures/", "Figure_2.png"),
       plot=Beta_div, height = 10, width = 12, dpi = 300)



