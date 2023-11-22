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

# Keep only sediment samples and remove ASVs with 0 reads
ps_sed <-  subset_samples(ps, Sample.type =="Sediment")
ps_sed <- prune_taxa(taxa_sums(ps_sed)>0, ps_sed)
ps_sed
# 120 samples 18851 taxa

# Simplify sample names 
sample_names(ps_sed) <- sub("Meso.*","",sample_names(ps_sed))


###______#####
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
  facet_grid(~Time,drop=T,scales="free")+
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




