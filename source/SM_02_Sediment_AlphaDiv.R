### Packages#####
library(here)
source(here("source", "libraries.R"))
library(phyloseq)
library(viridis)
library(metacoder)
library(dplyr)
library(ggpp)

#### Data ####

load(here("Rdata","ps_obj.RData"))

ps_sed <-  subset_samples(ps, Sample.type =="Sediment")
ps_sed <- prune_taxa(taxa_sums(ps_sed)>0, ps_sed)
ps_sed
# 120 samples 18851 taxa

# Simplify sample names 
sample_names(ps_sed) <- sub("Meso.*","",sample_names(ps_sed))


###______#####
### ALPHA DIV ####
### *** Raw abundance ####
alpha.div<- plot_richness(ps_sed, x="Time", color="Time", measures=c("Observed","Chao1","Simpson", "Shannon"))

alpha.div.plot <- alpha.div+
                  theme_bw()+
                  theme(title = element_text(size=20),
                        axis.text.x =element_blank(),
                        axis.text.y =element_text(size=12),
                        axis.title.x = element_blank(), 
                        axis.title.y = element_text(size=15),
                        axis.ticks.x = element_blank(),
                        legend.text=element_text(size=20))+
                  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),colour="grey50")+
                  #adds mean and standard deviation 
                  stat_pwc(method="wilcox_test",hide.ns=T,
                           p.adjust.method="BH",label="p.adj.signif",
                           tip.length = 0, step.increase = 0.03,vjust=1)+
                  # Computes pairwise wilcoxon rank-sum test in each facet
                  scale_colour_viridis(option="magma",discrete=T,name="Time")

### Checking stats ### 
####  Observed  
kruskal.test(value ~ Time, data = Rich.div) 
# p-value = 1.736e-10

pairwise.wilcox.test(Rich.div$value, Rich.div$Time,
                     p.adjust.method = "BH")
#        D4      D13     D20     D40    
#  D13 0.13    -       -       -      
#  D20 0.47    0.43    -       -      
#  D40 0.16    0.93    0.49    -      
#  D66 1.7e-07 1.7e-07 1.7e-07 1.7e-07

####  CHAO1 
kruskal.test(value ~ Time, data = Chao.div) 
# p-value = 1.864e-10

pairwise.wilcox.test(Chao.div$value, Chao.div$Time,
                     p.adjust.method = "BH")
#       D4      D13     D20     D40    
#   D13 0.14    -       -       -      
#   D20 0.48    0.44    -       -      
#   D40 0.17    0.89    0.48    -      
#   D66 7.5e-10 1.9e-09 1.4e-09 2.3e-09

####  SHANNON 
kruskal.test(value ~ Time, data = Shannon.div) 
# p-value = 1.527e-08

pairwise.wilcox.test(Shannon.div$value, Shannon.div$Time,
                     p.adjust.method = "BH")

#      D4      D13     D20     D40    
#  D13 0.14    -       -       -      
#  D20 0.47    0.42    -       -      
#  D40 0.47    0.23    0.83    -      
#  D66 1.9e-07 5.7e-07 2.5e-07 1.9e-07


####  SIMPSON 
kruskal.test(value ~ Time, data = Simpson.div) 
# p-value = 8.003e-06

pairwise.wilcox.test(Simpson.div$value, Simpson.div$Time,
                     p.adjust.method = "BH")
#       D4      D13     D20     D40    
#  D13 0.25516 -       -       -      
#  D20 0.65333 0.53300 -       -      
#  D40 0.95928 0.18609 0.65333 -      
#  D66 6.2e-05 0.00018 9.1e-05 6.2e-05



### *** Rarefied Abundance ####

# Rarefaction curves
rar_obj<- as.data.frame(t(otu_table(ps_sed)))
rarecurve(rar_obj, step=100, cex=0.5)

min(sample_sums(ps_sed))
max(sample_sums(ps_sed))

readsumsdf = data.frame(nreads = sample_sums(ps_sed), 
                        type = "Samples")

sdat <- as.data.frame(sample_data(ps_sed))

reads <- cbind(readsumsdf,sdat)

ggplot(reads, aes(x = rownames(reads), y = nreads,fill=Time))+
  theme_bw()+
  geom_bar(stat = "identity",color='black',size=0.2)+
  facet_grid(~Time,drop=T,scales="free")+
  theme(axis.text.x =element_text(angle=45,size=5))+
  scale_y_continuous(expand = expansion(c(0, 0.05)))+
  scale_fill_viridis(option="magma",discrete=T,name="Time")

ggsave(here("Results","Figures","Sed_Read#.pdf"),device='pdf',height = 7.5, width = 10.5)


### Rarefaction 
ps_sed_sub <- prune_samples(sample_sums(ps_sed)>10000,ps_sed)
ps_sed_sub <- prune_taxa(taxa_sums(ps_sed_sub)>0, ps_sed_sub)
# Loss of 5 D4 samples (8,15,27,30,31) and 13 ASVs only present in those samples

set.seed(12)
ps.rar <- rarefy_even_depth(ps_sed_sub)
min(sample_sums(ps.rar))
max(sample_sums(ps.rar))
# Rarefied to 20816 reads,  36 ASVs removed 



alpha.div.rar<- plot_richness(ps.rar, x="Time", color="Time", measures=c("Observed","Chao1","Simpson", "Shannon"))

alpha.div.rar.plot <- alpha.div.rar+
                      theme_bw()+
                      theme(title = element_text(size=20),
                            axis.text.x =element_blank(),
                            axis.text.y =element_text(size=12),
                            axis.title.x = element_blank(), 
                            axis.title.y = element_text(size=15),
                            axis.ticks.x = element_blank(),
                            legend.text=element_text(size=20))+
                      stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),colour="grey50")+
                      stat_pwc(method="wilcox_test",hide.ns=T,
                               p.adjust.method="BH",label="p.adj.signif",
                               tip.length = 0, step.increase = 0.03,vjust=1)+
                      scale_colour_viridis(option="magma",discrete=T,name="Time")

### *** Figure + save ####
Alpha_div_sed <- ggarrange(alpha.div.plot,alpha.div.rar.plot,
                           nrow=2,
                           labels="AUTO",
                           common.legend = T,
                           legend="bottom")

ggsave(here("Results","Figures","Sed_AlphaDiv.pdf"),device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results","Figures","Sed_AlphaDiv.png"),device='png',height = 7.5, width = 10.5)

####  Final figure  ######
alpha.div<- plot_richness(ps_sed, x="Time", color="Time", measures=c("Observed", "Simpson", "Shannon"))
# removed Chao1 as the results were very similar with observed richness 
alpha.div.plot <- alpha.div+
  theme_bw()+
  theme(title = element_text(size=20),
        axis.text.x =element_blank(),
        axis.text.y =element_text(size=12),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=15),
        axis.ticks.x = element_blank(),
        legend.text=element_text(size=20),
        strip.text = element_text(size=15),
        strip.background = element_rect(fill='white'))+
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),colour="grey50")+
  stat_pwc(method="wilcox_test",hide.ns=T,
           p.adjust.method="BH",label="p.adj.signif",
           tip.length = 0, step.increase = 0.03,vjust=1)+
  scale_colour_viridis(option="magma",discrete=T,name="Time")

# Mean and sd 
mean_alpha <- alpha.div$data %>%
  group_by(Time,variable) %>%
  summarise(mean = mean(value), sd = sd(value))

# Figure save
ggsave(here("Results","Figures","SupFig_AlphaDiv_Sed.pdf"),device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results","Figures","SupFig_AlphaDiv_Sed.png"),device='png',height = 7.5, width = 10.5)
