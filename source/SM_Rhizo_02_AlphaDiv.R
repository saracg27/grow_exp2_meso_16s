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


