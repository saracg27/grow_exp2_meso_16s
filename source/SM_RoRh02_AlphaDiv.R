### Packages#####
library(here)
source(here("source", "libraries.R"))
library(phyloseq)
library(viridis)
library(dplyr)
library(ggpp)

#### Data ####

load(here("Rdata","ps_obj.RData"))

ps_plant <-  subset_samples(ps, Sample.type !="Sediment")
ps_plant <- prune_taxa(taxa_sums(ps_plant)>0, ps_plant)
ps_plant
# 162 samples 6319 taxa


###______#####
### ALPHA DIV ####
### Raw abundance ####

#### ... Plant type #####
alpha.div<- plot_richness(ps_plant, x="Plant.type", color="Plant.type", measures=c("Observed","Simpson", "Shannon"))

plot <- select(alpha.div$data,c("samples","Time","Sample.type","Plant.type","Temperature","variable","value"))
plot$Time <- as.factor(plot$Time)
plot$Sample.type <- factor(plot$Sample.type,levels=c("Rhizosphere","Roots"))


ggplot(plot,aes(x=variable,y=value,shape=Plant.type))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(x=variable,y=value,shape=Plant.type),position=position_dodge(width=0.75),size=2.5,fill='black')+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank())+
  stat_pwc(method="wilcox_test",
           hide.ns=F,
           label="p.adj.signif",
           tip.length = 0, 
           vjust=0.25)+
  scale_shape_manual(values=c(22,23),name="Plant")+
  facet_wrap(Sample.type~variable,scales="free")+
  ylab("Alpha diversity score")

ggsave(here("Results","Figures","RootRhizo_AlphaDiv.pdf"),device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results","Figures","RootRhizo_AlphaDiv.png"),device='png',height = 7.5, width = 10.5)

#### ... Time #####
alpha.dim.time<- plot_richness(ps_plant, x="Time", color="Plant.type", measures=c("Observed","Simpson", "Shannon"))

plot <- select(alpha.div$data,c("samples","Time","Sample.type","Plant.type","Temperature","variable","value"))

library(ggh4x)
ggplot(plot,aes(x=variable,y=value,shape=Plant.type,color=Time))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(x=variable,y=value,shape=Plant.type),position=position_dodge(width=0.75),size=2,fill='black')+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank())+
  stat_pwc(method="wilcox_test",
           p.adj.method="fdr",
           hide.ns=T,
           label="p.adj.signif",
           tip.length = 0, 
           vjust=0.7,
           size = 0.5,
           label.size=4.2)+
  scale_colour_viridis(option="magma",discrete=T,name="Time")+
  scale_shape_manual(values=c(22,23),name="Plant")+
  facet_nested(variable~Sample.type+Plant.type,scales="free",independent="all",drop=T)+
  ylab("Alpha diversity score")

ggsave(here("Results","Figures","RootRhizo_Time_AlphaDiv.pdf"),device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results","Figures","RootRhizo_Time_AlphaDiv.png"),device='png',height = 7.5, width = 10.5)



#### ... Temperature #####
alpha.dim.time<- plot_richness(ps_plant, x="Temperature", color="Plant.type", measures=c("Observed","Simpson", "Shannon"))

plot <- select(alpha.div$data,c("samples","Time","Sample.type","Plant.type","Temperature","variable","value"))

library(ggh4x)
ggplot(plot,aes(x=variable,y=value,shape=Plant.type,color=Temperature))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(x=variable,y=value,shape=Plant.type),position=position_dodge(width=0.75),size=2,fill='black')+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank())+
  stat_pwc(method="wilcox_test",
           p.adj.method="fdr",
           hide.ns=T,
           label="p.adj.signif",
           tip.length = 0, 
           vjust=0.7,
           size = 0.5,
           label.size=4.2)+
  scale_color_manual(values=c("#0000FE","#B02223"),name="Temperature")+
  scale_shape_manual(values=c(22,23),name="Plant")+
  facet_nested(variable~Sample.type+Plant.type,scales="free",independent="all",drop=T)+
  ylab("Alpha diversity score")

ggsave(here("Results","Figures","RootRhizo_Temp_AlphaDiv.pdf"),device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results","Figures","RootRhizo_Temp_AlphaDiv.png"),device='png',height = 7.5, width = 10.5)



#### Rarefied root ####

load(here("Rdata","ps_obj.RData"))

ps_root <-  subset_samples(ps, Sample.type =="Roots")
ps_root <- prune_taxa(taxa_sums(ps_root)>0, ps_root)

# Rarefaction curves
rar_obj<- as.data.frame(t(otu_table(ps_root)))
rarecurve(rar_obj, step=100, cex=0.5)


### Number of reads per sample
min(sample_sums(ps_root))
max(sample_sums(ps_root))

readsumsdf = data.frame(nreads = sample_sums(ps_root), 
                        type = "Samples")

sdat <- data.frame(sample_data(ps_root))

reads <- cbind(readsumsdf,sdat)

ggplot(reads, aes(x = rownames(reads), y = nreads,fill=Time))+
  theme_bw()+
  geom_bar(stat = "identity",color='black',size=0.2)+
  facet_grid(~Time,drop=T,scales="free")+
  theme(axis.text.x =element_text(angle=45,size=5))+
  scale_y_continuous(expand = expansion(c(0, 0.05)))+
  scale_fill_viridis(option="magma",discrete=T,name="Time")

### Rarefaction 
set.seed(12)
ps.rar <- rarefy_even_depth(ps_root) # Loss of 707 OTUs
min(sample_sums(ps.rar))
max(sample_sums(ps.rar))
# Rarefied to 1339 reads, loss of 707 ASVs


alpha.div.rar<- plot_richness(ps.rar, x="Time", color="Time", measures=c("Observed","Simpson", "Shannon"))

plot <- select(alpha.div.rar$data,c("samples","Time","Sample.type","Plant.type","Temperature","variable","value"))
plot$Time <- as.factor(plot$Time)
plot$Sample.type <- factor(plot$Sample.type,levels=c("Rhizosphere","Roots"))

ggplot(plot,aes(x=variable,y=value,shape=Plant.type))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(x=variable,y=value,shape=Plant.type),position=position_dodge(width=0.75),size=2.5,fill='black')+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank())+
  stat_pwc(method="wilcox_test",
           hide.ns=F,
           label="p.adj.signif",
           tip.length = 0, 
           vjust=0.25)+
  scale_shape_manual(values=c(22,23),name="Plant")+
  facet_wrap(Sample.type~variable,scales="free")+
  ylab("Alpha diversity score")

