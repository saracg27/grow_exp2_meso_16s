### Packages#####
library(here)
source(here("source", "libraries.R"))
library(phyloseq)
library(viridis)
library(metacoder)
library(dplyr)
library(ggpp)

#### Data ####

load(here("Rdata","ps_water_obj.RData"))


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

ggsave(here("Results","Figures","Water_SType_AlphaDiv.pdf"),device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results","Figures","Water_SType_AlphaDiv.png"),device='png',height = 7.5, width = 10.5)


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

ggsave(here("Results","Figures","Water_Temp_AlphaDiv.pdf"),device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results","Figures","Water_Temp_AlphaDiv.png"),device='png',height = 7.5, width = 10.5)

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


ggsave(here("Results","Figures","Water_TimeSType_AlphaDiv.pdf"),device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results","Figures","Water_TimeSType_AlphaDiv.png"),device='png',height = 7.5, width = 10.5)



###Rarefied abundance ####

# Rarefaction curves
rar_obj<- as.data.frame(t(otu_table(ps)))
rarecurve(rar_obj, step=100, cex=0.5)


### Number of reads per sample
min(sample_sums(ps)) #17392
max(sample_sums(ps)) #54723

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

# Not same number of samples.. 
sdat <- data.frame(sample_data(ps))
sum_sdat <- summary(sdat)
data.frame(sum_sdat[1:7,6])


### Rarefaction 
set.seed(12)
ps.rar <- rarefy_even_depth(ps) #
min(sample_sums(ps.rar))
max(sample_sums(ps.rar))
# Rarefied to 17392 reads, no loss occured


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

#Similar to non-rarefied


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

#Similar to non-rarefied

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


# Significant difference between D10 and D35 in Observed - Scirpus 
# not observed in rarefied dataset.
#Otherwise results are similar with some changes in pvalue significance


