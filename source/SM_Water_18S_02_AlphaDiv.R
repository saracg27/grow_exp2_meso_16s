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

load(here("Rdata","ps_18S_water_obj.RData"))


###______#####
### ALPHA DIV ####
###Raw abundance ####

#####Sample type ####
alpha.div<- plot_richness(ps, x="Sample_type", measures=c("Observed","Simpson", "Shannon"))

# Simplified object for plotting purposes 
plot <- select(alpha.div$data,c("samples","Time","Sample_type","Temperature","variable","value"))

# Alpha div plot with Sample type as a grouping factor
ggplot(plot,aes(x=variable,y=value,shape=Sample_type))+
  geom_boxplot(outlier.shape = NA)+
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

ggsave(here("Results","Figures","18S_Water_SType_AlphaDiv.pdf"),device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results","Figures","18S_Water_SType_AlphaDiv.png"),device='png',height = 7.5, width = 10.5)


#####Temperature ####
alpha.div<- plot_richness(ps, x="Temperature", measures=c("Observed","Simpson", "Shannon"))

# Simplified object for plotting purposes 
 
plot <- select(alpha.div$data,c("samples","Time","Sample_type","Temperature","variable","value"))


# Alpha div plot with Temperature as a grouping factor

ggplot(plot,aes(x=variable,y=value,color=Temperature))+
  geom_boxplot(outlier.shape = NA)+
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

ggsave(here("Results","Figures","18S_Water_Temp_AlphaDiv.pdf"),device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results","Figures","18S_Water_Temp_AlphaDiv.png"),device='png',height = 7.5, width = 10.5)

#####Time x Sample type #####
alpha.div <- plot_richness(ps, x="Time", measures=c("Observed","Simpson", "Shannon"))

# Simplified object for plotting purposes 
plot <- select(alpha.div$data,c("samples","Time","Sample_type","Temperature","variable","value"))

# Alpha div plot with Sample type and Time as a grouping factor
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


ggsave(here("Results","Figures","18S_Water_TimeSType_AlphaDiv.pdf"),device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results","Figures","18S_Water_TimeSType_AlphaDiv.png"),device='png',height = 7.5, width = 10.5)



###Rarefied abundance ####

# Rarefaction curves
rar_obj<- as.data.frame(t(otu_table(ps)))
rarecurve(rar_obj, step=100, cex=0.5)

### Number of reads per sample
min(sample_sums(ps)) #17392
max(sample_sums(ps)) #54723


# Number of reads per sample
readsumsdf = data.frame(nreads = sample_sums(ps), 
                        type = "Samples")

sdat <- data.frame(sample_data(ps))

# Merge sample data and number of reads per sample
reads <- cbind(readsumsdf,sdat)

# Plot number of reads per sample with time as a grouping factor
ggplot(reads, aes(x = rownames(reads), y = nreads,fill=Time))+
  theme_bw()+
  geom_bar(stat = "identity",color='black',size=0.2)+
  facet_grid(~Time,drop=T,scales="free")+
  theme(axis.text.x =element_text(angle=45,size=5))+
  scale_y_continuous(expand = expansion(c(0, 0.05)))+
  scale_fill_viridis(option="magma",discrete=T,name="Time")


### Rarefaction 
set.seed(12)
ps.rar <- rarefy_even_depth(ps) #
min(sample_sums(ps.rar))
max(sample_sums(ps.rar))
# Rarefied to 17392 reads, no ASV loss occurred


#####Sample type ####
alpha.div<- plot_richness(ps.rar, x="Sample_type", measures=c("Observed","Simpson", "Shannon"))

# Simplified object for plotting purposes
plot <- select(alpha.div$data,c("samples","Time","Sample_type","Temperature","variable","value"))

# Alpha div plot with Sample type as a grouping factor
ggplot(plot,aes(x=variable,y=value,shape=Sample_type))+
  geom_boxplot(outlier.shape = NA)+
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
alpha.div<- plot_richness(ps.rar, x="Temperature", measures=c("Observed","Simpson", "Shannon"))

# Simplified object for plotting purposes
plot <- select(alpha.div$data,c("samples","Time","Sample_type","Temperature","variable","value"))


# Alpha div plot with Sample type as a grouping factor
ggplot(plot,aes(x=variable,y=value,color=Temperature))+
  geom_boxplot(outlier.shape = NA)+
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

#####Time x Sample type #####
alpha.div <- plot_richness(ps.rar, x="Time", measures=c("Observed","Simpson", "Shannon"))

# Simplified object for plotting purposes
plot <- select(alpha.div$data,c("samples","Time","Sample_type","Temperature","variable","value"))

# Alpha div plot with Sample type and Time as a grouping factor
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


# Change in significance level but not new significant difference

