### Packages#####
# Load required packages
library(here)
source(here("source", "libraries.R"))
library(phyloseq)
library(viridis)
library(metacoder)
library(dplyr)
library(ggrepel)
library(ggpp)


#### ASV Data ####
# load phyloseq object
load(here("Rdata","ps_obj.RData"))

# Subset sediment samples and remove ASVs with null abundance
ps_sed <-  subset_samples(ps, Sample.type =="Sediment")
ps_sed <- prune_taxa(taxa_sums(ps_sed)>0, ps_sed)
# 18851 taxa in 120 samples 


# Keeping D62 as sediment samples for chemistry analysis 
# were taken at the end of the experiment
#
# Sediment samples from day 62 
ps_sed_62 <- subset_samples(ps_sed, Time =="D62")
ps_sed_62 <- prune_taxa(taxa_sums(ps_sed_62)>0, ps_sed_62)
# 13322 taxa in 24 samples 

# Order agglomeration
ps_sed_62order <- tax_glom(ps_sed_62,taxrank = "Order",NArm = F)
# 304 orders  in 24 samples

#### RCLR transformation #####
ps.order.rclr <- microbiome::transform(ps_sed_62order, "rclr") 
# Retrieve abundance table
ASV_order_rclr <- as.data.frame(otu_table(ps.order.rclr))


#### ~~~~~~~~######
#### >>> Sediment <<< #####
#### Sediment chemistry data ####
sed_data <- read.table(file = "Data/Env_measures/env_final_sediment_exp2_20220225.csv", 
                         dec = ".", sep = ",")
colnames(sed_data) <- sed_data[1,] # Rename colnames 
sed_data <- sed_data[-1,-2] # Removing first row (colnames)  and Temp column


####  Merge datasets ####

# Retrieve sample data from phyloseq object
sdata <- data.frame(sample_data(ps_sed_62order))
colnames(sdata)
# Keep only relevant taxa for right_join
sdata <- select(sdata,c("CFL_ID_Sediments","Mesocosm..","Plant.type","Temperature","Greenhouse"))

# Merge the datasets
new_sdata <- sdata %>% right_join(sed_data, by=c("Mesocosm..","Greenhouse","Plant.type"))
colnames(new_sdata)

# Keep relevant columns
sed_rda_data <- select(new_sdata,c("CFL_ID_Sediments",
                                     "Plant.type",
                                     "Mesocosm..",
                                     "Temperature.x",
                                     "pH",
                                     "Electrical Conductivity_ds/m",
                                     "SAR",                         
                                     "% Saturation",                 
                                     "Calcium_meq/L",                
                                     "Calcium_mg/kg" ,              
                                     "Magnesium_meq/L",              
                                     "Magnesium_mg/kg",             
                                     "Sodium_meq/L",                
                                     "Sodium_mg/kg",                 
                                     "Potassium_meq/L",              
                                     "Potassium_mg/kg",             
                                     "Chloride_meq/L",              
                                     "Chloride_mg/kg",               
                                     "Sulfate(SO4)_meq/L",          
                                     "Sulfate(SO4)_mg/kg",          
                                     "Sulfate-S_meq/L",             
                                     "Sulfate-S_mg/kg")) 
                                     #"TGR_T/ac"))  # below detection limit        


# Setting ID as rownames and removing column
rownames(sed_rda_data) <- sed_rda_data[,1]
sed_rda_data <- sed_rda_data[,-1]

# Formating the data 9numbers to numeric , characters to factor)
sed_rda_data[,4:21] <- sapply(sed_rda_data[,4:21], as.numeric)
sed_rda_data[sapply(sed_rda_data, is.character)] <- lapply(sed_rda_data[sapply(sed_rda_data, is.character)], as.factor) # did it work? Check with str(meta)

# Standardizing the data
sed_rda_data[,4:21] <- decostand(sed_rda_data[,4:21] , method = "standardize")

# Checking standaridization
# round(apply(sed_rda_data[,4:21], 2, mean), 1) # Variables are now centered around a mean of 0
# apply(sed_rda_data[,4:21], 2, sd) # and scaled to have a standard deviation of 1

heatmap(abs(cor(sed_rda_data[,4:21])), 
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("topright", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))
# lots of correlation


####  RDA #####
ASV_object <- ASV_order_rclr

RDA <- vegan::rda(t(ASV_object) ~ ., data = sed_rda_data[,4:21])

RsquareAdj(RDA) 
#adj.r.squared = 0.1613894 for order agglomerated ASVs - Robust Ait

# Ordistep to keep only relevant variables
fwd.sel <- ordiR2step(rda(t(ASV_object) ~ 1, data = sed_rda_data[,4:21]), 
                      scope = formula(RDA), # Complete model 
                      direction = "forward",
                      R2scope = TRUE, # limited by the R2 of the complete modele 
                      #pstep = 1000,
                      trace = TRUE) # TRUE shows the selection process
fwd.sel$call
# robust Ait
# pH  for order agglomerated ASVs


RDA_signif <- vegan::rda(t(ASV_object) ~ pH , data = sed_rda_data[,4:21])
RsquareAdj(RDA_signif) 
#0.03302412 for Rclr orde glom
anova.cca(RDA_signif, permutations = 999) 

sqrt(vif.cca(RDA_signif)) # sqrt(vif())>2 is considered highly collinear

###### Scaling 1 #####
res.rda<-summary(RDA_signif,scaling = 1) # SCALING 1 
head(res.rda)
coord.asv.rda<-as.data.frame(res.rda$species) # ASV coordinates
coord.asv.rda2 <- cbind(paste0("ASV",seq(1:nrow(coord.asv.rda))),coord.asv.rda)
coord.asv.rda3 <- coord.asv.rda2
row.names(coord.asv.rda3) <- coord.asv.rda2[,1]
coord.asv.rda3 <- coord.asv.rda3[,-1]

coord.sites.rda<-as.data.frame(res.rda$sites)# Site coordinates
rownames(coord.sites.rda) <- substr(rownames(coord.sites.rda), 0, 3) #Rename sample names


coord.env <- as.data.frame(res.rda$biplot) # Environmental data coordinates
coord.env.mul <- coord.env*5 # Multiplication factor for plotting purposes


## Annotations for R2 and P in graph
df.annotations <- data.frame(
  label = paste(paste0("~italic(R)^{2} == ", round(RsquareAdj(RDA_signif)$adj.r.squared,3)),"~~",
                paste0("~italic(p) == ",round(anova.cca(RDA_signif, permutations = 1000)$'Pr(>F)'[1],3))))


### Labelling "extreme" taxa
 Order_annotation <-  subset(coord.asv.rda,abs(coord.asv.rda$RDA1)>1 | abs(coord.asv.rda$PC1)>1)
 tax <- as.data.frame(tax_table(ps.order.rclr))
 Interest <- subset(tax,rownames(tax)%in%rownames(Order_annotation))
 Interest <- cbind(Interest,Order_annotation)  


Rda.plot.rclr.S1<-ggplot(data=coord.sites.rda, aes(x=RDA1, y=PC1))+
    theme_bw()+
    geom_point(aes(shape=new_sdata$Plant.type,color=new_sdata$Temperature.x),size=2.8)+
    geom_point(data=coord.asv.rda,aes(x=RDA1, y=PC1),colour="purple4",shape=4, alpha=0.7)+
    geom_text_repel(data=Interest,label=Interest$Order,colour="purple4",
                    size=3.5,
                    min.segment.length = 0,
                    arrow=arrow(angle = 30, length = unit(0.01, "inches"),ends = "first", type = "open"),
                    fontface = "italic")+ 
    geom_hline(yintercept=0, linetype="dotted") +  
    geom_vline(xintercept=0, linetype="dotted") +
    geom_segment(data= coord.env.mul, aes(x=0, xend=RDA1, y=0, yend=PC1), 
                 color="tomato4", arrow=arrow(length=unit(0.01,"npc")))+
    geom_text_repel(data=coord.env.mul, aes(label=rownames(coord.env.mul)),
                    color="tomato4",force=2, size=4,fontface="bold")+
    labs(x = paste0("RDA1 (",round(res.rda$cont$importance[2,1]*100,2),"%)"),#Variance explained by each axis
         y = paste0("PC1 (",round(res.rda$cont$importance[2,2]*100,2),"%)"),#Variance explained by each axis
         title="RDA constrained by sediment chemistry",
         subtitle="Robust Aitchison transformed D62 abundance matrix as response variables - Scaling 1")+
    geom_label_npc(data= df.annotations , 
                   aes(npcx = "left", npcy = "bottom", label = label),
                   parse=T,size=6)+
    scale_shape_manual(name="Sample type",values=c(21,22,23),lables)+
    scale_color_manual(name="Temperature",values=c("blue","red"),labels=c("Cold","Warm"))

Rda.plot.rclr.S1


ggsave(here("Results","Figures","RDA_S1_D62_EnvSediment_RCLR_OrderGlom.pdf"),device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results","Figures","RDA_S1_D62_EnvSediment_RCLR_OrderGlom.png"),device='png',height = 7.5, width = 10.5)


###### Scaling 2 #####
res.rda<-summary(RDA_signif,scaling = 2) # SCALING 2 
head(res.rda)
coord.asv.rda<-as.data.frame(res.rda$species) # ASV coordinates
coord.asv.rda2 <- cbind(paste0("ASV",seq(1:nrow(coord.asv.rda))),coord.asv.rda)
coord.asv.rda3 <- coord.asv.rda2
row.names(coord.asv.rda3) <- coord.asv.rda2[,1]
coord.asv.rda3 <- coord.asv.rda3[,-1]

coord.sites.rda<-as.data.frame(res.rda$sites)# Site coordinates
rownames(coord.sites.rda) <- substr(rownames(coord.sites.rda), 0, 3) #Rename sample names


coord.env <- as.data.frame(res.rda$biplot) # Environmental data coordinates


## Annotations for R2 and P in graph
df.annotations <- data.frame(
  label = paste(paste0("~italic(R)^{2} == ", round(RsquareAdj(RDA_signif)$adj.r.squared,3)),"~~",
                paste0("~italic(p) == ",round(anova.cca(RDA_signif, permutations = 1000)$'Pr(>F)'[1],3))))


### Labelling "extreme" taxa
Order_annotation <-  subset(coord.asv.rda,abs(coord.asv.rda$RDA1)>1 | abs(coord.asv.rda$PC1)>1)
tax <- as.data.frame(tax_table(ps.order.rclr))
Interest <- subset(tax,rownames(tax)%in%rownames(Order_annotation))
Interest <- cbind(Interest,Order_annotation)  



Rda.plot.rclr.S2<-ggplot(data=coord.sites.rda, aes(x=RDA1, y=PC1)) +
  theme_bw()+
  geom_point()+
  geom_text_repel(label=rownames(coord.sites.rda))+
  geom_point(data=coord.asv.rda,aes(x=RDA1, y=PC1),colour="red",shape=4, alpha=0.2)+
  geom_text_repel(data=Interest,label=Interest$Order,colour="red",size=2,fontface = "italic")+ 
  geom_hline(yintercept=0, linetype="dotted") +  
  geom_vline(xintercept=0, linetype="dotted") +
  geom_segment(data= coord.env, aes(x=0, xend=RDA1, y=0, yend=PC1), 
               color="blue", arrow=arrow(length=unit(0.01,"npc")))+
  geom_text_repel(data=coord.env, aes(label=rownames(coord.env)),
                  color="blue", size=3.5)+
  geom_label_npc(data= df.annotations , 
                 aes(npcx = "left", npcy = "bottom", label = label),
                 parse=T,size=6)+
  labs(x = paste0("RDA1 (",round(res.rda$cont$importance[2,1]*100,2),"%)"),#Variance explained by each axis
       y = paste0("PC1 (",round(res.rda$cont$importance[2,2]*100,2),"%)"),#Variance explained by each axis
       title="RDA constrained by sediment chemistry",
       subtitle="Robust Aitchison transformed abundance matrix as response variables - Scaling 2")

Rda.plot.rclr.S2


ggsave(here("Results","Figures","RDA_S2_D62_EnvSediment_RCLR_OrderGlom.pdf"),device='pdf')
ggsave(here("Results","Figures","RDA_S2_D62_EnvSediment_RCLR_OrderGlom.png"),device='png')




##### ~~~~~~~~######
##### #### Warning #####
#
## Constraining sediment communities with water chemistry data might not be relevent
#
##### >>> WATER <<< #####
#### Water chemistry data ####
#water_data <- read.table(file = "Data/Env_measures/env_final_water_exp2_20220225.csv", 
#                         dec = ".", sep = ",")
#colnames(water_data) <- water_data[1,] # Rename colnames (similar column name in original table)
#water_data <- water_data[-1,-c(2,11)] # removing first row (colnames) 
## and one of the Electrical conductivity columns and Temp
#colnames(water_data)
#colnames(water_data)[4] <- "Mesocosm.."# Change name to match the colname in the sample data below 
#
#
#####  Merge datasets ####
#
#sdata <- data.frame(sample_data(ps_sed_62)) # Extract dataset from phyloseq object
#
## check the colnames and keep the usefull ones
#colnames(sdata)
#sdata <- select(sdata,c("CFL_ID_Sediments","Mesocosm..","Plant.type","Temperature","Greenhouse"))
#
## Merge the sample data and water chemistry data using several matching columns 
#new_sdata <- sdata %>% right_join(water_data, by=c("Mesocosm..","Greenhouse","Plant.type"))
#colnames(new_sdata)
#
#### Subset resulting merged dataset to keep relevant columns
#water_rda_data <- select(new_sdata,c("CFL_ID_Sediments",
#                                     "Plant.type",
#                                     "Mesocosm..",
#                                     "Temperature.x",
#                                     "pH",
#                                     #"Temperature of observed pH", #Not relevant
#                                     "Electrical Conductivity",   
#                                     "Calcium",   
#                                     "Magnesium",                 
#                                     "Sodium",                     
#                                     "Potassium",                 
#                                     #"Iron",  #below detection limit                     
#                                     #"Manganese", #below detection limit              
#                                     "Chloride",                   
#                                     #"Nitrate - N",  # most below detection limit              
#                                     #"Nitrite - N",   #below detection limit             
#                                     #"Nitrate and Nitrite - N",   # most below detection limit 
#                                     "Sulfate (SO4)",              
#                                     #"Hydroxide",    #below detection limit             
#                                     #"Carbonate",    #below detection limit              
#                                     "Bicarbonate",               
#                                     #"P-Alkalinity",   #below detection limit          
#                                     "T-Alkalinity",            
#                                     "Total Dissolved Solids",    
#                                     "Hardness",                  
#                                     "Ionic Balance"))
#
#
#
###### Water data ####
#
## Setting ID as rownames and removing ID column
#rownames(water_rda_data) <- water_rda_data[,1]
#water_rda_data <- water_rda_data[,-1]
#
## Formatting the chemistry data
#water_rda_data[,4:16] <- sapply(water_rda_data[,4:16], as.numeric)
#water_rda_data[sapply(water_rda_data, is.character)] <- lapply(water_rda_data[sapply(water_rda_data, is.character)], as.factor) # did it work? Check with str(meta)
#
## Standardizing the chemistry data
#water_rda_data[,4:16] <- decostand(water_rda_data[,4:16] , method = "standardize")
#
## Checking standaridization
## round(apply(water_rda_data[,4:16], 2, mean), 1) # Variables are now centered around a mean of 0
## apply(water_rda_data[,4:16], 2, sd) # and scaled to have a standard deviation of 1
#
## Heatmap to chekc colinearity between variables
#heatmap(abs(cor(water_rda_data[,4:16])), 
#        # Compute pearson correlation (note they are absolute values)
#        col = rev(heat.colors(6)), 
#        Colv = NA, Rowv = NA)
#legend("topright", 
#       title = "Absolute Pearson R",
#       legend =  round(seq(0,1, length.out = 6),1),
#       y.intersp = 0.7, bty = "n",
#       fill = rev(heat.colors(6)))
## T alkanity and bicarbonate are correlated
## hardeness and dissolved solids as well
#
#
#
##### Water data RDA #####
#
## Select the abundance dataset to use
#ASV_object <- ASV_order_rclr
#
## Compute the rda using chemistry data as constraining variable 
#RDA <- vegan::rda(t(ASV_object) ~ ., data = water_rda_data[,4:16])
#
#
#RsquareAdj(RDA) # adj.r.squared = 0.04788229 for un-agglomerated ASVs
## adj.r.squared = 0.1173067 for genus agglomerated ASVs
## adj.r.squared = 0.1081093 for order agglomerated ASVs
#
## $adj.r.squared = 0.1200128 for order agglomerated ASVs on RCLR transformed data
#
## Forward selection of variables with the most explanating power 
#fwd.sel <- ordiR2step(rda(t(ASV_object) ~ 1, data = water_rda_data[,4:16]), # modèle le plus simple
#                      scope = formula(RDA), # modèle "complet"
#                      direction = "forward",
#                      R2scope = TRUE, # limité par le R2 du modèle "complet"
#                      #pstep = 1000,
#                      trace = TRUE) # mettre TRUE pour voir le processus du sélection!
#fwd.sel$call
### Aitchison
## Chloride for un-aglomerated ASVs
## Chloride + Potassium + Bicarbonate for genus agglomerated ASVs
## Chloride + Potassium for order agglomerated ASVs
#
### Robust Aitchison
## Chloride + Potassium +`Electrical Conductivity` for order agglomerated ASVs
#
## Checking the significance of variables
#anova.cca(RDA, permutations = 1000)
#anova.cca(RDA, permutations = 1000,by="term")
## Calcium, Sodium, ph?
## Calcium Sodium for genus agglomerated
## Calcium Sodium for order agglomerated
#
#
## Reducing the number of variables keeping the ones identified by ordiplot
#RDA_signif <- vegan::rda(t(ASV_object) ~ Chloride + Potassium +`Electrical Conductivity`, data = water_rda_data[,4:16])
#RsquareAdj(RDA_signif)
#anova.cca(RDA_signif, permutations = 1000,by="term")
#
## Variance inflation factors function to identify useles constraining variables
#sqrt(vif.cca(RDA_signif)) # sqrt(vif())>2 is considered highly collinear
#
#
###### Scaling 1 #####
#res.rda<-summary(RDA_signif,scaling = 1) # SCALING 1 
#head(res.rda) # Check
#
## ASV coordinates
#coord.asv.rda<-as.data.frame(res.rda$species) 
#
## Change rownames from ASV sequence to ASV# for eventual plotting purposes
#coord.asv.rda2 <- cbind(paste0("ASV",seq(1:nrow(coord.asv.rda))),coord.asv.rda) # Change ASV sequence by ASV#
#coord.asv.rda3 <- coord.asv.rda2
#row.names(coord.asv.rda3) <- coord.asv.rda2[,1]
#coord.asv.rda3 <- coord.asv.rda3[,-1]
#
## Site coordinates
#coord.sites.rda<-as.data.frame(res.rda$sites)
##Renaming sample names for plotting purposes, keeping only the 3 first characters
#rownames(coord.sites.rda) <- substr(rownames(coord.sites.rda), 0, 3) 
#
## Environmental data coordinates
#coord.env <- as.data.frame(res.rda$biplot)
#coord.env.mul <- coord.env*5 # Multiplication factor for plotting purposes
#
#
### Annotations for R2 and P in graph
#df.annotations <- data.frame(
#    label = paste(paste0("~italic(R)^{2} == ", round(RsquareAdj(RDA_signif)$adj.r.squared,3)),"~~",
#                  paste0("~italic(p) == ",round(anova.cca(RDA_signif, permutations = 1000)$'Pr(>F)'[1],3))))
#
#
#### Labelling "extreme" taxa for plotting purposes
#Order_annotation <-  subset(coord.asv.rda,abs(coord.asv.rda$RDA1)>1 | abs(coord.asv.rda$RDA2)>1)
#tax <- as.data.frame(tax_table(ps.order.rclr))
#Interest <- subset(tax,rownames(tax)%in%rownames(Order_annotation))
#Interest <- cbind(Interest,Order_annotation)  
#
#
## RDA plot # 
#Rda.plot.rclr.S1<-ggplot(data=coord.sites.rda, aes(x=RDA1, y=RDA2)) +
#    theme_bw()+
#    geom_point()+
#    geom_text_repel(label=rownames(coord.sites.rda))+
#    geom_point(data=coord.asv.rda,aes(x=RDA1, y=RDA2),colour="red",shape=4, alpha=0.2)+
#    geom_text_repel(data=Interest,label=Interest$Order,colour="red",size=2,fontface = "italic")+ 
#    geom_hline(yintercept=0, linetype="dotted") +  
#    geom_vline(xintercept=0, linetype="dotted") +
#    geom_segment(data= coord.env.mul, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
#                 color="blue", arrow=arrow(length=unit(0.01,"npc")))+
#    geom_text_repel(data=coord.env.mul, aes(label=rownames(coord.env.mul)),
#                    color="blue", size=3.5)+
#    labs(x = paste0("RDA1 (",round(res.rda$cont$importance[2,1]*100,2),"%)"),#Variance explained by each axis
#         y = paste0("RDA2 (",round(res.rda$cont$importance[2,2]*100,2),"%)"),#Variance explained by each axis
#         title="RDA constrained by water chemistry",
#         subtitle="Robust Aitchison transformed abundance matrix as response variables - Scaling 1")+
#    geom_label_npc(data= df.annotations , 
#                   aes(npcx = "right", npcy = "bottom", label = label),
#                   parse=T,size=6)
#
#Rda.plot.rclr.S1
#
#
#ggsave(here("Results","Figures","RDA_S1_D62_EnvWater_RCLR_OrderGlom.pdf"),device='pdf')
#ggsave(here("Results","Figures","RDA_S1_D62_EnvWater_RCLR_OrderGlom.png"),device='png')
#
#
###### Scaling 2 #####
#res.rda<-summary(RDA_signif,scaling = 2) # SCALING 2 
#head(res.rda)
#res.rda$cont$importance
#coord.asv.rda<-as.data.frame(res.rda$species) # ASV coordinates
#coord.asv.rda2 <- cbind(paste0("ASV",seq(1:nrow(coord.asv.rda))),coord.asv.rda)
#coord.asv.rda3 <- coord.asv.rda2
#row.names(coord.asv.rda3) <- coord.asv.rda2[,1]
#coord.asv.rda3 <- coord.asv.rda3[,-1]
#
#coord.sites.rda<-as.data.frame(res.rda$sites)# Site coordinates
#rownames(coord.sites.rda) <- substr(rownames(coord.sites.rda), 0, 3) #Rename sample names
#
#
#coord.env <- as.data.frame(res.rda$biplot) # Environmental data coordinates
#
#
### Annotations for R2 and P in graph
#df.annotations <- data.frame(
#    label = paste(paste0("~italic(R)^{2} == ", round(RsquareAdj(RDA_signif)$adj.r.squared,3)),"~~",
#                  paste0("~italic(p) == ",round(anova.cca(RDA_signif, permutations = 1000)$'Pr(>F)'[1],3))))
#
#
#### Labelling "extreme" taxa
#Order_annotation <-  subset(coord.asv.rda,abs(coord.asv.rda$RDA1)>1 | abs(coord.asv.rda$RDA2)>1)
#tax <- as.data.frame(tax_table(ps.order.rclr))
#Interest <- subset(tax,rownames(tax)%in%rownames(Order_annotation))
#Interest <- cbind(Interest,Order_annotation)  
#
#
#
#Rda.plot.rclr.S2<-ggplot(data=coord.sites.rda, aes(x=RDA1, y=RDA2)) +
#    theme_bw()+
#    geom_point()+
#    geom_text_repel(label=rownames(coord.sites.rda))+
#    geom_point(data=coord.asv.rda,aes(x=RDA1, y=RDA2),colour="red",shape=4, alpha=0.2)+
#    geom_text_repel(data=Interest,label=Interest$Order,colour="red",size=2,fontface = "italic")+ 
#    geom_hline(yintercept=0, linetype="dotted") +  
#    geom_vline(xintercept=0, linetype="dotted") +
#    geom_segment(data= coord.env, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
#                 color="blue", arrow=arrow(length=unit(0.01,"npc")))+
#    geom_text_repel(data=coord.env, aes(label=rownames(coord.env)),
#                    color="blue", size=3.5)+
#    labs(x = paste0("RDA1 (",round(res.rda$cont$importance[2,1]*100,2),"%)"),#Variance explained by each axis
#         y = paste0("RDA2 (",round(res.rda$cont$importance[2,2]*100,2),"%)"),#Variance explained by each axis
#         title="RDA constrained by water chemistry",
#         subtitle="Robust Aitchison transformed abundance matrix as response variables - Scaling 2")+
#    geom_label_npc(data= df.annotations , 
#                   aes(npcx = "right", npcy = "bottom", label = label),
#                   parse=T,size=6)
#
#Rda.plot.rclr.S2
#
#
#ggsave(here("Results","Figures","RDA_S2_D62_EnvWater_RCLR_OrderGlom.pdf"),device='pdf')
#ggsave(here("Results","Figures","RDA_S2_D62_EnvWater_RCLR_OrderGlom.png"),device='png')
