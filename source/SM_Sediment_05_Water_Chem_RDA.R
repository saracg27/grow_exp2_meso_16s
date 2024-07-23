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
ps_sed_62 <-  subset_samples(ps, Sample.type =="Sediment" & Time =="D62")
ps_sed_62 <- prune_taxa(taxa_sums(ps_sed_62)>0, ps_sed_62)
# 13322 taxa in 24 samples 

#### Separation by temperature #####

ps_warm <- subset_samples(ps_sed_62, Temperature =="20°C day 10°C night")
ps_warm <- prune_taxa(taxa_sums(ps_warm)>0, ps_warm)
str(sample_data(ps_warm)) 

ps_cold<- subset_samples(ps_sed_62, Temperature =="10°C day 5°C night")
ps_cold <- prune_taxa(taxa_sums(ps_cold)>0, ps_cold)
str(sample_data(ps_cold)) 


#### Order agglomeration ####
ps_warm_order <- tax_glom(ps_warm,taxrank = "Order",NArm = F)
# 282 orders  in 12 samples
ps_cold_order <- tax_glom(ps_cold,taxrank = "Order",NArm = F)
# 270 orders  in 12 samples


#### RCLR transformation #####
ps_warm_rclr <- microbiome::transform(ps_warm_order, "rclr") 
# Retrieve abundance table
ASV_warm_rclr <- as.data.frame(otu_table(ps_warm_rclr))

ps_cold_rclr <- microbiome::transform(ps_cold_order, "rclr") 
# Retrieve abundance table
ASV_cold_rclr <- as.data.frame(otu_table(ps_cold_order))


#### ~~~~~~~~######

#### Sediment chemistry data ####
sed_data <- read.table(file = "Data/Env_measures/env_final_sediment_exp2_20220225.csv", 
                         dec = ".", sep = ",")
colnames(sed_data) <- sed_data[1,] # Rename colnames 
sed_data <- sed_data[-1,-2] # Removing first row (colnames)  and Temp column


#### WARM #####
#### | Merge datasets ####

# Retrieve sample data from phyloseq object
sdata <- data.frame(sample_data(ps_warm_order))
colnames(sdata)
# Keep only relevant taxa for right_join
sdata <- select(sdata,c("CFL_ID_Sediments","Mesocosm..","Plant.type","Temperature","Greenhouse"))

# Merge the datasets
new_sdata <- sdata %>% inner_join(sed_data, by=c("Mesocosm..","Greenhouse","Plant.type"))
colnames(new_sdata)

# Keep relevant columns
sed_rda_data <- select(new_sdata,c("CFL_ID_Sediments",
                                     "Plant.type",
                                     "Mesocosm..",
                                     "Temperature.x",
                                     "pH",
                                     #"Electrical Conductivity_ds/m", #correlated to other variables
                                     #"SAR",# Correlated to Chloride                     
                                     "% Saturation",                 
                                     #"Calcium_meq/L",                
                                     "Calcium_mg/kg" ,              
                                     #"Magnesium_meq/L",              
                                     "Magnesium_mg/kg",             
                                     #"Sodium_meq/L",                
                                     "Sodium_mg/kg",                 
                                     #"Potassium_meq/L",              
                                     "Potassium_mg/kg",             
                                     #"Chloride_meq/L",              
                                     "Chloride_mg/kg",               
                                     #"Sulfate(SO4)_meq/L",          
                                     "Sulfate(SO4)_mg/kg"))      
                                     #"Sulfate-S_meq/L",   # Correlated to Sulfate SO4         
                                     #"Sulfate-S_mg/kg")) 
                                     #"TGR_T/ac"))  # below detection limit        


# Setting ID as rownames and removing column
rownames(sed_rda_data) <- sed_rda_data[,1]
sed_rda_data <- sed_rda_data[,-1]

# Formating the data 9numbers to numeric , characters to factor)
sed_rda_data[,4:11] <- sapply(sed_rda_data[,4:11], as.numeric)
sed_rda_data[sapply(sed_rda_data, is.character)] <- lapply(sed_rda_data[sapply(sed_rda_data, is.character)], as.factor) # did it work? Check with str(meta)

# Standardizing the data
sed_rda_data[,4:11] <- decostand(sed_rda_data[,4:11] , method = "standardize")

# Checking standaridization
# round(apply(sed_rda_data[,4:11], 2, mean), 1) # Variables are now centered around a mean of 0
# apply(sed_rda_data[,4:11], 2, sd) # and scaled to have a standard deviation of 1

heatmap(abs(cor(sed_rda_data[,4:11])), 
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("topright", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))
# Still some correlations (chloride and sodium ; sulfate and magnesium)


#### | RDA #####
ASV_object <- ASV_warm_rclr

RDA <- vegan::rda(t(ASV_object) ~ ., data = sed_rda_data[,4:11])

RsquareAdj(RDA) 
#adj.r.squared = 0.32 for order agglomerated ASVs - Robust Ait
anova.cca(RDA, permutations = 999) 

sqrt(vif.cca(RDA)) # sqrt(vif())>2 is considered highly collinear


# Ordistep to keep only relevant variables
fwd.sel <- ordiR2step(rda(t(ASV_object) ~ 1, data = sed_rda_data[,4:11]), 
                      scope = formula(RDA), # Complete model 
                      direction = "both",
                      R2scope = TRUE, # limited by the R2 of the complete modele 
                      #pstep = 1000,
                      trace = TRUE) # TRUE shows the selection process
fwd.sel

?ordiR2step
# robust Ait
# pH  for order agglomerated ASVs


RDA_signif <- vegan::rda(t(ASV_object) ~ pH + 
                             #`% Saturation` + 
                             `Potassium_mg/kg` + 
                             `Magnesium_mg/kg`,
                         data = sed_rda_data[,4:11])

sqrt(vif.cca(RDA_signif)) 
# sqrt(vif())>2 is considered highly collinear

RsquareAdj(RDA_signif) 
#0.03302412 for Rclr orde glom
anova.cca(RDA_signif, permutations = 999) 

sqrt(vif.cca(RDA_signif)) 

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
 Order_annotation <-  subset(coord.asv.rda,abs(coord.asv.rda$RDA1)>1 | abs(coord.asv.rda$RDA2)>1)
 tax <- as.data.frame(tax_table(ps_warm_rclr))
 Interest <- subset(tax,rownames(tax)%in%rownames(Order_annotation))
 Interest <- cbind(Interest,Order_annotation)  


Rda.plot.rclr.S1<-ggplot(data=coord.sites.rda, aes(x=RDA1, y=RDA2))+
    theme_bw()+
    geom_point(aes(shape=new_sdata$Plant.type,color=new_sdata$Temperature.x),size=2.8)+
    geom_point(data=coord.asv.rda3,aes(x=RDA1, y=RDA2),colour="purple4",shape=4, alpha=0.7)+
    geom_text_repel(data=Interest,label=Interest$Order,colour="purple4",
                    size=3.5,
                    min.segment.length = 0,
                    #arrow=arrow(angle = 30, length = unit(0.01, "inches"),ends = "first", type = "open"),
                    fontface = "italic")+ 
    geom_hline(yintercept=0, linetype="dotted") +  
    geom_vline(xintercept=0, linetype="dotted") +
    geom_segment(data= coord.env.mul, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
                 color="tomato4", arrow=arrow(length=unit(0.01,"npc")))+
    geom_text_repel(data=coord.env.mul, aes(label=rownames(coord.env.mul)),
                    color="tomato4",force=2, size=4,fontface="bold")+
    labs(x = paste0("RDA1 (",round(res.rda$cont$importance[2,1]*100,2),"%)"),#Variance explained by each axis
         y = paste0("RDA2 (",round(res.rda$cont$importance[2,2]*100,2),"%)"),#Variance explained by each axis
         title="RDA constrained by sediment chemistry",
         subtitle="Robust Aitchison transformed D62 Warm samples - Scaling 1")+
    geom_label_npc(data= df.annotations , 
                   aes(npcx = "left", npcy = "bottom", label = label),
                   parse=T,size=6)+
    scale_shape_manual(name="Sample type",values=c(21,22,23))+
    scale_color_manual(name="Temperature",values=c("red"),labels=c("Warm"))

Rda.plot.rclr.S1


#ggsave(here("Results_W&C","Figures","Sed_RDA_S1_D62_Warm_OrderGlom.pdf"),device='pdf',height = 7.5, width = 10.5)
#ggsave(here("Results_W&C","Figures","Sed_RDA_S1_D62_Warm_OrderGlom.png"),device='png',height = 7.5, width = 10.5)


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
coord.env.mul <- coord.env*5 # Multiplication factor for plotting purposes

### Labelling "extreme" taxa
Order_annotation <-  subset(coord.asv.rda,abs(coord.asv.rda$RDA1)>1 | abs(coord.asv.rda$RDA2)>1)
tax <- as.data.frame(tax_table(ps_warm_rclr))
Interest <- subset(tax,rownames(tax)%in%rownames(Order_annotation))
Interest <- cbind(Interest,Order_annotation)  


Rda.plot.rclr.S2<-ggplot(data=coord.sites.rda, aes(x=RDA1, y=RDA2))+
    theme_bw()+
    geom_point(aes(shape=new_sdata$Plant.type,color=new_sdata$Temperature.x),size=2.8)+
    geom_point(data=coord.asv.rda3,aes(x=RDA1, y=RDA2),colour="purple4",shape=4, alpha=0.7)+
    geom_text_repel(data=Interest,label=Interest$Order,colour="purple4",
                    size=3.5,
                    min.segment.length = 0,
                    #arrow=arrow(angle = 30, length = unit(0.01, "inches"),ends = "first", type = "open"),
                    fontface = "italic")+ 
    geom_hline(yintercept=0, linetype="dotted") +  
    geom_vline(xintercept=0, linetype="dotted") +
    geom_segment(data= coord.env.mul, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
                 color="tomato4", arrow=arrow(length=unit(0.01,"npc")))+
    geom_text_repel(data=coord.env.mul, aes(label=rownames(coord.env.mul)),
                    color="tomato4",force=2, size=4,fontface="bold")+
    labs(x = paste0("RDA1 (",round(res.rda$cont$importance[2,1]*100,2),"%)"),#Variance explained by each axis
         y = paste0("RDA2 (",round(res.rda$cont$importance[2,2]*100,2),"%)"),#Variance explained by each axis
         title="RDA constrained by sediment chemistry",
         subtitle="Robust Aitchison transformed D62 Warm samples - Scaling 2")+
    geom_label_npc(data= df.annotations , 
                   aes(npcx = "left", npcy = "bottom", label = label),
                   parse=T,size=6)+
    scale_shape_manual(name="Sample type",values=c(21,22,23))+
    scale_color_manual(name="Temperature",values=c("red"),labels=c("Warm"))

Rda.plot.rclr.S2

#ggsave(here("Results_W&C","Figures","Sed_RDA_S2_D62_Warm_OrderGlom.pdf"),device='pdf',height = 7.5, width = 10.5)
#ggsave(here("Results_W&C","Figures","Sed_RDA_S2_D62_Warm_OrderGlom.png"),device='png',height = 7.5, width = 10.5)

RDA_sed_warm_chem <- ggarrange(Rda.plot.rclr.S1,Rda.plot.rclr.S2,common.legend = T,legend = "bottom")

ggsave(here("Results_W&C","Figures","RDA_sed_warm_chem.pdf"),device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results_W&C","Figures","RDA_sed_warm_chem.png"),device='png',height = 7.5, width = 10.5)


#### COLD #####
#### | Merge datasets ####

# Retrieve sample data from phyloseq object
sdata <- data.frame(sample_data(ps_cold_order))
colnames(sdata)
# Keep only relevant taxa for right_join
sdata <- select(sdata,c("CFL_ID_Sediments","Mesocosm..","Plant.type","Temperature","Greenhouse"))

# Merge the datasets
new_sdata <- sdata %>% inner_join(sed_data, by=c("Mesocosm..","Greenhouse","Plant.type"))
colnames(new_sdata)

# Keep relevant columns
sed_rda_data <- select(new_sdata,c("CFL_ID_Sediments",
                                   "Plant.type",
                                   "Mesocosm..",
                                   "Temperature.x",
                                   "pH",
                                   #"Electrical Conductivity_ds/m", #correlated to other variables
                                   #"SAR",# Correlated to Chloride                     
                                   "% Saturation",                 
                                   #"Calcium_meq/L",                
                                   "Calcium_mg/kg" ,              
                                   #"Magnesium_meq/L",              
                                   "Magnesium_mg/kg",             
                                   #"Sodium_meq/L",                
                                   "Sodium_mg/kg",                 
                                   #"Potassium_meq/L",              
                                   "Potassium_mg/kg",             
                                   #"Chloride_meq/L",              
                                   "Chloride_mg/kg",               
                                   #"Sulfate(SO4)_meq/L",          
                                   "Sulfate(SO4)_mg/kg"))      
#"Sulfate-S_meq/L",   # Correlated to Sulfate SO4         
#"Sulfate-S_mg/kg")) 
#"TGR_T/ac"))  # below detection limit        


# Setting ID as rownames and removing column
rownames(sed_rda_data) <- sed_rda_data[,1]
sed_rda_data <- sed_rda_data[,-1]

# Formating the data 9numbers to numeric , characters to factor)
sed_rda_data[,4:11] <- sapply(sed_rda_data[,4:11], as.numeric)
sed_rda_data[sapply(sed_rda_data, is.character)] <- lapply(sed_rda_data[sapply(sed_rda_data, is.character)], as.factor) # did it work? Check with str(meta)

# Standardizing the data
sed_rda_data[,4:11] <- decostand(sed_rda_data[,4:11] , method = "standardize")

# Checking standaridization
# round(apply(sed_rda_data[,4:11], 2, mean), 1) # Variables are now centered around a mean of 0
# apply(sed_rda_data[,4:11], 2, sd) # and scaled to have a standard deviation of 1

heatmap(abs(cor(sed_rda_data[,4:11])), 
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("topright", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))
# Still some correlations (chloride and sodium ; sulfate and magnesium)


#### | RDA #####
ASV_object <- ASV_cold_rclr

RDA <- vegan::rda(t(ASV_object) ~ ., data = sed_rda_data[,4:11])

RsquareAdj(RDA) 
#adj.r.squared = 0.10 for order agglomerated ASVs - Robust Ait
anova.cca(RDA, permutations = 999) 

sqrt(vif.cca(RDA)) # sqrt(vif())>2 is considered highly collinear


# Ordistep to keep only relevant variables
fwd.sel <- ordiR2step(rda(t(ASV_object) ~ 1, data = sed_rda_data[,4:11]), 
                      scope = formula(RDA), # Complete model 
                      direction = "both",
                      R2scope = TRUE, # limited by the R2 of the complete modele 
                      #pstep = 1000,
                      trace = TRUE) # TRUE shows the selection process



# robust Ait
# pH  for order agglomerated ASVs


RDA_signif <- vegan::rda(t(ASV_object) ~  + 
                             `Sulfate(SO4)_mg/kg` + 
                             `Potassium_mg/kg` + 
                             `Calcium_mg/kg`,
                         data = sed_rda_data[,4:11])

sqrt(vif.cca(RDA_signif)) 
# sqrt(vif())>2 is considered highly collinear

RsquareAdj(RDA_signif) 
#0.03302412 for Rclr orde glom
anova.cca(RDA_signif, permutations = 999) 

sqrt(vif.cca(RDA_signif)) 

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
Order_annotation <-  subset(coord.asv.rda,abs(coord.asv.rda$RDA1)>1 | abs(coord.asv.rda$RDA2)>1)
tax <- as.data.frame(tax_table(ps_cold_rclr))
Interest <- subset(tax,rownames(tax)%in%rownames(Order_annotation))
Interest <- cbind(Interest,Order_annotation)  


Rda.plot.rclr.S1<-ggplot(data=coord.sites.rda, aes(x=RDA1, y=RDA2))+
    theme_bw()+
    geom_point(aes(shape=new_sdata$Plant.type,color=new_sdata$Temperature.x),size=2.8)+
    geom_point(data=coord.asv.rda3,aes(x=RDA1, y=RDA2),colour="purple4",shape=4, alpha=0.7)+
    geom_text_repel(data=Interest,label=Interest$Order,colour="purple4",
                    size=3.5,
                    min.segment.length = 0,
                    #arrow=arrow(angle = 30, length = unit(0.01, "inches"),ends = "first", type = "open"),
                    fontface = "italic")+ 
    geom_hline(yintercept=0, linetype="dotted") +  
    geom_vline(xintercept=0, linetype="dotted") +
    geom_segment(data= coord.env.mul, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
                 color="tomato4", arrow=arrow(length=unit(0.01,"npc")))+
    geom_text_repel(data=coord.env.mul, aes(label=rownames(coord.env.mul)),
                    color="tomato4",force=2, size=4,fontface="bold")+
    labs(x = paste0("RDA1 (",round(res.rda$cont$importance[2,1]*100,2),"%)"),#Variance explained by each axis
         y = paste0("RDA2 (",round(res.rda$cont$importance[2,2]*100,2),"%)"),#Variance explained by each axis
         title="RDA constrained by sediment chemistry",
         subtitle="Robust Aitchison transformed D62 cold samples - Scaling 1")+
    geom_label_npc(data= df.annotations , 
                   aes(npcx = "left", npcy = "bottom", label = label),
                   parse=T,size=6)+
    scale_shape_manual(name="Sample type",values=c(21,22,23))+
    scale_color_manual(name="Temperature",values=c("red"),labels=c("cold"))

Rda.plot.rclr.S1


#ggsave(here("Results_W&C","Figures","Sed_RDA_S1_D62_cold_OrderGlom.pdf"),device='pdf',height = 7.5, width = 10.5)
#ggsave(here("Results_W&C","Figures","Sed_RDA_S1_D62_cold_OrderGlom.png"),device='png',height = 7.5, width = 10.5)


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
coord.env.mul <- coord.env*5 # Multiplication factor for plotting purposes

### Labelling "extreme" taxa
Order_annotation <-  subset(coord.asv.rda,abs(coord.asv.rda$RDA1)>1 | abs(coord.asv.rda$RDA2)>1)
tax <- as.data.frame(tax_table(ps_cold_rclr))
Interest <- subset(tax,rownames(tax)%in%rownames(Order_annotation))
Interest <- cbind(Interest,Order_annotation)  


Rda.plot.rclr.S2<-ggplot(data=coord.sites.rda, aes(x=RDA1, y=RDA2))+
    theme_bw()+
    geom_point(aes(shape=new_sdata$Plant.type,color=new_sdata$Temperature.x),size=2.8)+
    geom_point(data=coord.asv.rda3,aes(x=RDA1, y=RDA2),colour="purple4",shape=4, alpha=0.7)+
    geom_text_repel(data=Interest,label=Interest$Order,colour="purple4",
                    size=3.5,
                    min.segment.length = 0,
                    #arrow=arrow(angle = 30, length = unit(0.01, "inches"),ends = "first", type = "open"),
                    fontface = "italic")+ 
    geom_hline(yintercept=0, linetype="dotted") +  
    geom_vline(xintercept=0, linetype="dotted") +
    geom_segment(data= coord.env.mul, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
                 color="tomato4", arrow=arrow(length=unit(0.01,"npc")))+
    geom_text_repel(data=coord.env.mul, aes(label=rownames(coord.env.mul)),
                    color="tomato4",force=2, size=4,fontface="bold")+
    labs(x = paste0("RDA1 (",round(res.rda$cont$importance[2,1]*100,2),"%)"),#Variance explained by each axis
         y = paste0("RDA2 (",round(res.rda$cont$importance[2,2]*100,2),"%)"),#Variance explained by each axis
         title="RDA constrained by sediment chemistry",
         subtitle="Robust Aitchison transformed D62 cold samples - Scaling 2")+
    geom_label_npc(data= df.annotations , 
                   aes(npcx = "left", npcy = "bottom", label = label),
                   parse=T,size=6)+
    scale_shape_manual(name="Sample type",values=c(21,22,23))+
    scale_color_manual(name="Temperature",values=c("red"),labels=c("cold"))

Rda.plot.rclr.S2

#ggsave(here("Results_W&C","Figures","Sed_RDA_S2_D62_cold_OrderGlom.pdf"),device='pdf',height = 7.5, width = 10.5)
#ggsave(here("Results_W&C","Figures","Sed_RDA_S2_D62_cold_OrderGlom.png"),device='png',height = 7.5, width = 10.5)

RDA_sed_cold_chem <- ggarrange(Rda.plot.rclr.S1,Rda.plot.rclr.S2,common.legend = T,legend = "bottom")

ggsave(here("Results_W&C","Figures","RDA_sed_cold_chem.pdf"),device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results_W&C","Figures","RDA_sed_cold_chem.png"),device='png',height = 7.5, width = 10.5)



