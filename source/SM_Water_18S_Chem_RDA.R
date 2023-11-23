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


### ASV Data ####
# load phyloseq object
load(here("Rdata","ps_18S_water_obj.RData"))

### Water chemistry data ####
water_data <- read.table(file = "Data/Env_measures/env_final_water_exp2_20220225.csv", 
                         dec = ".", sep = ",")
colnames(water_data) <- water_data[1,] # Rename colnames (similar column name in original table)
water_data <- water_data[-1,-c(2,11)] # removing first row (colnames) 
# and one of the Electrical conductivity columns and Temp
colnames(water_data)
colnames(water_data)[4] <- "Mesocosm"# Change name to match the colname in the sample data below 



#### ~~~~~~~~######

### D77 ######
ps_77 <- subset_samples(ps, Time =="D77")
ps_77 <- prune_taxa(taxa_sums(ps_77)>0, ps_77)
# 24 samples 285 taxa

#### RCLR transformation #####

ps_77order <- tax_glom(ps_77,taxrank="Order",NArm=F) #109 orders
ps.order.rclr <- microbiome::transform(ps_77order, "rclr") 

####  Merge datasets ####
sdata <- data.frame(sample_data(ps_77order)) # Extract dataset from phyloseq object

colnames(sdata)

# We need to have matching columns between sdata and water data 
sdata$Mesocosm <- sub("\\_.*", "", sdata$alias)
sdata$Mesocosm <- sub("M", "Mesocosm_",sdata$Mesocosm)

sdata$Greenhouse <- str_sub(sdata$Mesocosm,-1,-1)
sdata$Mesocosm <- str_sub(sdata$Mesocosm,end=-2)

sdata$Temperature <- sub("10C day - 5C night","cold",sdata$Temperature)
sdata$Temperature <- sub("20C day - 10C night","warm",sdata$Temperature)

colnames(sdata)[which(names(sdata) == "Sample_type")] <- "Plant.type"


# Match factor name with sdata
water_data$Plant.type<- sub("no plant", "No_plant", water_data$Plant.type)

## Just take the necessary columns
sdata <- select(sdata,c("Mesocosm","Plant.type","Temperature","Greenhouse"))

# Merge the sample data and water chemistry data using several matching columns 
new_sdata <- sdata %>% right_join(water_data, by=c("Mesocosm","Greenhouse","Plant.type","Temperature"))
colnames(new_sdata)

### Subset resulting merged dataset to keep relevant columns
water_rda_data <- select(new_sdata,c("Mesocosm",
                                     "Plant.type",
                                     "Temperature",
                                     "pH",
                                     #"Temperature of observed pH", #Not relevant
                                     "Electrical Conductivity",   
                                     "Calcium",   
                                     "Magnesium",                 
                                     "Sodium",                     
                                     "Potassium",                 
                                     #"Iron",  #below detection limit                     
                                     #"Manganese", #below detection limit              
                                     "Chloride",                   
                                     #"Nitrate - N",  # most below detection limit              
                                     #"Nitrite - N",   #below detection limit             
                                     #"Nitrate and Nitrite - N",   # most below detection limit 
                                     "Sulfate (SO4)",              
                                     #"Hydroxide",    #below detection limit             
                                     #"Carbonate",    #below detection limit              
                                     "Bicarbonate",               
                                     #"P-Alkalinity",   #below detection limit          
                                     "T-Alkalinity",            
                                     "Total Dissolved Solids",    
                                     "Hardness",                  
                                     "Ionic Balance"))



#### Formatting data ####

# Formatting the chemistry data
water_rda_data[,4:16] <- sapply(water_rda_data[,4:16], as.numeric)
water_rda_data[sapply(water_rda_data, is.character)] <- lapply(water_rda_data[sapply(water_rda_data, is.character)], as.factor) 
# did it work? Check with str(water_rda_data)
# Standardizing the chemistry data
water_rda_data[,4:16] <- decostand(water_rda_data[,4:16] , method = "standardize")

# Checking standaridization
# round(apply(water_rda_data[,4:16], 2, mean), 1) # Variables are now centered around a mean of 0
# apply(water_rda_data[,4:16], 2, sd) # and scaled to have a standard deviation of 1

# Heatmap to chekc colinearity between variables
heatmap(abs(cor(water_rda_data[,4:16])), 
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("topright", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))
# T alkanity and bicarbonate are correlated
# hardeness and dissolved solids as well







#### RDA #####
ASV_order_rclr <- data.frame(otu_table(ps.order.rclr))

# Select the abundance dataset to use
ASV_object <- ASV_order_rclr

# Compute the rda using chemistry data as constraining variable 
RDA <- vegan::rda(t(ASV_object) ~ ., data = water_rda_data[,4:16])


RsquareAdj(RDA) # adj.r.squared = x for un-agglomerated ASVs
# adj.r.squared = x for genus agglomerated ASVs
# adj.r.squared = 0.15 for order agglomerated ASVs


# Forward selection of variables with the most explanating power 
fwd.sel <- ordiR2step(rda(t(ASV_object) ~ 1, data = water_rda_data[,4:16]), # modèle le plus simple
                      scope = formula(RDA), # modèle "complet"
                      direction = "forward",
                      R2scope = TRUE, # limité par le R2 du modèle "complet"
                      #pstep = 1000,
                      trace = TRUE) # mettre TRUE pour voir le processus du sélection!
fwd.sel$call

# Does not find a more parsimonous model

# Checking the significance of variables
anova.cca(RDA, permutations = 1000)
anova.cca(RDA, permutations = 1000,by="term")
# Magnesium and Calcium  and Total Dissolved Solid for order agglomerated

# Variance inflation factors function to identify useless constraining variables
sqrt(vif.cca(RDA)) # sqrt(vif())>2 is considered highly collinear


RDA_signif <- vegan::rda(t(ASV_object) ~ pH  + Potassium + Chloride + Magnesium  , data = water_rda_data[,4:16])
RsquareAdj(RDA_signif)
sqrt(vif.cca(RDA_signif))
anova.cca(RDA_signif, permutations = 1000,by="term")
anova.cca(RDA_signif, permutations = 1000)

##### Scaling 1 #####
res.rda<-summary(RDA_signif,scaling = 1) # SCALING 1 
head(res.rda) # Check

# ASV coordinates
coord.asv.rda<-as.data.frame(res.rda$species) 

# Change rownames from ASV sequence to ASV# for eventual plotting purposes
coord.asv.rda2 <- cbind(paste0("ASV",seq(1:nrow(coord.asv.rda))),coord.asv.rda) # Change ASV sequence by ASV#
coord.asv.rda3 <- coord.asv.rda2
row.names(coord.asv.rda3) <- coord.asv.rda2[,1]
coord.asv.rda3 <- coord.asv.rda3[,-1]

# Site coordinates
coord.sites.rda<-as.data.frame(res.rda$sites)
#Renaming sample names for plotting purposes
rownames(coord.sites.rda) <- gsub("\\..*","",rownames(coord.sites.rda))


# Environmental data coordinates
coord.env <- as.data.frame(res.rda$biplot)
coord.env.mul <- coord.env*5 # Multiplication factor for plotting purposes


## Annotations for R2 and P in graph
df.annotations <- data.frame(
    label = paste(paste0("~italic(R)^{2} == ", round(RsquareAdj(RDA_signif)$adj.r.squared,3)),"~~",
                  paste0("~italic(p) == ",round(anova.cca(RDA_signif, permutations = 1000)$'Pr(>F)'[1],3))))


### Labelling "extreme" taxa for plotting purposes
Order_annotation <-  subset(coord.asv.rda,abs(coord.asv.rda$RDA1)>1 | abs(coord.asv.rda$RDA2)>1)
tax <- as.data.frame(tax_table(ps.order.rclr))
Interest <- subset(tax,rownames(tax)%in%rownames(Order_annotation))
Interest <- cbind(Interest,Order_annotation)  


# RDA plot # 
Rda.plot.rclr.S1<-ggplot(data=coord.sites.rda, aes(x=RDA1, y=RDA2))+
    theme_bw()+
    geom_point(aes(shape=new_sdata$Plant.type,color=new_sdata$Temperature),size=2.8)+
    #geom_text_repel(label=rownames(coord.sites.rda))+
    geom_point(data=coord.asv.rda,aes(x=RDA1, y=RDA2),colour="purple4",shape=4, alpha=0.7)+
    geom_text_repel(data=Interest,label=Interest$Order,colour="purple4",
                    size=3.5,
                    min.segment.length = 0,
                    arrow=arrow(angle = 30, length = unit(0.01, "inches"),ends = "first", type = "open"),
                    fontface = "italic")+ 
    geom_hline(yintercept=0, linetype="dotted") +  
    geom_vline(xintercept=0, linetype="dotted") +
    geom_segment(data= coord.env.mul, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
                 color="tomato4", arrow=arrow(length=unit(0.01,"npc")))+
    geom_text_repel(data=coord.env.mul, aes(label=rownames(coord.env.mul)),
                    color="tomato4",force=2, size=4,fontface="bold")+
    labs(x = paste0("RDA1 (",round(res.rda$cont$importance[2,1]*100,2),"%)"),#Variance explained by each axis
         y = paste0("RDA2 (",round(res.rda$cont$importance[2,2]*100,2),"%)"),#Variance explained by each axis
         title="RDA constrained by water chemistry",
         subtitle="Robust Aitchison transformed abundance matrix as response variables - Scaling 1")+
    geom_label_npc(data= df.annotations , 
                   aes(npcx = "right", npcy = "bottom", label = label),
                   parse=T,size=6)+ 
    scale_shape_manual(name="Sample type",values=c(21,22,23))+
    scale_color_manual(name="Temperature",values=c("blue","red"),labels=c("Cold","Warm"))

Rda.plot.rclr.S1

ggsave(here("Results","Figures","18S_Water_RDA_S1_D77_WaterChem_RCLR_OrderGlom.pdf"),device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results","Figures","18S_Water_RDA_S1_D77_WaterChem_RCLR_OrderGlom.png"),device='png',height = 7.5, width = 10.5)


##### Scaling 2 #####
res.rda<-summary(RDA_signif,scaling = 2) # SCALING 2 
head(res.rda)
res.rda$cont$importance
coord.asv.rda<-as.data.frame(res.rda$species) # ASV coordinates
coord.asv.rda2 <- cbind(paste0("ASV",seq(1:nrow(coord.asv.rda))),coord.asv.rda)
coord.asv.rda3 <- coord.asv.rda2
row.names(coord.asv.rda3) <- coord.asv.rda2[,1]
coord.asv.rda3 <- coord.asv.rda3[,-1]

coord.sites.rda<-as.data.frame(res.rda$sites)# Site coordinates
rownames(coord.sites.rda) <- gsub("\\..*","",rownames(coord.sites.rda))


coord.env <- as.data.frame(res.rda$biplot) # Environmental data coordinates


## Annotations for R2 and P in graph
df.annotations <- data.frame(
    label = paste(paste0("~italic(R)^{2} == ", round(RsquareAdj(RDA_signif)$adj.r.squared,3)),"~~",
                  paste0("~italic(p) == ",round(anova.cca(RDA_signif, permutations = 1000)$'Pr(>F)'[1],3))))


### Labelling "extreme" taxa
Order_annotation <-  subset(coord.asv.rda,abs(coord.asv.rda$RDA1)>0.5 | abs(coord.asv.rda$RDA2)>0.5)
tax <- as.data.frame(tax_table(ps.order.rclr))
Interest <- subset(tax,rownames(tax)%in%rownames(Order_annotation))
Interest <- cbind(Interest,Order_annotation)  



Rda.plot.rclr.S2<-ggplot(data=coord.sites.rda, aes(x=RDA1, y=RDA2))+ 
    theme_bw()+
    geom_point(aes(shape=new_sdata$Plant.type,color=new_sdata$Temperature),size=2.8)+
    #geom_text_repel(label=rownames(coord.sites.rda))+
    geom_point(data=coord.asv.rda,aes(x=RDA1, y=RDA2),colour="purple4",shape=4, alpha=0.7)+
    geom_text_repel(data=Interest,label=Interest$Order,colour="purple4",
                    size=3.5,
                    min.segment.length = 0,
                    arrow=arrow(angle = 30, length = unit(0.01, "inches"),ends = "first", type = "open"),
                    fontface = "italic")+ 
    geom_hline(yintercept=0, linetype="dotted") +  
    geom_vline(xintercept=0, linetype="dotted") +
    geom_segment(data= coord.env.mul, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
                 color="tomato4", arrow=arrow(length=unit(0.01,"npc")))+
    geom_text_repel(data=coord.env.mul, aes(label=rownames(coord.env.mul)),
                    color="tomato4",force=2, size=4,fontface="bold")+
    labs(x = paste0("RDA1 (",round(res.rda$cont$importance[2,1]*100,2),"%)"),#Variance explained by each axis
         y = paste0("RDA2 (",round(res.rda$cont$importance[2,2]*100,2),"%)"),#Variance explained by each axis
         title="RDA constrained by water chemistry",
         subtitle="Robust Aitchison transformed abundance matrix as response variables - Scaling 2")+
    geom_label_npc(data= df.annotations , 
                   aes(npcx = "right", npcy = "bottom", label = label),
                   parse=T,size=6)+ 
    scale_shape_manual(name="Sample type",values=c(21,22,23))+
    scale_color_manual(name="Temperature",values=c("blue","red"),labels=c("Cold","Warm"))


Rda.plot.rclr.S2


ggsave(here("Results","Figures","18S_Water_RDA_S2_D77_WaterChem_RCLR_OrderGlom.pdf"),device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results","Figures","18S_Water_RDA_S2_D77_WaterChem_RCLR_OrderGlom.png"),device='png',height = 7.5, width = 10.5)

### did not run d62 ##
### 
### 
### 
### 
### 
###¸¸¸¸¸¸¸¸¸¸¸¸¸####

### D62 ######
ps_62 <- subset_samples(ps, Time =="D62")
ps_62 <- prune_taxa(taxa_sums(ps_62)>0, ps_62)
# 23 samples 1244 taxa

#### RCLR transformation #####

ps_62order <- tax_glom(ps_62,taxrank="Order",NArm=F)
ps.order.rclr <- microbiome::transform(ps_62order, "rclr") 

####  Merge datasets ####
sdata <- data.frame(sample_data(ps_62order)) # Extract dataset from phyloseq object

colnames(sdata)

# We need to have matching columns between sdata and water data 
sdata$Mesocosm <- sub("\\_.*", "", sdata$alias)
sdata$Mesocosm <- sub("M", "Mesocosm_",sdata$Mesocosm)

sdata$Greenhouse <- str_sub(sdata$Mesocosm,-1,-1)
sdata$Mesocosm <- str_sub(sdata$Mesocosm,end=-2)

sdata$Temperature <- sub("10C day - 5C night","cold",sdata$Temperature)
sdata$Temperature <- sub("20C day - 10C night","warm",sdata$Temperature)

colnames(sdata)[which(names(sdata) == "Sample_type")] <- "Plant.type"

# Drop sample Mesocosm-10E-16 to match number of samples in sdata
water_data <- water_data[-which(water_data$ID == "Mesocosm-10E-16"),]

# Match factor name with sdata
water_data$Plant.type<- sub("no plant", "No_plant", water_data$Plant.type)

## Just take the necessary columns
sdata <- select(sdata,c("Mesocosm","Plant.type","Temperature","Greenhouse"))

# Merge the sample data and water chemistry data using several matching columns 
new_sdata <- sdata %>% right_join(water_data, by=c("Mesocosm","Greenhouse","Plant.type","Temperature"))
colnames(new_sdata)

### Subset resulting merged dataset to keep relevant columns
water_rda_data <- select(new_sdata,c("Mesocosm",
                                     "Plant.type",
                                     "Temperature",
                                     "pH",
                                     #"Temperature of observed pH", #Not relevant
                                     "Electrical Conductivity",   
                                     "Calcium",   
                                     "Magnesium",                 
                                     "Sodium",                     
                                     "Potassium",                 
                                     #"Iron",  #below detection limit                     
                                     #"Manganese", #below detection limit              
                                     "Chloride",                   
                                     #"Nitrate - N",  # most below detection limit              
                                     #"Nitrite - N",   #below detection limit             
                                     #"Nitrate and Nitrite - N",   # most below detection limit 
                                     "Sulfate (SO4)",              
                                     #"Hydroxide",    #below detection limit             
                                     #"Carbonate",    #below detection limit              
                                     "Bicarbonate",               
                                     #"P-Alkalinity",   #below detection limit          
                                     "T-Alkalinity",            
                                     "Total Dissolved Solids",    
                                     "Hardness",                  
                                     "Ionic Balance"))



#### Formatting data ####

# Formatting the chemistry data
water_rda_data[,4:16] <- sapply(water_rda_data[,4:16], as.numeric)
water_rda_data[sapply(water_rda_data, is.character)] <- lapply(water_rda_data[sapply(water_rda_data, is.character)], as.factor) 
# did it work? Check with str(water_rda_data)
# Standardizing the chemistry data
water_rda_data[,4:16] <- decostand(water_rda_data[,4:16] , method = "standardize")

# Checking standaridization
# round(apply(water_rda_data[,4:16], 2, mean), 1) # Variables are now centered around a mean of 0
# apply(water_rda_data[,4:16], 2, sd) # and scaled to have a standard deviation of 1

# Heatmap to chekc colinearity between variables
heatmap(abs(cor(water_rda_data[,4:16])), 
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("topright", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))
# T alkanity and bicarbonate are correlated
# hardeness and dissolved solids as well







#### RDA #####
ASV_order_rclr <- data.frame(otu_table(ps.order.rclr))

# Select the abundance dataset to use
ASV_object <- ASV_order_rclr

# Compute the rda using chemistry data as constraining variable 
RDA <- vegan::rda(t(ASV_object) ~ ., data = water_rda_data[,4:16])


RsquareAdj(RDA) # adj.r.squared = x for un-agglomerated ASVs
# adj.r.squared = x for genus agglomerated ASVs
# adj.r.squared = 0.1323821 for order agglomerated ASVs


# Forward selection of variables with the most explanating power 
fwd.sel <- ordiR2step(rda(t(ASV_object) ~ 1, data = water_rda_data[,4:16]), # modèle le plus simple
                      scope = formula(RDA), # modèle "complet"
                      direction = "forward",
                      R2scope = TRUE, # limité par le R2 du modèle "complet"
                      #pstep = 1000,
                      trace = TRUE) # mettre TRUE pour voir le processus du sélection!
fwd.sel$call
## Robust Aitchison
# Chloride + Potassium + pH for order agglomerated ASVs

# Checking the significance of variables
anova.cca(RDA, permutations = 1000)
anova.cca(RDA, permutations = 1000,by="term")

# Calcium and pH for order agglomerated


# Reducing the number of variables keeping the ones identified by ordiplot
RDA_signif <- vegan::rda(t(ASV_object) ~ Chloride + pH+ Potassium , data = water_rda_data[,4:16])
RsquareAdj(RDA_signif) #adj rsquared = 0.09 
anova.cca(RDA_signif, permutations = 1000,by="term")

# Variance inflation factors function to identify useles constraining variables
sqrt(vif.cca(RDA_signif)) # sqrt(vif())>2 is considered highly collinear


##### Scaling 1 #####
res.rda<-summary(RDA_signif,scaling = 1) # SCALING 1 
head(res.rda) # Check

# ASV coordinates
coord.asv.rda<-as.data.frame(res.rda$species) 

# Change rownames from ASV sequence to ASV# for eventual plotting purposes
coord.asv.rda2 <- cbind(paste0("ASV",seq(1:nrow(coord.asv.rda))),coord.asv.rda) # Change ASV sequence by ASV#
coord.asv.rda3 <- coord.asv.rda2
row.names(coord.asv.rda3) <- coord.asv.rda2[,1]
coord.asv.rda3 <- coord.asv.rda3[,-1]

# Site coordinates
coord.sites.rda<-as.data.frame(res.rda$sites)
#Renaming sample names for plotting purposes
rownames(coord.sites.rda) <- gsub("\\..*","",rownames(coord.sites.rda))


# Environmental data coordinates
coord.env <- as.data.frame(res.rda$biplot)
coord.env.mul <- coord.env*5 # Multiplication factor for plotting purposes


## Annotations for R2 and P in graph
df.annotations <- data.frame(
    label = paste(paste0("~italic(R)^{2} == ", round(RsquareAdj(RDA_signif)$adj.r.squared,3)),"~~",
                  paste0("~italic(p) == ",round(anova.cca(RDA_signif, permutations = 1000)$'Pr(>F)'[1],3))))


### Labelling "extreme" taxa for plotting purposes
Order_annotation <-  subset(coord.asv.rda,abs(coord.asv.rda$RDA1)>1 | abs(coord.asv.rda$RDA2)>1)
tax <- as.data.frame(tax_table(ps.order.rclr))
Interest <- subset(tax,rownames(tax)%in%rownames(Order_annotation))
Interest <- cbind(Interest,Order_annotation)  


# RDA plot # 
Rda.plot.rclr.S1<-ggplot(data=coord.sites.rda, aes(x=RDA1, y=RDA2))+
    theme_bw()+
    geom_point(aes(shape=new_sdata$Plant.type,color=new_sdata$Temperature),size=2.8)+
    #geom_text_repel(label=rownames(coord.sites.rda))+
    geom_point(data=coord.asv.rda,aes(x=RDA1, y=RDA2),colour="purple4",shape=4, alpha=0.2)+
    geom_text_repel(data=Interest,label=Interest$Order,colour="purple4",size=3.5,fontface = "italic")+ 
    geom_hline(yintercept=0, linetype="dotted") +  
    geom_vline(xintercept=0, linetype="dotted") +
    geom_segment(data= coord.env.mul, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
                 color="tomato4", arrow=arrow(length=unit(0.01,"npc")))+
    geom_text_repel(data=coord.env.mul, aes(label=rownames(coord.env.mul)),
                    color="tomato4", size=3.8)+
    labs(x = paste0("RDA1 (",round(res.rda$cont$importance[2,1]*100,2),"%)"),#Variance explained by each axis
         y = paste0("RDA2 (",round(res.rda$cont$importance[2,2]*100,2),"%)"),#Variance explained by each axis
         title="RDA constrained by water chemistry",
         subtitle="Robust Aitchison transformed abundance matrix as response variables - Scaling 1")+
    geom_label_npc(data= df.annotations , 
                   aes(npcx = "right", npcy = "bottom", label = label),
                   parse=T,size=6)+ 
    scale_shape_manual(name="Sample type",values=c(21,22,23))+
    scale_color_manual(name="Temperature",values=c("blue","red"),labels=c("Cold","Warm"))


Rda.plot.rclr.S1


ggsave(here("Results","Figures","16S_Water_RDA_S1_D62_WaterChem_RCLR_OrderGlom.pdf"),device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results","Figures","16S_Water_RDA_S1_D62_WaterChem_RCLR_OrderGlom.png"),device='png',height = 7.5, width = 10.5)


##### Scaling 2 #####
res.rda<-summary(RDA_signif,scaling = 2) # SCALING 2 
head(res.rda)
res.rda$cont$importance
coord.asv.rda<-as.data.frame(res.rda$species) # ASV coordinates
coord.asv.rda2 <- cbind(paste0("ASV",seq(1:nrow(coord.asv.rda))),coord.asv.rda)
coord.asv.rda3 <- coord.asv.rda2
row.names(coord.asv.rda3) <- coord.asv.rda2[,1]
coord.asv.rda3 <- coord.asv.rda3[,-1]

coord.sites.rda<-as.data.frame(res.rda$sites)# Site coordinates
rownames(coord.sites.rda) <- gsub("\\..*","",rownames(coord.sites.rda))


coord.env <- as.data.frame(res.rda$biplot) # Environmental data coordinates


## Annotations for R2 and P in graph
df.annotations <- data.frame(
    label = paste(paste0("~italic(R)^{2} == ", round(RsquareAdj(RDA_signif)$adj.r.squared,3)),"~~",
                  paste0("~italic(p) == ",round(anova.cca(RDA_signif, permutations = 1000)$'Pr(>F)'[1],3))))


### Labelling "extreme" taxa
Order_annotation <-  subset(coord.asv.rda,abs(coord.asv.rda$RDA1)>0.7 | abs(coord.asv.rda$RDA2)>0.7)
tax <- as.data.frame(tax_table(ps.order.rclr))
Interest <- subset(tax,rownames(tax)%in%rownames(Order_annotation))
Interest <- cbind(Interest,Order_annotation)  



Rda.plot.rclr.S2<-ggplot(data=coord.sites.rda, aes(x=RDA1, y=RDA2)) +
    theme_bw()+
    geom_point(aes(shape=new_sdata$Plant.type,color=new_sdata$Temperature),size=2.8)+
    #geom_text_repel(label=rownames(coord.sites.rda))+
    geom_point(data=coord.asv.rda,aes(x=RDA1, y=RDA2),colour="purple4",shape=4, alpha=0.2)+
    geom_text_repel(data=Interest,label=Interest$Order,colour="purple4",size=3.5,fontface = "italic")+ 
    geom_hline(yintercept=0, linetype="dotted") +  
    geom_vline(xintercept=0, linetype="dotted") +
    geom_segment(data= coord.env.mul, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
                 color="tomato4", arrow=arrow(length=unit(0.01,"npc")))+
    geom_text_repel(data=coord.env.mul, aes(label=rownames(coord.env.mul)),
                    color="tomato4", size=3.8)+
    labs(x = paste0("RDA1 (",round(res.rda$cont$importance[2,1]*100,2),"%)"),#Variance explained by each axis
         y = paste0("RDA2 (",round(res.rda$cont$importance[2,2]*100,2),"%)"),#Variance explained by each axis
         title="RDA constrained by water chemistry",
         subtitle="Robust Aitchison transformed abundance matrix as response variables - Scaling 2")+
    geom_label_npc(data= df.annotations , 
                   aes(npcx = "right", npcy = "bottom", label = label),
                   parse=T,size=6)+ 
    scale_shape_manual(name="Sample type",values=c(21,22,23))+
    scale_color_manual(name="Temperature",values=c("blue","red"),labels=c("Cold","Warm"))


Rda.plot.rclr.S2


ggsave(here("Results","Figures","16S_Water_RDA_S2_D62_WaterChem_RCLR_OrderGlom.pdf"),device='pdf',height = 7.5, width = 10.5)
ggsave(here("Results","Figures","16S_Water_RDA_S2_D62_WaterChem_RCLR_OrderGlom.png"),device='png',height = 7.5, width = 10.5)



























