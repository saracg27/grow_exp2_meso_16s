load(file = here::here("RData", "neg_controls.RData"))


# -------------------------------------------------------
# DATA WRANGLING
# -------------------------------------------------------
## Order metadata ----
meta_sorted <- meta[order(row.names(meta)),] 
str(meta_sorted) # 298 in 11 variables
## Order counts ----
counts_sorted <- counts_asv_neg_control[,order(colnames(counts_asv_neg_control))]
row.names(meta_sorted) == colnames(counts_sorted) # three last in counts_sorted not the same: taxomy
## Aggregate ----
counts_sorted <- aggregate(. ~ taxonomy, data = counts_sorted, sum)
## Get taxonomy as row names ---- 
rownames(counts_sorted) <- counts_sorted$taxonomy
## Get only samples, no tax ----
tax <- counts_sorted$taxonomy
counts_sorted$taxonomy <- NULL
## Transpose ----
counts_t <- t(counts_sorted)
## Keep only samples in metadata that are present in counts_sorted ----
meta_sorted <- subset(meta_sorted, rownames(meta_sorted) %in% colnames(counts_sorted)) 
## Check row names order ----
row.names(meta_sorted) == row.names(counts_t) # all TRUE

## Rearrange factor levels in meta_sorted
levels(meta_sorted$Plant.type)
meta_sorted$Plant.type <- factor(meta_sorted$Plant.type, 
                                 levels = c("no plant",
                                            "Scirpus",
                                            "Triglochin",
                                            "Neg control extraction rhizosphere",
                                            "Neg control extraction roots",
                                            "Neg control extraction sediment",
                                            "Neg control PCR rhizosphere",
                                            "Neg control PCR roots",
                                            "Neg control PCR sediment"),
                                 labels = c("no plant",
                                            "Scirpus",
                                            "Triglochin",
                                            "ctrl_ext_rh",
                                            "ctrl_ext_rt",
                                            "ctrl_ext_sd",
                                            "ctrl_pcr_rh",
                                            "ctrl_pcr_rt",
                                            "ctrl_pcr_sd"))
## Clean names ----
genus_names <- colnames(counts_t)
genus_names <- gsub(".*f__", "", genus_names) 
genus_names <- gsub(".g__", "", genus_names) 
genus_names
genus_names <- sub("Beijerinckiaceae;Methylobacterium-Methylorubrum; s__Methylobacterium_aerolatum", "Methylobacterium aerolatum", genus_names)
genus_names <- sub("KD4-96;KD4-96; s__uncultured_Chloroflexi", "uncultured Chloroflexi", genus_names)
genus_names <- sub("uncultured;uncultured; s__uncultured_Rhodobacteraceae", "uncultured Rhodobacteraceae", genus_names)
genus_names <- sub("uncultured;uncultured; s__uncultured_Acidobacteriales", "uncultured Acidobacteriales", genus_names)
genus_names <- sub("Chitinophagaceae;Sediminibacterium; s__Sediminibacterium_sp.", "Sediminibacterium", genus_names)
genus_names <- sub("Gemmatimonadaceae;Roseisolibacter; s__uncultured_bacterium", "Roseisolibacter", genus_names)
genus_names <- sub("Polyangiaceae;Polyangium; s__metagenome", "Polyangium", genus_names)
genus_names <- sub("Beijerinckiaceae;uncultured", "uncultured Beijerinckiaceae", genus_names)
genus_names <- sub("Caulobacteraceae;uncultured", "uncultured Caulobacteraceae", genus_names)
genus_names <- sub("Micropepsaceae;uncultured", "uncultured Micropepsaceae", genus_names)
genus_names <- sub("Pirellulaceae;uncultured", "uncultured Pirellulaceae", genus_names)
genus_names <- sub("Gemmatimonadaceae;uncultured", "uncultured Gemmatimonadaceae", genus_names)
genus_names <- sub("Beijerinckiaceae;Methylobacterium-Methylorubrum", "Methylobacterium-Methylorubrum", genus_names)
genus_names <- sub("Nitrososphaeraceae;Nitrososphaeraceae; s__uncultured_euryarchaeote", "uncultured Nitrososphaeraceae", genus_names)
genus_names <- sub("Xiphinematobacteraceae;Candidatus_Xiphinematobacter; s__uncultured_Verrucomicrobia", "Candidatus Xiphinematobacter", genus_names)
genus_names <- sub("Obscuribacteraceae;Candidatus_Obscuribacter; s__uncultured_bacterium", "Candidatus Obscuribacter", genus_names)
genus_names <- sub(".*;", "", genus_names)
genus_names

colnames(counts_t) <- genus_names

## Bind treatments and class abundances ----
abund_treat <- cbind(meta_sorted,counts_t)
dim(abund_treat) # 297 - 51

## Gather in long format for ggplot ----
long <- gather(abund_treat,Taxa,RelAbund,12:ncol(abund_treat)) 
levels(factor(long$Taxa))
head(long) # It looks ok!


# -------------------------------------------------------
# STACKED BARS RAW COUNTS
# -------------------------------------------------------
pal40 <- c("#F0F8FF", "#EEDFCC", "#8B8378", "#7FFFD4", "#458B74", "#C1CDCD", "#838B8B", "#FFE4C4", "#000000", 
           "#0000FF", "#8A2BE2", "#A52A2A", "#FF4040", "#DEB887", "#5F9EA0", "#98F5FF", "#7AC5CD", "#7FFF00", 
           "#458B00", "#FF7F24", "#6495ED", "#00FFFF", "#FFB90F", "#CAFF70", "#BF3EFF", "#E9967A", "#8FBC8F", 
           "#FF1493", "#CD1076", "#00BFFF", "#104E8B", "#FFD700", "#F0FFF0", "#FFB5C5", "#CDBE70", "#FFA07A", 
           "#20B2AA", "#FFFF00", "#EE82EE", "#BFBFBF")
## Create the stack bar chart
stack <- ggplot(long, aes(fill = Taxa, y = RelAbund, x = Plant.type)) + 
    geom_bar(stat = "identity") +
    labs( y = "Raw counts of ASVs found in controls" ) + 
    theme_minimal() +
    theme(strip.text.x = element_text(face = "bold", size = 14),
          axis.text.x = element_text(size = 12, angle = 360, hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold", ),
          plot.title = element_text(face = "bold", size = 16),
          legend.title = element_text(size = 14, color = "black", face = "bold"),
          legend.text = element_text(size = 8, face = "bold", color ="black"),
          legend.position = "bottom",
          axis.title = element_text(size = 16, face = "bold"),
          rect = element_rect(fill = "white")) +
    scale_fill_manual(values = pal40, guide = guide_legend(label.theme = element_text(face = "italic", size = 12))) +
    scale_color_manual(values = pal40, guide = guide_legend(label.theme = element_text(face = "italic", size = 12))) +
    scale_x_discrete(name = "treatments")

stack

ggsave(file = here::here("output", "figures", "Raw_abund_control_asvs.pdf"), stack, width = 15, height = 7, units = "in", dpi = 300,bg = "white")

# -------------------------------------------------------
# STACKED BARS RELATIVE ABUNDANCE
# -------------------------------------------------------
stack <- ggplot(long, aes(fill = Taxa, y = RelAbund, x = Plant.type)) + 
    geom_bar(stat = "identity", position = "fill") +
    labs( y = "Relative abundance of ASVs found in controls" ) + 
    theme_minimal() +
    theme(strip.text.x = element_text(face = "bold", size = 14),
          axis.text.x = element_text(size = 12, angle = 360, hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold", ),
          plot.title = element_text(face = "bold", size = 16),
          legend.title = element_text(size = 14, color = "black", face = "bold"),
          legend.text = element_text(size = 8, face = "bold", color ="black"),
          legend.position = "bottom",
          axis.title = element_text(size = 16, face = "bold"),
          rect = element_rect(fill = "white")) +
    scale_fill_manual(values = pal40, guide = guide_legend(label.theme = element_text(face = "italic", size = 12))) +
    scale_color_manual(values = pal40, guide = guide_legend(label.theme = element_text(face = "italic", size = 12))) +
    #facet_grid(.~ Plant.type,  scales = "free_x", space = "free_x") +
    #facet_wrap(. ~ Temperature) +
    #scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
    scale_x_discrete(name = "treatments")

stack

ggsave(file = here::here("output", "figures", "all_rel_abund_control_asv.pdf"), stack, width = 15, height = 7, units = "in", dpi = 300,bg = "white")
# -------------------------------------------------------
# STACKED BARS MOST ABUNDANT RAW COUNTS
# -------------------------------------------------------
## Only most abundant: >300 counts 
## Remove low abundances
colSums(counts_t) 
abund <- counts_t[, colSums(counts_t) > 300] # 22 groups
colnames(abund)
abund_treat <- cbind(meta_sorted,abund)
long <- gather(abund_treat,Taxa,RelAbund,12:ncol(abund_treat)) 
head(long) # It looks ok!

pal3 <- c("#F2D696","#F59C9B","#71ACD6","#7B9E81","#CDC0B0","#968EC2","#B2D2E8",
          "#F7C6EC", "#C7DBC5","#FAFAD2","#C25B57","#A2C5F2", "#C6F7D0","#B88988",
          "#C3C5FA","#8CBCFA","#DBADD3","#E0EEEE","#FFD39B", "#E7D0F5","#C3F3FA","#008B8B")


stack <- ggplot(long, aes(fill = Taxa, y = RelAbund, x = Plant.type)) + 
    geom_bar(stat = "identity") +
    labs( y = "Counts of most abundant ASVs found in controls" ) + 
    theme_minimal() +
    theme(strip.text.x = element_text(face = "bold", size = 14),
          axis.text.x = element_text(size = 12, angle = 360, hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold", ),
          plot.title = element_text(face = "bold", size = 16),
          legend.title = element_text(size = 14, color = "black", face = "bold"),
          legend.text = element_text(size = 12, face = "italic", color ="black"),
          legend.position = "right",
          legend.box = "vertical",
          axis.title = element_text(size = 16, face = "bold"),
          rect = element_rect(fill = "white")) +
    scale_fill_manual(values = pal3, guide = guide_legend(label.theme = element_text(face = "italic", size = 12))) +
    scale_color_manual(values = pal3, guide = guide_legend(label.theme = element_text(face = "italic", size = 12))) +
    #facet_grid(.~ Plant.type,  scales = "free_x", space = "free_x") +
    #facet_wrap(. ~ Temperature) +
    #scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
    scale_x_discrete(name = "treatments") +
    guides(fill=guide_legend(ncol =1))

stack

ggsave(file = here::here("output", "figures", "Raw_abund_control_asvs_300plus.pdf"), stack, width = 12, height = 7, units = "in", dpi = 300,bg = "white")
# ggsave(file = here::here("output", "figures", "relAbund_phylum_16Swater.eps"), stack_phylum, width = 10, height = 7, units = "in", dpi = 300,bg = "white")

### It seems that there is not too many counts in the controls. We can proceed without correcting abundances in samples
