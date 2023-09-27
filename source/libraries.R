# Common libraries

# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install_packages('BiocManager')

BiocManager::install('EnhancedVolcano')
BiocManager::install("ComplexHeatmap")

#install.packages("pacman")
pacman::p_load(pacman, colorspace, cowplot, DESeq2, dunn.test, glue,
               ggdist, ggforce, ggpubr, ggtext, here, labdsv, MASS, psych, 
               tidyverse, roperators, vegan, EnhancedVolcano, ComplexHeatmap)
