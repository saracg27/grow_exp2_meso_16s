## Packages
library(here)
source(here("source", "libraries.R"))

## Metadata
meta <- read.csv(file = here("data","clean", "metadata.tsv"), dec = ".", header = T, row.names = 1, sep = "\t", comment.char = "") #load
str(meta)
## Convert character vectors to factors
meta[sapply(meta, is.character)] <- lapply(meta[sapply(meta, is.character)], as.factor) # did it work? Check with str(meta)
rownames(meta) <- stringr::str_replace_all(rownames(meta), '[-]', '.')

## Raw filtered counts
raw <- read.csv(file = here("data", "clean", "feature_table_filtered.tsv"), dec = ".", sep = "\t", header = T, row.names = 1, comment.char = "") #20078 objs in 299 vars
names(raw) <- sub("X", "", names(raw))
#cols_raw <- as.data.frame(names(raw)) #inspecting cols show a Sequence column and a  taxonomy column.
raw$Sequence <- NULL
dim(raw)

## Identify controls in meta
controls_meta <- meta[meta$Experiment != "Exp2",]
controls <- rownames(controls_meta)
# Get only controls in raw data
controls_raw <- raw[,colnames(raw) %in% rownames(controls_meta)]
# Identify ASV that are ONLY present in controls
controls_raw <- controls_raw[rowSums(controls_raw) > 0,]
# Get subset of ASVs from controls in the raw dataset
asv_neg_control <- rownames(controls_raw)
counts_asv_neg_control <- raw[rownames(raw) %in% asv_neg_control, ] # working dataset

# tidy environment
rm(raw)
save.image(file = here::here("RData", "neg_controls.RData"))
