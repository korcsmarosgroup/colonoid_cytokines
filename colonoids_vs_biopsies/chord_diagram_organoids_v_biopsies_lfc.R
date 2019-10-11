# Script to create chord diagram to compare degs from organoids and biopsies
# fdr 0.01 and lfc 1.5 biopsy data after ensembl conversion
#
# Input:
#       all_org_files - File paths for the preprocessed (during 'network_generation.py') differential expression data for each 
#                       of the cytokine treated organoids. With adj p value filter applied.
#       biop_file - Excel file of preprocessed biopsy differential expression (provided by Polchronis)
#       cols - Text file containing the colour (hex code) for each cytokine group eg. 'Il13,tnfa'. Group in column 'category' and hex code in column 'colour'.
#
# Output: 
#     Text files saved of biopsy degs which are in the colonoid datasets (1 for UC and 1 for CD)
#     Script creates chord diagrams (one for UC and one for CD) but they must be saved manually.

##### Set up #####

library(circlize)
library(tidyr)
library(stringr)
library(readxl)
library(org.Hs.eg.db)
library(dplyr) # Must load after Biomart

setwd("../input_data/differential_expression_datasets/colonoids_processed/")
#all_org_files <- list.files(pattern = "_degs_adjp_0.01_lfc_1.txt$")
all_org_files <- list.files(pattern = "_degs_adjp_0.01.txt$")
cols <- "../../plot_colours_all.txt"
biop_file <- "../biopsies/GSE16879DEGsUCandcCDvsHCGEO2R.xlsx"

colours <- read.delim(cols, header = TRUE)

##### Process organoid DEGs #####

# Import deg lists
idx <- 0
for (file in all_org_files) {
  
  idx <- idx+1
  
  # Open file
  degs <- read.delim(file, sep = "\t",header=TRUE,na.strings=c("","NA"))
  
  # Remove columns, filter for lfc >= 1.5
  #degs <- degs  %>% filter((logFC >1.5) | (logFC < -1.5)) %>% dplyr::select(ENSEMBL.ID, SYMBOL)
  degs <- degs %>% dplyr::select(ENSEMBL.ID, SYMBOL)
  
  # Get cytokine name
  #cytokine <- str_replace(file, "_degs_adjp_0.01_lfc_1.txt", "")
  cytokine <- str_replace(file, "_degs_adjp_0.01.txt", "")
  degs <- degs %>% mutate(cytokine = cytokine) 
  
  #join degs together
  if (idx == 1){
    all_degs <- degs
  } else {
    all_degs <- full_join(all_degs, degs, by=c("ENSEMBL.ID", "SYMBOL"))
  }
  
}

# Join the cytokine columns together
all_degs2 <- unite(all_degs, cytokines, 3:6, sep="; ") %>%
    mutate(cytokines = str_replace_all(cytokines, '; NA,?','')) %>%
    mutate(cytokines = str_replace_all(cytokines, 'NA,?; ',''))

# Count numbers of DEGs in each category (collapse)
degs_count <- all_degs2 %>% 
  group_by(cytokines) %>% 
  summarise(org_count = n())

rm(all_degs, degs,all_org_files)

##### Process biopsy DEGs #####

# GEO file
cd_degs <- read_excel(biop_file, sheet ="cCDvsHC")
uc_degs <- read_excel(biop_file, sheet ="UCvsHC")
# Seperate gene symbols
cd_degs <- separate_rows(cd_degs,Gene.symbol,sep="///")
uc_degs <- separate_rows(uc_degs,Gene.symbol,sep="///")
# Filter for p adj <= 0.01 and lfc >= 1.5, get gene names
#cd_degs <- cd_degs %>% filter(adj.P.Val <= 0.01 & (logFC >= 1.5 | logFC <= -1.5)) %>% dplyr::select(c(Gene.symbol)) %>% mutate(in_cd = "TRUE") %>% unique() %>% drop_na()
#uc_degs <- uc_degs %>% filter(adj.P.Val <= 0.01 & (logFC >= 1.5 | logFC <= -1.5)) %>% dplyr::select(c(Gene.symbol)) %>% mutate(in_uc = "TRUE") %>% unique() %>% drop_na()
cd_degs <- cd_degs %>% filter(adj.P.Val <= 0.01 ) %>% dplyr::select(c(Gene.symbol)) %>% mutate(in_cd = "TRUE") %>% unique() %>% drop_na()
uc_degs <- uc_degs %>% filter(adj.P.Val <= 0.01 ) %>% dplyr::select(c(Gene.symbol)) %>% mutate(in_uc = "TRUE") %>% unique() %>% drop_na()


# Convert gene symbol to ensembl ID
row.names(cd_degs) <- cd_degs$Gene.symbol
cd_degs$Ensembl <- mapIds(org.Hs.eg.db, keys = row.names(cd_degs), keytype = "SYMBOL", column="ENSEMBL",multiVals="first")
row.names(uc_degs) <- uc_degs$Gene.symbol
uc_degs$Ensembl <- mapIds(org.Hs.eg.db, keys = row.names(uc_degs), keytype = "SYMBOL", column="ENSEMBL",multiVals="first")
cd_degs <- cd_degs %>% na.omit() %>% dplyr::select(c(Ensembl, in_cd)) %>% unique()
uc_degs <- uc_degs %>% na.omit() %>% dplyr::select(c(Ensembl, in_uc)) %>% unique()

##### Join organoid and biopsy data #####

# Join onto the all_degs2 table based on the ensemblid
# Extract only the rows with data in  uc column (then crohns col)
all_degs_uc <- left_join(all_degs2, uc_degs, by = c("ENSEMBL.ID"="Ensembl")) %>% filter(in_uc == TRUE)
all_degs_cd <- left_join(all_degs2, cd_degs, by = c("ENSEMBL.ID"="Ensembl")) %>% filter(in_cd == TRUE)

# Save these lists of biopsy degs which are in the colonoid datasets
#write.table(all_degs_uc, "uc_biopsy_degs_in_colonoid_data_adjp_0.01_lfc_1.5.txt", sep = "\t", row.names = F, quote = F)
#write.table(all_degs_cd, "cd_biopsy_degs_in_colonoid_data_adjp_0.01_lfc_1.5.txt", sep = "\t", row.names = F, quote = F)
write.table(all_degs_uc, "uc_biopsy_degs_in_colonoid_data_adjp_0.01.txt", sep = "\t", row.names = F, quote = F)
write.table(all_degs_cd, "cd_biopsy_degs_in_colonoid_data_adjp_0.01.txt", sep = "\t", row.names = F, quote = F)


# Collapse the table by the cytokines column with count
uc_count <- all_degs_uc %>% group_by(cytokines) %>% summarize(uc_count=n())
cd_count <- all_degs_cd %>% group_by(cytokines) %>% summarize(cd_count=n())

# Join this resulting table to the degs_count table so 3 columns: cytokine, count in organoids, count in uc (or crohns). This will be the category, start count and end count.
# Add cytokine cat again in capitals
# Reorder columns
all_uc_count <- degs_count %>% left_join(uc_count, by = c("cytokines")) %>%
  mutate("CYTOKINE" = toupper(cytokines)) %>% dplyr::select(c("cytokines", "CYTOKINE", "org_count", "uc_count"))
# Replace na with 0
all_uc_count$uc_count <- as.double(all_uc_count$uc_count)
all_uc_count$uc_count <- all_uc_count$uc_count %>% replace_na(0.1)
all_uc_count$org_count <- as.double(all_uc_count$org_count)
all_uc_count$org_count <- all_uc_count$org_count %>% replace_na(0.1)

all_cd_count <- degs_count %>% left_join(cd_count, by = c("cytokines")) %>%
  mutate("CYTOKINE" = toupper(cytokines)) %>% dplyr::select(c("cytokines", "CYTOKINE", "org_count", "cd_count"))
all_cd_count$cd_count <- as.double(all_cd_count$cd_count)
all_cd_count$cd_count <- all_cd_count$cd_count %>% replace_na(0.1)
all_cd_count$org_count <- as.double(all_cd_count$org_count)
all_cd_count$org_count <- all_cd_count$org_count %>% replace_na(0.1)

# Filter so that both categories are more than 5
#all_uc_count_filt <- all_uc_count %>% filter(uc_count > 5 & org_count > 5)
#all_cd_count_filt <- all_cd_count %>% filter(cd_count > 5 & org_count > 5)

# Final format: cytokine_cat(source), Cytokine_cat_capitals(target), start_count, end_count
rm(cd_degs, uc_degs, all_degs_uc, all_degs_cd, uc_count, cd_count)

#### Plot data #####

# Set up colours properly
# Convert to ';' separators if necessary
colours <- colours %>% mutate(category = gsub(",","; ", category))
colours_v <- setNames(as.character(colours$colour), colours$category)

# Get required length of gaps so that the proportions are the same between the figures - need to make everything scale to 360
#uc_total <- sum(all_uc_count$org_count) + sum(all_uc_count$uc_count) # divide by 360 
#cd_total <- sum(all_cd_count$org_count) + sum(all_cd_count$cd_count)
all_uc_count <- all_uc_count %>% mutate(uc_count_scale = uc_count/1.5) %>% mutate(org_count_scale = org_count/1.5) %>% select(-c(uc_count, org_count))
all_cd_count <- all_cd_count %>% mutate(cd_count_scale = cd_count/1.5) %>% mutate(org_count_scale = org_count/1.5) %>% select(-c(cd_count, org_count))
uc_total <- sum(all_uc_count$org_count_scale) + sum(all_uc_count$uc_count_scale)
cd_total <- sum(all_cd_count$org_count_scale) + sum(all_cd_count$cd_count_scale)
uc_gap <- (360 - uc_total)/2 
cd_gap <- (360 - cd_total)/2

## All data

# UC
# Get list of categories
list_cyto_uc <- all_uc_count[['cytokines']]
list_cap_cyto_uc <- rev(all_uc_count[['CYTOKINE']])
len_uc <- length(list_cyto_uc) -1

circos.clear()
circos.par(start.degree = 307, gap.after=c(rep(0, len_uc), uc_gap, rep(0, len_uc), uc_gap))
chordDiagram(all_uc_count, grid.col=colours_v, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.07),
             annotationTrackHeight = c(0.07, 1), order = c(list_cyto_uc, list_cap_cyto_uc))


# Crohns
list_cyto_cd <- all_cd_count[['cytokines']]
list_cap_cyto_cd <- rev(all_cd_count[['CYTOKINE']])
len_cd <- length(list_cyto_cd) -1

circos.clear()
circos.par(start.degree = 301, gap.after=c(rep(0, len_cd), cd_gap, rep(0, len_cd), cd_gap))
chordDiagram(all_cd_count, grid.col=colours_v, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.07),
             annotationTrackHeight = c(0.07, 1), order = c(list_cyto_cd, list_cap_cyto_cd))

#for(si in get.all.sector.index()) {
#  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
#  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
#  circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 1, 
#              facing = "inside", niceFacing = TRUE, col = "black",cex=1.5)
#}
