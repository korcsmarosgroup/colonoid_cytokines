# Investigate the UC/CD biopsy DEGs which are also in the colonoid data. 
# What is the difference between UC and CD, and what are the functions of the cytokine specific biopsy DEGs?
#
# Input: 
#       uc/cd_degs_f - file path to text files of biopsy degs which are in the colonoid datasets, labelled by which cytokine they are targeted by in teh colonoid datasets.
#
# Output:
#       1. Tab delimited text file summarising how many genes are uc specific, cd specific or shared for each cytokine category (of the biopsy degs in the colonoid dataset)
#       2. Tab delimited text file with the gene IDs for each category
#       3. Venn diagrams of the above data - one for each cytokine category
#       4. Tab delimited text file summarising the significant Reactome pathways in each category (cytokine, uc/cd overlap)

##### Setup #####

library(dplyr)
library(VennDiagram)
library(ReactomePA)
library(clusterProfiler)
library(tidyr)
library(ggplot2)
library(scales)
library(ggthemes)

uc_degs_f <- "../input_data/differential_expression_datasets/biopsies/uc_biopsy_degs_in_colonoid_data_adjp_0.01.txt"
cd_degs_f <- "../input_data/differential_expression_datasets/biopsies/cd_biopsy_degs_in_colonoid_data_adjp_0.01.txt"

uc_degs <- read.csv(uc_degs_f, sep = "\t")
cd_degs <- read.csv(cd_degs_f, sep = "\t")

##### Count overlaps #####

# Join the datasets on every matching column
uc_cd_degs <- full_join(uc_degs, cd_degs)

# How many genes of each category of cytokines are shared in uc and cd?
# Add column to indicate if gene in both uc and cd
uc_cd_degs2 <- uc_cd_degs %>% mutate(uc_and_cd = ifelse(!is.na(in_uc) & !is.na(in_cd), 1, 0))
# Add column to indicate if gene in uc and not cd
uc_cd_degs2 <- uc_cd_degs2 %>% mutate(uc_not_cd = ifelse(in_uc == "TRUE" & is.na(in_cd), 1, 0))
# Add column to indicate if gene in cd and not uc
uc_cd_degs2 <- uc_cd_degs2 %>% mutate(cd_not_uc = ifelse(in_cd == "TRUE" & is.na(in_uc), 1, 0))
# Collapse by cytokine and count the number of genes in each category
uc_cd_degs_col <- uc_cd_degs2 %>% 
  select(-c(in_uc, in_cd, SYMBOL, ENSEMBL.ID)) %>%
  group_by(cytokines) %>%
  summarise_all(sum)

# Save the summary file showing the overlaps
write.table(uc_cd_degs_col, "../uc_and_cd_biopsy_degs_in_colonoid_data_adjp_0.01.txt", sep = "\t", quote= F, row.names = F)
# Save summary with the genes in the overlaps
write.table(uc_cd_degs, "../genes_uc_and_cd_biopsy_degs_in_colonoid_data_adjp_0.01.txt", sep = "\t", quote= F, row.names = F)

##### Plot overlaps #####

# Reformat data into long format
uc_cd_degs_col1 <- uc_cd_degs_col %>% select(cytokines,  "cCD not UC"= cd_not_uc, "cCD & UC" = uc_and_cd, "UC not cCD" =uc_not_cd)
uc_cd_degs_long <- gather(uc_cd_degs_col1, condition, number_of_genes, "cCD not UC":"UC not cCD", factor_key=TRUE)

# Plot as horizontal stacked bar chart
degs_plot <- ggplot(data = uc_cd_degs_long, aes(x = cytokines, y = number_of_genes, fill = condition)) + 
  geom_col() + coord_flip() + scale_fill_colorblind() + #geom_text(aes(label = number_of_genes)) +
  theme(text = element_text(size=20), axis.text = element_text(color = "black"), plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Cytokine group") + ylab("Number of genes") + labs(fill = "Condition") 

# Save plot
ggsave("../plot_uc_and_cd_biopsy_degs_in_colonoid_data_adjp_0.01.pdf", degs_plot, height = 5, width = 10, device = "pdf")

##### Plot Venns #####

# Iterate the cytokine categories
for (i in 1:nrow(uc_cd_degs_col)){
  
  uc_all <- uc_cd_degs_col[i,2] + uc_cd_degs_col[i,3]
  cd_all <- uc_cd_degs_col[i,2] + uc_cd_degs_col[i,4]
  cyto <- uc_cd_degs_col[[i,1]]
    
  pdf(paste0("../", cyto, "_biopsy_degs_in_colonoid_data_adjp_0.01.pdf"),width=7,height=7) 
  draw.pairwise.venn(area1 = uc_all, area2 = cd_all, area3 = uc_cd_degs_col[i,4], cross.area = uc_cd_degs_col[i,2], 
                   category = c("UC DEGs", "CD DEGs"), lty = "blank",
                   fontface = rep("plain", 3), fontfamily = rep("sans", 3), cex = rep(6, 3),
                   cat.fontface = rep("plain", 2), cat.fontfamily = rep("sans", 2), cat.cex = rep(2,2),
                   fill = c("skyblue", "orange"))
  dev.off()
}

##### Functional analysis #####

# Table for all the results
all_reactome <- data.frame()

# Iterate data in the uc_cd_degs_col df, where more than 5, carry out Reactome functional analysis
for (i in 1:nrow(uc_cd_degs_col)) {
  for (j in 2:ncol(uc_cd_degs_col)) {
    if (uc_cd_degs_col[i,j] >= 5){
      
      # Get cytokine name
      cyto <- uc_cd_degs_col[[i,1]]
      # Get category
      cat <- colnames(uc_cd_degs_col)[j]
      
      # Get gene list
      gene_ids <- uc_cd_degs2 %>% filter((UQ(as.symbol(cat)) == 1)&(cytokines == cyto))
      
      # Convert to Entrez IDs
      gene_ent <- bitr(gene_ids$ENSEMBL.ID, fromType = "ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
      
      # Reactome pathway enrichment
      genes_reactome <- enrichPathway(gene=gene_ent$ENTREZID,pvalueCutoff=0.1,readable=T)
      # Get reactome results
      genes_reactome2 <- as.data.frame(genes_reactome)
      genes_reactome2 <- mutate(genes_reactome2, Cytokine = cyto, Category = cat)
      
      # Append results to previous results
      all_reactome <- rbind(all_reactome, genes_reactome2)
    }
  } 
}

# Save results
write.table(all_reactome, "../reactome_uc_and_cd_biopsy_degs_in_colonoid_data_adjp_0.01.txt", sep = "\t", quote= F, row.names = F)
