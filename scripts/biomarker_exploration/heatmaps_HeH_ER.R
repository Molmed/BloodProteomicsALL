#setwd("~/Phd_Uppsala_University/2021/Proteomics_Projects/MariaG/Olink/Scripts/")

library("RColorBrewer")
library("pheatmap")
library('ComplexHeatmap')
library(circlize)
x <- read.csv('outpath/Filtered_gex_all_patients.csv',  header= TRUE, check.names = FALSE)
pheno <- read.csv("outpath/Clinical_all_data_n=352_AL.csv", row.names = 1, header= TRUE, check.names = FALSE)

sum(colnames(x)[2:length(colnames(x))] == rownames(pheno))
# Keep only HeH and ETV6::RUNX1 from BCP-ALL, Keep Controls of course

pheno <- pheno[
  pheno$Final_Subtype %in% c("HeH", "t(12;21)") |
    pheno$group %in% c("T-ALL", "AML", "Control", "T controls", "B controls"),
]

dim(pheno)

table(pheno$group)


x2 <- cbind(
  x[, 1, drop = FALSE],
  x[, rownames(pheno), drop = FALSE]
)

dim(x2) # the extra column is the gene name

sum(colnames(x2)[2:length(colnames(x2))] == rownames(pheno))



proteins <- read.csv("outpath/protein_info.csv", header= TRUE, check.names = FALSE)

pheno$group2 <- pheno$group

names_dict = c('Control' = 'CD34+',
               'B controls' = 'CD19+',
               'T controls' = 'CD3+',
               'AML'= 'AML',
               'BCP-ALL' = 'BCP-ALL',
               'T-ALL' = 'T-ALL')

pheno$group2 <- names_dict[pheno$group2]
pheno$group2[pheno$Final_Subtype %in% 'HeH'] <- 'HeH'
pheno$group2[pheno$Final_Subtype %in% 't(12;21)'] <- 'ER'
table(pheno$group)
table(pheno$group2)
sum(colnames(x2)[2:length(colnames(x2))] == rownames(pheno))
sum(x2$gene == proteins$gene)
dim(x2)
#remove the gene name column
dim(x2[2:length(colnames(x2))])
# scale by the features - which in our case are the rows
df <- t(scale(t(x2[2:length(colnames(x2))])))
#df <- scale(x2[2:length(colnames(x2))])
dim(df)
min(df)
max(df)


cmap <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

p1 <- quantile(df, probs = 0.01)  # 1st percentile
p99 <- quantile(df, probs = 0.99)  # 99th percentile

# Define the custom colormap using percentiles
cmap <- colorRamp2(
  c(p1, 0, p99),                                # Map colors to percentiles
  colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(3)  # Colors
)

cmap <- colorRamp2(
  seq(-4, 4, length.out = 100),  # Map to -4 to 4
  colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
)


# Define colors for qualitative and continuous variables
comp_col = list(
  comparison = c('Leukemia vs controls' = '#2a9d8f',
  'AML vs controls' = '#e9c46a',
  'BCP-ALL vs controls'='#f4a261',
  'T-ALL vs controls'='#e76f51')
 
)

# Create the row annotation
la <- rowAnnotation(
  comparison = proteins$comparison,
  col = comp_col
)


immune_col <- list (immunophenotype = c(
  'AML' = '#e9c46a',
  'HeH'='#f4a261',
  'ER'='#b2560D',
  'T-ALL'='#e76f51',
  'CD34+' = '#f4cae4',
  'CD19+' = 'grey',
  'CD3+' = 'lightblue'
))


# Create column annotations for all column names
ta <- HeatmapAnnotation(
  immunophenotype = pheno$group2,  # Column names as annotations
  col = immune_col              # Use generated colors
)

#pdf(file = "../Figures/Supp_FigureS5/HeatMap_Transcriptomics_BCP_ALL_ER_HeH_only.pdf", width = 10, height = 10)
# Create the heatmap with row splitting and no clustering
Heatmap(
  as.matrix(df), # Select numerical columns for heatmap
  name = "Log2 GEX",
  row_split = proteins$comparison,
  row_labels = proteins$gene,

  #cluster_rows = FALSE,               # Disable row clustering
  row_names_gp = gpar(fontsize = 8),  # Adjust font size
  left_annotation = la,
  top_annotation = ta,# Add annotation on the left
  show_column_names = FALSE,
  
  
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  row_title = NULL,border = TRUE,column_title = NULL
  
  
)
#dev.off()
pdf(file = "HeatMap_Transcriptomics_immuno_groups_BCP_ALL_ER_HeH_only.pdf", width = 10, height = 10)
print(
Heatmap(
  as.matrix(df), # Select numerical columns for heatmap
  name = "Log2 GEX",
  row_split = proteins$comparison,
  column_split = pheno$group2,
  row_labels = proteins$gene,
  #cluster_rows = FALSE,               # Disable row clustering
  row_names_gp = gpar(fontsize = 8),  # Adjust font size
  left_annotation = la,
  top_annotation = ta,# Add annotation on the left
  show_column_names = FALSE,

  show_column_dend = FALSE,
  show_row_dend = FALSE,
  row_title = NULL,border = TRUE,column_title = NULL
  
)
)

dev.off()

