#This is a script to analyze and make figures for the differental gene expression analysis
#remove all data
rm(list = ls())
####_______________________________________________________________________________________________________________________####
# Load the necessary libraries
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tibble)
#_______________________________________________________________________________________________________________________

# Load the data
#gene_list

ch_1_qtl <- read.csv("/Users/thilanka_ranaweera/Documents/potato_colb/DGE/chr1QTL_genes_VER_v1.csv", header = TRUE, sep = ",")
ch_11_qtl <- read.csv("/Users/thilanka_ranaweera/Documents/potato_colb/DGE/chr11QTL_genes_VER_v1.csv", header = TRUE, sep = ",")

#differential_expression
#DM1S1 comparisons
s1_p1 <- read.csv("/Users/thilanka_ranaweera/Documents/potato_colb/DGE/Fold_change_values/differential_expression_between_F2_S1_P1_DM1S1.csv", header = TRUE, sep = ",", row.names = 1)
s2_p1 <- read.csv("/Users/thilanka_ranaweera/Documents/potato_colb/DGE/Fold_change_values/differential_expression_between_F2_S2_P1_DM1S1.csv", header = TRUE, sep = ",", row.names = 1)
s3_p1 <- read.csv("/Users/thilanka_ranaweera/Documents/potato_colb/DGE/Fold_change_values/differential_expression_between_F2_S3_P1_DM1S1.csv", header = TRUE, sep = ",", row.names = 1)
#DM1S2 comparisons
s1_p2 <- read.csv("/Users/thilanka_ranaweera/Documents/potato_colb/DGE/Fold_change_values/differential_expression_between_F2_S1_P2_VER.csv", header = TRUE, sep = ",", row.names = 1)
s2_p2 <- read.csv("/Users/thilanka_ranaweera/Documents/potato_colb/DGE/Fold_change_values/differential_expression_between_F2_S2_P2_VER.csv", header = TRUE, sep = ",", row.names = 1)
s3_p3 <- read.csv("/Users/thilanka_ranaweera/Documents/potato_colb/DGE/Fold_change_values/differential_expression_between_F2_S3_P2_VER.csv", header = TRUE, sep = ",", row.names = 1)

head(s1_p1)
#_______________________________________________________________________________________________________________________
#wrangle data

#get the diiferentlly expressed genes for ch1_qtl genes

#empty dataframe with ch1_qtl genes as row names
ch1_qtl_diff_ex <- data.frame(matrix(NA, nrow = nrow(ch_1_qtl),ncol = 1))
rownames(ch1_qtl_diff_ex) <- ch_1_qtl$ID

#fill the dataframe with the differentially expressed genes
#left join the ch1_qtl genes with the significantly differentially expressed genes
ch1_qtl_diff_ex <- left_join(rownames_to_column(ch1_qtl_diff_ex), rownames_to_column(s1_p1[s1_p1[,6] < 0.05,2,drop = F]), by = "rowname")
ch1_qtl_diff_ex <- ch1_qtl_diff_ex[,-2]
rownames(ch1_qtl_diff_ex) <- ch1_qtl_diff_ex$rowname
#remove the rowname column
ch1_qtl_diff_ex <- ch1_qtl_diff_ex[,-1,drop = F]
colnames(ch1_qtl_diff_ex) <- "S1_DM1S1"
head(ch1_qtl_diff_ex)

#joining s1_p1
s2_p1 <- rownames_to_column(s2_p1[s2_p1[,6] < 0.05,2,drop = F])
colnames(s2_p1) <- c("rowname","S2_DM1S1")
ch1_qtl_diff_ex <- left_join(rownames_to_column(ch1_qtl_diff_ex), s2_p1, by = "rowname")
head(ch1_qtl_diff_ex)

#joining s3_p1
s3_p1 <- rownames_to_column(s3_p1[s3_p1[,6] < 0.05,2,drop = F])
colnames(s3_p1) <- c("rowname","S3_DM1S1")
ch1_qtl_diff_ex <- left_join(ch1_qtl_diff_ex, s3_p1, by = "rowname")
head(ch1_qtl_diff_ex)

#joining s1_p2
s1_p2 <- rownames_to_column(s1_p2[s1_p2[,6] < 0.05,2,drop = F])
colnames(s1_p2) <- c("rowname","S1_VER")
ch1_qtl_diff_ex <- left_join(ch1_qtl_diff_ex, s1_p2, by = "rowname")
head(ch1_qtl_diff_ex)

#joining s2_p2
s2_p2 <- rownames_to_column(s2_p2[s2_p2[,6] < 0.05,2,drop = F])
colnames(s2_p2) <- c("rowname","S2_VER")
ch1_qtl_diff_ex <- left_join(ch1_qtl_diff_ex, s2_p2, by = "rowname")
head(ch1_qtl_diff_ex)

#joining s3_p3
s3_p3 <- rownames_to_column(s3_p3[s3_p3[,6] < 0.05,2,drop = F])
colnames(s3_p3) <- c("rowname","S3_VER")
ch1_qtl_diff_ex <- left_join(ch1_qtl_diff_ex, s3_p3, by = "rowname")

rownames(ch1_qtl_diff_ex) <- ch1_qtl_diff_ex$rowname
#remove the rowname column
ch1_qtl_diff_ex <- ch1_qtl_diff_ex[,-1,drop = F]

#drop the rows that have all NA values
ch1_qtl_diff_ex <- ch1_qtl_diff_ex[rowSums(is.na(ch1_qtl_diff_ex)) < ncol(ch1_qtl_diff_ex),]
nrow(ch1_qtl_diff_ex)

#save the data
write.csv(ch1_qtl_diff_ex, "/Users/thilanka_ranaweera/Documents/potato_colb/DGE/Fold_change_values/ch1_qtl_diff_ex.csv")

#_______________________________________________________________________________________________________________________
#set random seed
set.seed(25101991)
#prepare the data for the heatmap
ch1_qtl_diff_ex_mod <- ch1_qtl_diff_ex
ch1_qtl_diff_ex_mod[is.na(ch1_qtl_diff_ex_mod)] <- 0

col_fun = colorRamp2(c(-5, -2, 0, 2, 5), c("deepskyblue3","skyblue2","white","hotpink1", "deeppink2"))
col_fun(seq(-2, 2))


HM <- Heatmap(as.matrix(ch1_qtl_diff_ex_mod), name = "log2FC", col = col_fun, row_km = 5,cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 5), show_row_names = TRUE, show_row_dend = FALSE,
              row_gap = unit(3.5, "mm"),column_gap = unit(3.5, "mm"), column_names_gp = gpar(fontsize = 8), border = TRUE, 
              row_title_rot = 0, column_title_rot = 90, column_title_gp = gpar(fontsize = 20), row_title_gp = gpar(fontsize = 10))

pdf("/Users/thilanka_ranaweera/Documents/potato_colb/DGE/Figures/Significantly_differentilly_expressed_genes_in_ch1_QTL.pdf",width=10,height=30)
HM
dev.off()


#####################################################################################################################################################################
#####################################################################################################################################################################
#wrangle data for ch11_qtl genes

#get the diiferentlly expressed genes for ch11_qtl genes

#empty dataframe with ch11_qtl genes as row names
ch11_qtl_diff_ex <- data.frame(matrix(NA, nrow = nrow(ch_11_qtl),ncol = 1))
rownames(ch11_qtl_diff_ex) <- ch_11_qtl$ID

#fill the dataframe with the differentially expressed genes
#left join the ch11_qtl genes with the significantly differentially expressed genes
ch11_qtl_diff_ex <- left_join(rownames_to_column(ch11_qtl_diff_ex), rownames_to_column(s1_p1[s1_p1[,6] < 0.05,2,drop = F]), by = "rowname")
ch11_qtl_diff_ex <- ch11_qtl_diff_ex[,-2]
rownames(ch11_qtl_diff_ex) <- ch11_qtl_diff_ex$rowname
#remove the rowname column
ch11_qtl_diff_ex <- ch11_qtl_diff_ex[,-1,drop = F]
colnames(ch11_qtl_diff_ex) <- "S1_DM1S1"
head(ch11_qtl_diff_ex)

#joining s2_p1
s2_p1 <- rownames_to_column(s2_p1[s2_p1[,6] < 0.05,2,drop = F])
colnames(s2_p1) <- c("rowname","S2_DM1S1")
ch11_qtl_diff_ex <- left_join(rownames_to_column(ch11_qtl_diff_ex), s2_p1, by = "rowname")
head(ch11_qtl_diff_ex)

#joining s3_p1
s3_p1 <- rownames_to_column(s3_p1[s3_p1[,6] < 0.05,2,drop = F])
colnames(s3_p1) <- c("rowname","S3_DM1S1")
ch11_qtl_diff_ex <- left_join(ch11_qtl_diff_ex, s3_p1, by = "rowname")
head(ch11_qtl_diff_ex)

#joining s1_p2
s1_p2 <- rownames_to_column(s1_p2[s1_p2[,6] < 0.05,2,drop = F])
colnames(s1_p2) <- c("rowname","S1_VER")
ch11_qtl_diff_ex <- left_join(ch11_qtl_diff_ex, s1_p2, by = "rowname")
head(ch11_qtl_diff_ex)

#joining s2_p2
s2_p2 <- rownames_to_column(s2_p2[s2_p2[,6] < 0.05,2,drop = F])
colnames(s2_p2) <- c("rowname","S2_VER")
ch11_qtl_diff_ex <- left_join(ch11_qtl_diff_ex, s2_p2, by = "rowname")
head(ch11_qtl_diff_ex)

#joining s3_p3
s3_p3 <- rownames_to_column(s3_p3[s3_p3[,6] < 0.05,2,drop = F])
colnames(s3_p3) <- c("rowname","S3_VER")
ch11_qtl_diff_ex <- left_join(ch11_qtl_diff_ex, s3_p3, by = "rowname")

rownames(ch11_qtl_diff_ex) <- ch11_qtl_diff_ex$rowname
#remove the rowname column
ch11_qtl_diff_ex <- ch11_qtl_diff_ex[,-1,drop = F]

#drop the rows that have all NA values
ch11_qtl_diff_ex <- ch11_qtl_diff_ex[rowSums(is.na(ch11_qtl_diff_ex)) < ncol(ch11_qtl_diff_ex),]
nrow(ch11_qtl_diff_ex)
head(ch11_qtl_diff_ex)

#save the data
write.csv(ch11_qtl_diff_ex, "/Users/thilanka_ranaweera/Documents/potato_colb/DGE/Fold_change_values/ch11_qtl_diff_ex.csv")

#_______________________________________________________________________________________________________________________
#set random seed
set.seed(25101991)
#prepare the data for the heatmap
ch11_qtl_diff_ex_mod <- ch11_qtl_diff_ex
ch11_qtl_diff_ex_mod[is.na(ch11_qtl_diff_ex_mod)] <- 0

HM_11 <- Heatmap(as.matrix(ch11_qtl_diff_ex_mod), name = "log2FC", col = col_fun, row_km = 5,cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 5), show_row_names = TRUE, show_row_dend = FALSE,
              row_gap = unit(3.5, "mm"),column_gap = unit(3.5, "mm"), column_names_gp = gpar(fontsize = 8), border = TRUE, 
              row_title_rot = 0, column_title_rot = 90, column_title_gp = gpar(fontsize = 20), row_title_gp = gpar(fontsize = 10))

pdf("/Users/thilanka_ranaweera/Documents/potato_colb/DGE/Figures/Significantly_differentilly_expressed_genes_in_ch11_QTL.pdf",width=10,height=30)
HM_11
dev.off()
