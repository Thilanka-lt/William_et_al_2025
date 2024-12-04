#This is a code that calculates the differential expression between the two groups of interest # nolint
############################################################################################################ # nolint
rm(list=ls())
setwd("./")
############################################################################################################ # nolint
#librarirs#
library("DESeq2")
library("ggplot2")
library("dplyr")
#library(ComplexHeatmap)
#library(circlize)
####################
#function to calculate the differential expression#
###################################################
do_dif_exp <- function(counts, Desing, reference){
    "this is a function that calculates the differential expression between the two groups of interest
    Input: 1. counts: counts dara frame
           2. desing: desing data frame (should have single column specifying condions, sample names should be the rownames)
           3. reference: the reference level for the factor of interest (e.g. control group) "

  #make the DESeqDataSet object
  #
  #check if the column name is "parent"
    if(!"parent" %in% colnames(Desing)){
        stop("The column name of the desing table should be parent")
    }
    #check if the rownames of the desing table are the same as the column names of the counts table
    if(!all(rownames(Desing) %in% colnames(counts))){
        stop("The rownames of the desing table should be the same as the column names of the counts table")
    }

  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = Desing,
                                design = ~ parent)
  #prefiltering
  #The DESeq2 package provides a function to filter out low count rows of the count matrix,
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  #determining the factor level
  #By default, R will choose a reference level for factors based on alphabetical order.
  # Then, if you never tell the DESeq2 functions which level you want to compare against 
  #(e.g. which level represents the control group), the comparisons will be based on the 
  #alphabetical order of the levels
  
  #using relevel, just specifying the reference level:
  dds$parent <- relevel(dds$parent, ref = reference)
  #Diffrential expression analysis
  dds <- DESeq(dds)
  res <- results(dds)
  return(res)
}
############################################################################################################ # nolint
#read in the data#
########################
##Counts data#######
#########################
#DM1S1
counts_DM1S1 <- read.csv('/Users/thilanka_ranaweera/Documents/potato_colb/DGE/Read_counts_DM1S1_mapping.txt',head=T,stringsAsFactors=F, row.names = 1,sep = "\t")
head(counts_DM1S1)
counts_DM1S1["Soltu.DM1S1.01G015930",]
#replce X in the colnames with DM1S1
colnames(counts_DM1S1) <- gsub("X","DM1S1_",colnames(counts_DM1S1))
nrow(counts_DM1S1)

#VER
counts_VER <- read.csv('/Users/thilanka_ranaweera/Documents/potato_colb/DGE/Read_counts_VER_mapping.txt',head=T,stringsAsFactors=F, row.names = 1,sep = "\t")
colnames(counts_VER) <- gsub("_VER","",colnames(counts_VER))
#replce X in the colnames with VER
colnames(counts_VER) <- gsub("X","VER_",colnames(counts_VER))
head(counts_VER)
nrow(counts_VER)
counts_VER["Solver.v1.01_VERG042990",]

designer_table <- read.table('/Users/thilanka_ranaweera/Documents/potato_colb/DGE/Metadata_will_thilanka_edit.txt',head=T,stringsAsFactors=F, sep = "\t")

#name conversion file
nm_convo <- read.table('/Users/thilanka_ranaweera/Documents/potato_colb/DM1S1_vs_VER_blast_top_hit.txt',head=F,stringsAsFactors=F, sep = "\t")
head(nm_convo)
nrow(nm_convo)
length(nm_convo$V1)
length(unique(nm_convo$V1))
length(unique(nm_convo$V2))

############################################################################################################
############################################################################################################
############################################################################################################
#clean the data

#remove .1 from gene names
nm_convo$V1 <- gsub("\\.[^.]*$", "", nm_convo$V1)
nm_convo$V2 <- gsub("\\.[^.]*$", "", nm_convo$V2)

#remove duplicated raws in V1 keeping just on row of the duplicated rows
nm_convo <- nm_convo[!duplicated(nm_convo$V1),]
nm_convo_new <- data.frame(matrix(ncol = 3,nrow = 0))
for (dm1s1_name in unique(nm_convo$V2)){
  #dm1s1_name <- 'Soltu.DM1S1.12G024340'
  sub_name_conv = nm_convo[nm_convo[,"V2"] == dm1s1_name,]
  if(nrow(sub_name_conv) > 1){
    #chage V2 names by adding _1, _2, _3, ... to the duplicated names
    for(i in 1:nrow(sub_name_conv)){
      sub_name_conv[i,"V2"] <- paste0(sub_name_conv[i,"V2"],"_",i)
    }
    nm_convo_new <- rbind(nm_convo_new,sub_name_conv)
  }
  else{
    nm_convo_new <- rbind(nm_convo_new,sub_name_conv)
  }
}
nrow(nm_convo_new)
nrow(nm_convo)
#order nm_convo based on V2 column
nm_convo_new <- nm_convo_new[order(nm_convo_new$V2),]
#get the matching rows in the counts data
#match the gene names in the counts data with the gene names in the nm_convo
counts_VER_mod <- counts_VER[nm_convo_new$V1,]
nrow(counts_VER_mod)

#get the V2 genes that has a underscore in them
#nm_convo_underscore <- nm_convo[nm_convo$V2 %in% grep("_",nm_convo_new$V2,value = T),2]
counts_DM1S1_mod <- counts_DM1S1
for(name in rownames(counts_DM1S1_mod)){
  #name ="Soltu.DM1S1.00G000040"
  #check if nm_convo_underscore has values starting with the name
  name_match <- grep(paste0("^",name),nm_convo_new$V2,value = T)
  name_row <- counts_DM1S1_mod[name,] 
  #if nor length of name_match not equal to 0
  if(length(name_match) > 1){
      print(name)
    #empty dataframe with rownames as the names in name_match
    sub_df <- data.frame(matrix(ncol = ncol(name_row),nrow = length(name_match)))
    colnames(sub_df) <- colnames(counts_DM1S1_mod)
    rownames(sub_df) <- name_match
    #fill each row of the with name_row[1,]
    for(i in 1:length(name_match)){
      sub_df[i,] <- name_row[1,]
      }

    #remove name_row
    counts_DM1S1_mod <- counts_DM1S1_mod[!rownames(counts_DM1S1_mod) %in% name,]
    #add the duplicated rows to the counts_DM1S1
    counts_DM1S1_mod <- rbind(counts_DM1S1_mod,sub_df)
  }

    

}
#look at 35357th row in the counts_DM1S1_mod
counts_DM1S1_mod["Soltu.DM1S1.01G015930",]

counts_DM1S1_mod_2 <- counts_DM1S1_mod
nrow(counts_DM1S1_mod_2)
#counts_DM1S1_mod_2 <- counts_DM1S1_mod_2[nm_convo$V2,]
#get the rownames that does not match with nm_convo$V2
#rownames(counts_DM1S1_mod_2[!rownames(counts_DM1S1_mod_2) %in% nm_convo_new$V2,])

counts_DM1S1_mod_2 <- counts_DM1S1_mod_2[nm_convo_new$V2,]
nrow(counts_DM1S1_mod_2)
nrow(counts_VER_mod)






rownames(counts_DM1S1_mod_2) <- nm_convo_new$V1[match(rownames(counts_DM1S1_mod_2), nm_convo_new$V2)]

counts_DM1S1_mod_2["Solver.v1.01_VERG042990",]


head(counts_DM1S1_mod_2)
head(counts_VER_mod)
nrow(counts_DM1S1_mod_2)
nrow(counts_VER_mod)
rownames(counts_DM1S1_mod_2)
rownames(counts_VER_mod)

#check if rownames of the two dataframes are the same
length(intersect(rownames(counts_DM1S1_mod_2),rownames(counts_VER_mod)))
#join two dataframes using join
counts_DM1S1_mod_2 <- tibble::rownames_to_column(counts_DM1S1_mod_2, "rowname")
counts_VER_mod <- tibble::rownames_to_column(counts_VER_mod, "rowname")
head(counts_DM1S1_mod_2)
head(counts_VER_mod)

counts_all <- left_join(counts_DM1S1_mod_2, counts_VER_mod, by = "rowname")
head(counts_all)
#make rowname the rownames of the dataframe
rownames(counts_all) <- counts_all$rowname
counts_all <- counts_all[,-1]
counts_all["Solver.v1.01_VERG042990",]
#double checking counts all for few genes
counts_DM1S1["Soltu.DM1S1.01G015800",]
counts_VER["Solver.v1.01_VERG042630",]
counts_all["Solver.v1.01_VERG042630",]
#save names conversion file
write.csv(nm_convo_new, file = "/Users/thilanka_ranaweera/Documents/potato_colb/DGE/gene_id_conversion_between_DM1S1_and_VER.csv",quote = F,row.names = F)
############################################################################################################
############################################################################################################
############################################################################################################

#calculating differential expression between F2_S2 vs F2_S1
designer_table_S2_S1 <- designer_table[designer_table[,2]== "F2_S2" | designer_table[,2]== "F2_S1" ,]
rownames(designer_table_S2_S1) <- designer_table_S2_S1[,1]
head(counts_DM1S1_mod_2)
head(counts_VER_mod)
#get counts for F2_S2 and F2_S1
counts_S2_S1 <- counts_all[,designer_table_S2_S1$experiment]
head(counts_S2_S1)


#calculate the differential expression
res_S2_S1 <- do_dif_exp(counts_S2_S1, designer_table_S2_S1, "F2_S1")
#save the differential expression results
write.csv(res_S2_S1, file = "/Users/thilanka_ranaweera/Documents/potato_colb/DGE/Fold_change_values/differential_expression_between_F2_S2_S1.csv",quote = F)

#calculating differential expression between F2_S2 vs F2_S3
designer_table_S2_S3 <- designer_table[designer_table[,2]== "F2_S2" | designer_table[,2]== "F2_S3" ,]
rownames(designer_table_S2_S3) <- designer_table_S2_S3[,1]
head(counts_DM1S1_mod_2)
head(counts_VER_mod)
#get counts for F2_S2 and F2_S3
counts_S2_S3 <- counts_all[,designer_table_S2_S3$experiment]
head(counts_S2_S3)

#calculate the differential expression
res_S2_S3 <- do_dif_exp(counts_S2_S3, designer_table_S2_S3, "F2_S2")
res_S2_S3["Solver.v1.01_VERG042990",]

#save the differential expression results
write.csv(res_S2_S3, file = "/Users/thilanka_ranaweera/Documents/potato_colb/DGE/Fold_change_values/differential_expression_between_F2_S2_S3.csv",quote = F)



#calculating differential expression between F2_S1 vs P1_DM1S1
designer_table_S1_DM1S1 <- designer_table[designer_table[,2]== "F2_S1" | designer_table[,2]== "P1_DM1S1" ,]
rownames(designer_table_S1_DM1S1) <- designer_table_S1_DM1S1[,1]

#get counts for F2_S1 and P1_DM1S1
counts_S1_DM1S1 <- counts_all[,designer_table_S1_DM1S1$experiment]

#calculate the differential expression
res_S1_DM1S1 <- do_dif_exp(counts_S1_DM1S1, designer_table_S1_DM1S1, "P1_DM1S1")
res_S1_DM1S1["Solver.v1.01_VERG042990",]
#save the differential expression results
write.csv(res_S1_DM1S1, file = "/Users/thilanka_ranaweera/Documents/potato_colb/DGE/Fold_change_values/differential_expression_between_F2_S1_P1_DM1S1.csv",quote = F)


#calculating differential expression between F2_S2 vs P1_DM1S1
designer_table_S2_DM1S1 <- designer_table[designer_table[,2]== "F2_S2" | designer_table[,2]== "P1_DM1S1" ,]
rownames(designer_table_S2_DM1S1) <- designer_table_S2_DM1S1[,1]

#get counts for F2_S2 and P1_DM1S1
counts_S2_DM1S1 <- counts_all[,designer_table_S2_DM1S1$experiment]

#calculate the differential expression
res_S2_DM1S1 <- do_dif_exp(counts_S2_DM1S1, designer_table_S2_DM1S1, "P1_DM1S1")
res_S2_DM1S1["Solver.v1.01_VERG042990",]
#save the differential expression results
write.csv(res_S2_DM1S1, file = "/Users/thilanka_ranaweera/Documents/potato_colb/DGE/Fold_change_values/differential_expression_between_F2_S2_P1_DM1S1.csv",quote = F)


#calculating differential expression between F2_S3 vs P1_DM1S1
designer_table_S3_DM1S1 <- designer_table[designer_table[,2]== "F2_S3" | designer_table[,2]== "P1_DM1S1" ,]

rownames(designer_table_S3_DM1S1) <- designer_table_S3_DM1S1[,1]

#get counts for F2_S3 and P1_DM1S1
counts_S3_DM1S1 <- counts_all[,designer_table_S3_DM1S1$experiment]

#calculate the differential expression
res_S3_DM1S1 <- do_dif_exp(counts_S3_DM1S1, designer_table_S3_DM1S1, "P1_DM1S1")
res_S3_DM1S1["Solver.v1.01_VERG042990",]
#save the differential expression results
write.csv(res_S3_DM1S1, file = "/Users/thilanka_ranaweera/Documents/potato_colb/DGE/Fold_change_values/differential_expression_between_F2_S3_P1_DM1S1.csv",quote = F)


#calculating differential expression between F2_S1 vs P2_VER
designer_table_S1_VER <- designer_table[designer_table[,2]== "F2_S1" | designer_table[,2]== "P2_VER" ,]
rownames(designer_table_S1_VER) <- designer_table_S1_VER[,1]
#get counts for F2_S1 and P2_VER
counts_S1_VER <- counts_all[,designer_table_S1_VER$experiment]
#calculate the differential expression
res_S1_VER <- do_dif_exp(counts_S1_VER, designer_table_S1_VER, "P2_VER")
res_S1_VER["Solver.v1.01_VERG042990",]
counts_S1_VER["Solver.v1.01_VERG042990",]
#save the differential expression results
write.csv(res_S1_VER, file = "/Users/thilanka_ranaweera/Documents/potato_colb/DGE/Fold_change_values/differential_expression_between_F2_S1_P2_VER.csv",quote = F)

#calculating differential expression between F2_S2 vs P2_VER
designer_table_S2_VER <- designer_table[designer_table[,2]== "F2_S2" | designer_table[,2]== "P2_VER" ,]
rownames(designer_table_S2_VER) <- designer_table_S2_VER[,1]
#get counts for F2_S2 and P2_VER
counts_S2_VER <- counts_all[,designer_table_S2_VER$experiment]
#calculate the differential expression
res_S2_VER <- do_dif_exp(counts_S2_VER, designer_table_S2_VER, "P2_VER")
res_S2_VER["Solver.v1.01_VERG042990",]
#save the differential expression results
write.csv(res_S2_VER, file = "/Users/thilanka_ranaweera/Documents/potato_colb/DGE/Fold_change_values/differential_expression_between_F2_S2_P2_VER.csv",quote = F)

#calculating differential expression between F2_S3 vs P2_VER
designer_table_S3_VER <- designer_table[designer_table[,2]== "F2_S3" | designer_table[,2]== "P2_VER" ,]
rownames(designer_table_S3_VER) <- designer_table_S3_VER[,1]
#get counts for F2_S3 and P2_VER
counts_S3_VER <- counts_all[,designer_table_S3_VER$experiment]
#calculate the differential expression
res_S3_VER <- do_dif_exp(counts_S3_VER, designer_table_S3_VER, "P2_VER")
res_S3_VER["Solver.v1.01_VERG042990",]
#save the differential expression results
write.csv(res_S3_VER, file = "/Users/thilanka_ranaweera/Documents/potato_colb/DGE/Fold_change_values/differential_expression_between_F2_S3_P2_VER.csv",quote = F)
############################################################################################################
############################################################################################################
############################################################################################################
#Visualizing the differential expression results


res_S2_DM1S1 <- read.csv("/Users/thilanka_ranaweera/Documents/potato_colb/DGE/Fold_change_values/differential_expression_between_F2_S2_P1_DM1S1.csv",head=T,stringsAsFactors=F, row.names = 1)
res_S1_DM1S1 <- read.csv("/Users/thilanka_ranaweera/Documents/potato_colb/DGE/Fold_change_values/differential_expression_between_F2_S1_P1_DM1S1.csv",head=T,stringsAsFactors=F, row.names = 1)
nrow(res_S2_DM1S1)
nrow(res_S1_DM1S1)

#get significant genes
res_S2_DM1S1_sig <- res_S2_DM1S1[res_S2_DM1S1$padj < 0.05,2, drop = F]
res_S1_DM1S1_sig <- res_S1_DM1S1[res_S1_DM1S1$padj < 0.05,2, drop = F]
nrow(res_S2_DM1S1_sig)
nrow(res_S1_DM1S1_sig)
#get commonn significant genes
common_sig_genes <- intersect(rownames(res_S2_DM1S1_sig),rownames(res_S1_DM1S1_sig))
#make rownames a column

res_S2_DM1S1_sig_com <- res_S2_DM1S1_sig[common_sig_genes, ,drop = F]
res_S1_DM1S1_sig_com <- res_S1_DM1S1_sig[common_sig_genes, , drop = F]

res_S2_DM1S1_sig_com <- tibble::rownames_to_column(as.data.frame(res_S2_DM1S1_sig_com), "rowname")
res_S1_DM1S1_sig_com <- tibble::rownames_to_column(as.data.frame(res_S1_DM1S1_sig_com), "rowname")

#ggs scatter plotbetween log2 fold change of  res_S2_DM1S1 and res_S1_DM1S1

flod_change_comparisons <- left_join(res_S2_DM1S1_sig_com, res_S1_DM1S1_sig_com, by = "rowname")
#make the rowname the rownames of the dataframe
rownames(flod_change_comparisons) <- flod_change_comparisons$rowname
flod_change_comparisons <- flod_change_comparisons[,-1]
head(flod_change_comparisons)
colnames(flod_change_comparisons) <- c("log2FoldChange_S2_DM1S1","log2FoldChange_S1_DM1S1")
#plot the log2 fold change of res_S2_DM1S1 and res_S1_DM1S1
p1 <- ggplot(flod_change_comparisons, aes(x = log2FoldChange_S2_DM1S1, y = log2FoldChange_S1_DM1S1)) +
  #set alpha to 0.5
  geom_point(alpha = 0.5) +
  #x axis and y axis limits
  xlim(c(-20,20)) + ylim(c(-20,20)) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  #abline on x axis
  geom_hline(yintercept = -1, color = "blue") + geom_hline(yintercept = 1, color = "blue") +
  #abline on y axis
  geom_vline(xintercept = -1, color = "blue") + geom_vline(xintercept = 1, color = "blue") +
  xlab("log2FoldChange_S2_DM1S1") +
  ylab("log2FoldChange_S1_DM1S1") +
  ggtitle("log2FoldChange_S2_DM1S1 vs log2FoldChange_S1_DM1S1")


#dothe same with res_S2_DM1S1 res_S3_DM1S1
res_S3_DM1S1 <- read.csv("/Users/thilanka_ranaweera/Documents/potato_colb/DGE/Fold_change_values/differential_expression_between_F2_S3_P1_DM1S1.csv",head=T,stringsAsFactors=F, row.names = 1)
nrow(res_S3_DM1S1)
#get significant genes
res_S3_DM1S1_sig <- res_S3_DM1S1[res_S3_DM1S1$padj < 0.05,2, drop = F]
nrow(res_S3_DM1S1_sig)
#get commonn significant genes
common_sig_genes <- intersect(rownames(res_S2_DM1S1_sig),rownames(res_S3_DM1S1_sig))
#make rownames a column

res_S3_DM1S1_sig_com <- res_S3_DM1S1_sig[common_sig_genes, ,drop = F]

res_S3_DM1S1_sig_com <- tibble::rownames_to_column(as.data.frame(res_S3_DM1S1_sig_com), "rowname")

#ggs scatter plotbetween log2 fold change of  res_S2_DM1S1 and res_S3_DM1S1

flod_change_comparisons <- left_join(res_S2_DM1S1_sig_com, res_S3_DM1S1_sig_com, by = "rowname")
#make the rowname the rownames of the dataframe
rownames(flod_change_comparisons) <- flod_change_comparisons$rowname
flod_change_comparisons <- flod_change_comparisons[,-1]
head(flod_change_comparisons)
colnames(flod_change_comparisons) <- c("log2FoldChange_S2_DM1S1","log2FoldChange_S3_DM1S1")
#plot the log2 fold change of res_S2_DM1S1 and res_S3_DM1S1

p2 <- ggplot(flod_change_comparisons, aes(x = log2FoldChange_S2_DM1S1, y = log2FoldChange_S3_DM1S1)) +
  #set alpha to 0.5
  geom_point(alpha = 0.5) +
  #x axis and y axis limits
  xlim(c(-20,20)) + ylim(c(-20,20)) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  #abline on x axis
  geom_hline(yintercept = -1, color = "blue") + geom_hline(yintercept = 1, color = "blue") +
  #abline on y axis
  geom_vline(xintercept = -1, color = "blue") + geom_vline(xintercept = 1, color = "blue") +
  xlab("log2FoldChange_S2_DM1S1") +
  ylab("log2FoldChange_S3_DM1S1") +
  ggtitle("log2FoldChange_S2_DM1S1 vs log2FoldChange_S3_DM1S1")


#draw beween res_S1_DM1S1 and res_S3_DM1S1

flod_change_comparisons <- left_join(res_S1_DM1S1_sig_com, res_S3_DM1S1_sig_com, by = "rowname")
#make the rowname the rownames of the dataframe
rownames(flod_change_comparisons) <- flod_change_comparisons$rowname
flod_change_comparisons <- flod_change_comparisons[,-1]
head(flod_change_comparisons)
colnames(flod_change_comparisons) <- c("log2FoldChange_S1_DM1S1","log2FoldChange_S3_DM1S1")
#plot the log2 fold change of res_S1_DM1S1 and res_S3_DM1S1

p3 <- ggplot(flod_change_comparisons, aes(x = log2FoldChange_S1_DM1S1, y = log2FoldChange_S3_DM1S1)) +
  #set alpha to 0.5
  geom_point(alpha = 0.5) +
  #x axis and y axis limits
  xlim(c(-20,20)) + ylim(c(-20,20)) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  #abline on x axis
  geom_hline(yintercept = -1, color = "blue") + geom_hline(yintercept = 1, color = "blue") +
  #abline on y axis
  geom_vline(xintercept = -1, color = "blue") + geom_vline(xintercept = 1, color = "blue") +
  xlab("log2FoldChange_S1_DM1S1") +
  ylab("log2FoldChange_S3_DM1S1") +
  ggtitle("log2FoldChange_S1_DM1S1 vs log2FoldChange_S3_DM1S1")




#draw all three plots in one pdf using grid.arrange
library(gridExtra)
pdf("/Users/thilanka_ranaweera/Documents/potato_colb/DGE/Figures/fold_change_comparisons_between_F2s_vs_DM1S1_parents.pdf", width = 28, height = 10)
grid.arrange(p1,p2,p3, ncol = 3)
dev.off()
