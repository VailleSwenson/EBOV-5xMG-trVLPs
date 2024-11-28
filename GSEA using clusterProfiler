##############################################################################################
### Gene set enrichment analysis using ClusterProfiler ###
##############################################################################################

#Based on the tutorial by BioStatsSquid: https://biostatsquid.com/pathway-enrichment-analysis-tutorial-clusterprofile

#Begin by cleaning up the RStudio environment. 
rm(list = ls(all.names = TRUE))

#Load the libraries. 
library(ggplot2)          #Version 3.5.1
library(pheatmap)         #Version 1.0.12
library(clusterProfiler)  #Version 4.10.1
library(enrichplot)       #Version 1.22.0
library('org.Mm.eg.db')   #Version 3.18.0

#Adjacency matrix to list function (geneset file to a form readable by clusterProfiler).
#The elements of the new list become the column names from the matrix.
#The values of the list entry (element) are the entries from that column. 
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)){
    pws.1[[pws]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.1)
}

#Set global variable of what genes are in the dataset based on WT DEG list: 
GenesinData <- read.csv("Brady_DEG_WTvsMock.csv")
GenesinData$gene <- GenesinData$Gene.name
GenesinData <- subset(GenesinData, select = c(gene))
GenesinData <- na.omit(GenesinData)
bg_genesIDs <- GenesinData$gene

#Read in the pathway file as a global variable: 
#Read in the selected .gmt file:
file <- "mh.all.v2024.1.Mm.symbols.gmt"   #This is the master mouse hallmark gene set. 

#Read in and parse information from the gmt file.
pwl2 <- read.gmt(file)

#Cut down the gmt file to only include genes included in the selected dataset:
bg_genes <- pwl2[pwl2$gene %in% bg_genesIDs,]

#Create a function to run clusterProfiler over the input:
pipeline <- function(df){
  #Create a duplicate column of the gene names called "ID"
  df$ID <- df$Gene.name
  
  #Add another column that has the absolute value of log2FC
  df$absoluteLog2FoldChange <- abs(df$log2FoldChange)
  #Input the log2FC cutoff: 
  logFoldChangeCutoff <- 0.585 
  
  #Replace df with df with another variable that differentiates what is up/down regulated and significant:s
  df <- df %>% mutate(diffexpressed = case_when(
    log2FoldChange > logFoldChangeCutoff & padj < 0.05 ~ 'UP', #Up-regulated
    log2FoldChange < logFoldChangeCutoff & absoluteLog2FoldChange > logFoldChangeCutoff & padj < 0.05 ~ 'DOWN', #Down-regulated
    TRUE ~ 'NO' #If the two conditions above are not true, the default value is "NO"
  ))
  
  #Remove non-significant genes from dataset:
  df <- df[df$diffexpressed != 'NO',]
  #Remove duplicate values
  unique(df$diffexpressed)
  
  #Split the data into the genes up-regulated or down-regulated. 
  deg_results_list <- split(df, df$diffexpressed)
  
  #Run ClusterProfiler:
  #Create a variable called 'res', which is the results of the enricher() function being applied to all elements of the "deg_results_list" list (using lapply function) 
  res <- lapply(names(deg_results_list), 
                function(x) enricher(gene = deg_results_list[[x]]$ID, minGSSize = 0, TERM2GENE = bg_genes))
  
  #Rename the names of the res lists to include whether the genes are up or down regulated using deg_results_list variable. 
  names(res) <- names(deg_results_list)
  #Create a new variable called 'res_df' that contains a list of names of pathways and the results (e.g. log2FC, etc.)
  res_df <- lapply(names(res), function(x) rbind(res[[x]]@result))
  
  #Rename the names of the res_df lists to be that of res (e.g. pathways) 
  names(res_df) <- names(res)
  res_df <- do.call(rbind, res_df)
  head(res_df)
  
  #Convert the enrichResults to a dataframe with the pathways
  res_df <- lapply(names(res), function(x) rbind(res[[x]]@result))
  names(res_df) <- names(res)
  res_df <- do.call(rbind, res_df)
  
  res_df <- res_df %>% mutate(minuslog10padj = -log10(p.adjust),
                              diffexpressed = gsub("\\..*","",rownames(res_df)))
  
  # Subset to those pathways that have p adj < cutoff and gene count > cutoff
  padj_cutoff <- 0.05
  genecount_cutoff <- 10
  target_pws <- unique(res_df$ID[res_df$p.adjust < padj_cutoff & res_df$Count > genecount_cutoff]) 
  res_df <- res_df[res_df$ID %in% target_pws, ]
  
  #Return the results
  return(res_df)
}


#Run the function over the datasets
MA <- read.csv("Brady_DEG_MAvsMock.csv")
MA_Pathways <- pipeline(MA)
WT <- read.csv("Brady_DEG_WTvsMock.csv")
WT_Pathways <- pipeline(WT)

#Pathway cut off:
padj_cutoff <- 0.05 
minGeneCount <- 10

MA_Pathways_MinGC <- MA_Pathways %>% filter(Count >= minGeneCount & p.adjust <= padj_cutoff)
WT_Pathways_MinGC <- WT_Pathways %>% filter(Count >= minGeneCount & p.adjust <= padj_cutoff)

#Filter the results into two variables - one that has the upregulated pathways, and one that has the downregulated pathways that are significant based on cutoff. 
MA_Pathways_Up <- MA_Pathways_MinGC %>% filter(diffexpressed == 'UP' & p.adjust <= padj_cutoff)
WT_Pathways_Up <- WT_Pathways_MinGC %>% filter(diffexpressed == 'UP' & p.adjust <= padj_cutoff)
MA_Pathways_DOWN <- MA_Pathways_MinGC %>% filter(diffexpressed == 'DOWN' & p.adjust <= padj_cutoff)
WT_Pathways_DOWN <- WT_Pathways_MinGC %>% filter(diffexpressed == 'DOWN' & p.adjust <= padj_cutoff)

#Add a column that defines what condition the data is from
MA_Pathways_Up$Condition <- "MA"
WT_Pathways_Up$Condition <- "WT"
MA_Pathways_DOWN$Condition <- "MA"
WT_Pathways_DOWN$Condition <- "WT"

#Merge the files together vertically using the pathway names:
trVLP_UP_DF <- rbind(MA_Pathways_Up, WT_Pathways_Up)
trVLP_DOWN_DF <- rbind(MA_Pathways_DOWN, WT_Pathways_DOWN)

#Save the pathway results to a .csv file, if desired: 
#write.csv(trVLP_UP_DF, "Name1.csv", row.names=FALSE)
#write.csv(trVLP_DOWN_DF, "Name2.csv", row.names=FALSE)

#########################################################################################################################
#Making a dotplot from Enricher() Dotplot():
enrichres <- new("enrichResult",
                 readable = FALSE,
                 result = trVLP_UP_DF,  ###Input here
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.1,
                 organism = "mouse",
                 ontology = "UNKNOWN",
                 gene = bg_genes$gene,
                 keytype = "UNKNOWN",
                 universe = unique(bg_genes$gene),
                 geneSets = bg_genes)
class(enrichres)
dotplot(enrichres, 
        x = "Condition",
        color = "p.adjust",
        showCategory=15, 
        split = "Condition",
        label_format = 50,
        title='trVLP Upregulated Pathways')
