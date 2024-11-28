##############################################################################################
### Total DEG Bar Graph ###
##############################################################################################
#Begin by cleaning up the RStudio environment. 
rm(list = ls(all.names = TRUE))

#Load the required libraries: 
library(ggplot2)    #Version 3.5.1
library(dplyr)      #Version 1.1.4

#Read in the DEG files (from DESeq2):
MA_df <- read.csv("Brady_DEG_MAvsMock.csv")
WT_df <- read.csv("Brady_DEG_WTvsMock.csv")

#Set two variables to represent the log2 fold change minimum values. 
l2FC_Up <- 0.585
l2FC_Down <- -0.585

#Create an additional column on each dataframe that will denote whether a gene meets significance criteria:
#The default is 'NO' - eg the gene does not meet significance criteria
MA_df$diffexpressed <- 'NO' #Set the default value to 'NO'
#Change the value of the diffexpressed column of a given gene if the gene meets significance criteria, and if so, note if it is upregulated or downregulated:
MA_df <- MA_df %>% mutate(diffexpressed = case_when(
  log2FoldChange > l2FC_Up & padj < 0.05 ~ 'UP', #Upregulated
  log2FoldChange < l2FC_Down & padj < 0.05 ~ 'DOWN', #Downregulated
))
#Repeat for the wildtype data:
WT_df$diffexpressed <- 'NO' 
WT_df <- WT_df %>% mutate(diffexpressed = case_when(
  log2FoldChange > l2FC_Up & padj < 0.05 ~ 'UP', #Upregulated
  log2FoldChange < l2FC_Down & padj < 0.05 ~ 'DOWN', #Downregulated
))

#Remove genes that don't meet significance criteria:
MA_df <- MA_df[MA_df$diffexpressed != 'NO',]
WT_df <- WT_df[WT_df$diffexpressed != 'NO',]

#Remove NA values
MA_df <- na.omit(MA_df)
WT_df <- na.omit(WT_df)

#Count the number of upregulated and downregulated genes: 
MA_UP <- as.integer(nrow(MA_df[MA_df$diffexpressed == "UP",]))
MA_DOWN <- as.integer(nrow(MA_df[MA_df$diffexpressed == "DOWN",]))
WT_UP <- as.integer(nrow(WT_df[WT_df$diffexpressed == "UP",]))
WT_DOWN <- as.integer(nrow(WT_df[WT_df$diffexpressed == "DOWN",]))

#Create DEG Table
DEGTable = setNames(data.frame(matrix(ncol=3,nrow=2)),c("Sample","Regulation","DEG"))
DEGTable[1,1] <- "WT"
DEGTable[2,1] <- "WT"
DEGTable[3,1] <- "MA"
DEGTable[4,1] <- "MA"

DEGTable[1,2] <- "Upregulated"
DEGTable[2,2] <- "Downregulated"
DEGTable[3,2] <- "Upregulated"
DEGTable[4,2] <- "Downregulated"

DEGTable[1,3] <- WT_UP
DEGTable[2,3] <- WT_DOWN
DEGTable[3,3] <- MA_UP
DEGTable[4,3] <- MA_DOWN

#Make the downregulated genes negative so that the bars go down on the bar graph!
DEGTable <- DEGTable %>%
  mutate(DEG = ifelse(Regulation == "Upregulated",DEG,-1*DEG))

#Create the break values for the graphs:
breaks_values <- pretty(DEGTable$DEG)

#Create the plot:
p <- DEGTable %>%
  ggplot(aes(x=factor(Sample,level=c('WT', 'MA')),y=DEG, fill=Regulation))+
  geom_bar(stat = "identity", width=0.5)+
  geom_text(aes(label=abs(DEG),
                vjust = ifelse(Regulation == "Upregulated",-0.5,1.5)))

p + 
  scale_y_continuous(
    breaks = breaks_values,
    limits = c(-1500,1500),
    labels=abs(breaks_values))+
  theme(
    panel.background = element_rect(fill='transparent'),
    axis.line = element_line(color = "black"))+
  geom_hline(yintercept=0)+
  scale_fill_manual(values=c('#D73027','#4575B4'), breaks=c("Upregulated", "Downregulated"))+
  labs(x = "Condition", y = "Number of DEGs")
