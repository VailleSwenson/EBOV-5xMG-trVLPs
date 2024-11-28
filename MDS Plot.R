######################################################################################################
### Multidimensional Scaling (MDS) Plot ###
######################################################################################################
#Begin by cleaning up the RStudio environment. 
rm(list = ls(all.names = TRUE))

#Load the required libraries.
library(cluster)  #Used version 2.1.6
library(limma)    #Used version 3.58.1

#Read in a .csv file containing the raw count values for all genes and all samples. 
counts <- read.csv("raw_counts.csv")

#Remove the columns that do not contain count data.
counts <- subset(counts, select = -c(1,11))

#Create a vector containing the desired shape for each replicate.
Shape <- c(15, 15, 15, 16, 16, 16, 17, 17, 17)

#Create a vector containing the desired color for each replicate. 
Color <- c("#44AA99", "#44AA99", "#44AA99", "#DDCC77", "#DDCC77", "#DDCC77", "#882255", "#882255", "#882255")

#Call the plotMDS() function from limma
plotMDS(counts, col=Color, pch=Shape, gene.selection = "common", xlab="Dimension 1", ylab="Dimension 2")

#Add a legend to the MDS plot.
legend("topright", col=c("#44AA99", "#DDCC77", "#882255"), pch=c(15,16,17), legend=c("MA", "Mock", "WT"))
