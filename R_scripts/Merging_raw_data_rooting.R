####----
# Code to merge rooting raw data spreadsheets
# and preparing data for statistical analyses
####-----

#MERGED FINAL DATASET
# Following cols: 
# - Block
# - Treatment
# - Replicate
# - Genotype
# - Individual
# - Callus
# - Root

###~~
#List all csv files (raw data from blocks)
###~~~
csv <- list.files(pattern='.csv')

###~~~
#Execute loop to process all files and merge them (see below)
###~~~
# Empty object that will contain all processed data
MERGE <- NULL
for(i in 1:length(csv)){
  ###
  #Read in sv file
  mat <- read.csv(csv[i])
  
  ###~~~
  #Create final matrix
  ###~~~
  
  #List of individuals
  indcsv <- LETTERS[seq( from = 1, to = 9 )]
  
  #Empty matrix
  FINALmat <- data.frame(matrix(ncol=7, nrow = nrow(mat)*length(indcsv)))
  colnames(FINALmat) <- c("Block","Treatment", "Replicate", "Genotype", "Individual", "Callus", "Root")
  
  ###~~~
  #Start populating FINALmat 
  ###~~~
  # Add Block, Treatment and Replicates
  #Get data for block, treatment (and replicates)
  blockTreatRep <- rep(as.vector(mat$X), 9)
  FINALmat$Block <- sapply(strsplit(blockTreatRep, split='_'), "[", 1)
  FINALmat$Treatment <- paste(sapply(strsplit(blockTreatRep, split='_'), "[", 2), sapply(strsplit(blockTreatRep, split='_'), "[", 3), sep='_')
  FINALmat$Replicate <- sapply(strsplit(blockTreatRep, split='_'), "[", 4)
  
  #Deal with Controls
  contr <- grep("Control_", FINALmat$Treatment)
  FINALmat$Replicate[contr] <- sapply(strsplit(FINALmat$Treatment[contr], split="_"), "[",2)
  FINALmat$Treatment[contr] <- sapply(strsplit(FINALmat$Treatment[contr], split="_"), "[",1)
  
  ###~~~
  # Fetch individual data
  ###~~~
  
  #Where are ind in mat?
  indcol <- match(indcsv, colnames(mat))
  
  #Fetch raw data for each individual in block
  OUT <- NULL
  for(j in 1:length(indcol)){
    #Extract info for each individual
    tmp <- mat[,indcol[j]:(indcol[j]+2)]
    colnames(tmp) <- c("Ind", "Callus", "Root")
    OUT <- rbind(OUT, tmp)
  }
  
  ###~~~
  #Add OUT to FINALmat
  ###~~~
  
  #Genotypes
  FINALmat$Genotype <- sapply(strsplit(as.vector(OUT$Ind), split='_'), "[", 1)
  #Ind
  FINALmat$Individual <- paste(sapply(strsplit(as.vector(OUT$Ind), split='_'), "[", 2), sapply(strsplit(as.vector(OUT$Ind), split='_'), "[", 3), sep='_')
  #Callus
  FINALmat$Callus <- OUT$Callus
  #Root
  FINALmat$Root <- OUT$Root
  
  ###~~~
  #Save FINALmat
  ###~~~
  write.csv(FINALmat, file=paste0("Processed_", csv[i]), row.names = F, quote = F)
  
  ###~~~
  #MERGE all csv
  ###~~~
  MERGE <- rbind(MERGE, FINALmat)
  
}

pdf("Rooting_results.pdf")
boxplot(MERGE$Root ~ MERGE$Treatment, xlab = "Treatment", ylab = "Number of roots")
dev.off()
###~~~
#Save MERGED files with all blocks
###~~~
write.csv(MERGE, file=paste0("Processed_", "blocks_rooting.csv"), row.names = F, quote = F)

###~~~
#ANOVAs
###~~~
#MERGE <- read.csv("Processed_blocks_rooting.csv")

MERGE <- data.frame(MERGE)

MERGE$IndGeno <- paste(MERGE$Genotype, MERGE$Individual, sep="_")

anovaOUT <- aov(Root ~ Treatment * IndGeno, data = MERGE)

summary(anovaOUT)

significant_vars <- c("Treatment", "IndGeno", "Treatment:IndGeno")
anovaOUT_tukey <- TukeyHSD(anovaOUT, significant_vars)


saveRDS(t2x_FvFm_aov1_tukey, "Final_Data/rds/t2x_FvFm_aov1_tukey.rds")

###---
# Are some individuals better performing?
###---

#Heatmap
require(RColorBrewer)
require(gplots)

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("green", "yellow", "red"))(n = 299)

mat_data <- table(MERGE$IndGeno, MERGE$Root)
colInd <- read.csv("Individuals.csv")
pdf("Heatmap_compare_ind_rooting.pdf")
out <- heatmap.2(mat_data,
          cellnote = mat_data,  # same data set for cell labels
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          notecex = .8, # cex of numbers in cells
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          xlab = 'N. roots per tip', #Label of x axis
          cexRow = .6, # cex of rows (= samples ID)
          margins =c(4,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row", # only draw a row dendrogram
          key=F,
          RowSideColors = as.vector(colInd$Group),
          Colv="NA") # turn off column clustering
dev.off()

Best <- subset(MERGE, MERGE$IndGeno %in% colInd$Individual[which(colInd$Group == "blue")])
anovaOUT2 <- aov(Root ~ Treatment * IndGeno, data = Best)

summary(anovaOUT2)

significant_vars <- c("Treatment", "IndGeno", "Treatment:IndGeno")

anovaOUT_tukey2 <- TukeyHSD(anovaOUT2, significant_vars)
