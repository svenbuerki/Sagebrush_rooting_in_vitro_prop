####----
# 1. PROCESS RAW DATA: Merge rooting raw data spreadsheets
# and prepare data for statistical analyses
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
csv <- list.files(pattern='block.csv')

#Rm processed files
toRm <- grep("Processed_", csv)
if(length(toRm) > 0){
  csv <- csv[-toRm]
}

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
  write.csv(FINALmat, file=paste0("Processed_data/Processed_", csv[i]), row.names = F, quote = F)
  
  ###~~~
  #MERGE all csv
  ###~~~
  MERGE <- rbind(MERGE, FINALmat)
  
}

#Finalize preping the data
MERGE <- data.frame(MERGE)

MERGE$IndGeno <- paste(MERGE$Genotype, MERGE$Individual, sep="_")

#Create a binary for Root
MERGE$RootBin <- MERGE$Root
MERGE$RootBin[MERGE$RootBin > 0] <- 1


###~~~
#Save MERGED files with all blocks
###~~~
write.csv(MERGE, file=paste0("Processed_data/Processed_", "blocks_rooting.csv"), row.names = F, quote = F)

####-----
# INFER SAMPLING STATS
####-----

#Create and populate table to summarize sampling (Table 1)
samp_mat <- matrix(ncol=9, nrow=1)
colnames(samp_mat) <- c("N. block", "N. treatments", "Treatments", "N. replicates", "N. genotypes", 
                        "N. individuals", "N. ind. per genotype", "Ind. G1", "Ind. G2")

#N block
samp_mat[,1] <- length(unique(MERGE$Block))
#N treat
samp_mat[,2] <- length(unique(MERGE$Treatment))
#Treat ID
samp_mat[,3] <- paste(sort(unique(MERGE$Treatment)), collapse = ", ")
#N replicates
samp_mat[,4] <- length(unique(MERGE$Replicate))
#N genotypes
samp_mat[,5] <- length(unique(MERGE$Genotype))
#N ind
samp_mat[,6] <- length(unique(MERGE$IndGeno))
#N ind per genotype
samp_mat[,7] <- paste(paste("G1:", list(table(MERGE$Genotype)/15)[[1]][1]), paste("G2:", list(table(MERGE$Genotype)/15)[[1]][2]), sep=', ')
#List ind G1
samp_mat[,8] <- paste(sort(unique(subset(MERGE$IndGeno, MERGE$Genotype == "G1"))), collapse=', ')
#List ind G2
samp_mat[,9] <- paste(sort(unique(subset(MERGE$IndGeno, MERGE$Genotype == "G2"))), collapse=', ')

#Write table
write.csv(samp_mat, file="Table_1_sampling_summary.csv", row.names = F)

###
#INFER Rooting stats
###

#What proportion of tips rooted?

100*(table(MERGE$RootBin)/length(MERGE$RootBin))
# Only 39% of tips rooted!
#0        1 
#60.59259 39.40741 

#Rooting by treatment
table(as.numeric(MERGE$RootBin), MERGE$Treatment)

#   Control IBA_05 IBA_1 NAA_05 NAA_1
#0     123     71    54     80    81
#1      12     64    81     55    54

#Table with freq, average, sd for each treatment
RootStats <- matrix(ncol=3, nrow=length(unique(MERGE$Treatment)))
colnames(RootStats) <- c("Freq. (0/1)","Avg.", "sd")

TreatTypes <- sort(unique(MERGE$Treatment))

rownames(RootStats) <- TreatTypes
for(i in 1:length(TreatTypes)){
  tmp <- subset(MERGE, MERGE$Treatment == TreatTypes[i])
  #Freq
  RootStats[i,1] <- paste(table(tmp$RootBin), collapse = "/")
  #Avg
  RootStats[i,2] <- round(mean(as.vector(tmp$RootBin)),2)
  #Sd
  RootStats[i,3] <- round(sd(as.vector(tmp$RootBin)),2)
  #quantile5%
  #quantile95%
  
}

x <- as.data.frame(RootStats)
library(ggplot2)

qplot(x    = Avg. , 
      y    = rownames(x),
      data = x) +
  
  geom_errorbar(aes( 
    ymin  = 0, 
    ymax  = 1, 
    width = 0.15))
#Plot
barplot(as.numeric(RootStats[,2]), names=rownames(RootStats), ylab="Average rooting", col='white', ylim = c(min(as.numeric(RootStats[,2])-as.numeric(RootStats[,3])), max(as.numeric(RootStats[,2])+as.numeric(RootStats[,3]))))

segments(y1=as.numeric(RootStats[,2])-as.numeric(RootStats[,3]), y0=as.numeric(RootStats[,2])+as.numeric(RootStats[,3]), x0=c(0.5,1.5,2.5,3.5,4.5), x1=c(0.5,1.5,2.5,3.5,4.5))

###~~~
#QUESTIONS
###~~~

###
#Q1: Is there a treatment that significantly increases rooting of sagebrush shoot tips?
# First, we have to test whether there is a significant treatment effect in explaining rooting data using ANOVAs.
# If yes, then we will classify treatments based on Tukey tests.
# Answer: Yes, IBA_1 is significantly outperforming the other treatments (but there is also high individual effects, see Q2)
#
#Q2: Are some individuals more efficent at rooting, which should be proritized to establish lines?
# First, we have to establish whether there is a significant individual effect in explaining the rooting data using ANOVAs. 
# If yes, then we will identify better performing individuals using:
# i) Tukey tests (pval: 0.01),
# ii) To qualify as top performer an ind. has to outperfom >=20% of individuals. 
# iii) Heatmap analysis on root data where individuals are sorted based on a clustering analysis.  
# Answer: Yes, as shown in Q1, there is a significant individual effect.
# Individuals are clustering into 3 groups (and there are 3 top performers: G2_b27_1, G2_b24_1, G2_b4_1): 
# - Black: No rooting independently of treatment
# - Pink: Some rooting, but majority of tips are not producing any roots.
# - Blue: Significant rooting with IBA_1 being the best treatment (where top performers are).
###

###
# Additional analyses:
# Observations: During the experiment, calli developed on shoot tips, which was unexpected. 
# --> This might be due to age of plants, but it requiers additional analyses. 
# More specifically, we are testing whether callus formation is treatment specific using ANOVAs (same model than with roots)
# See bottom of document for more details on analyses.

####----
# 2. STATISTICAL ANALYSES: ANOVAs and Tukey tests
####-----

###
#Q1: Is there a treatment that increases rooting of shoot tips?
###
#To answer that question, we are conducting ANOVA analyses followed by Tukey tests on significant variables.

#Our overall model is: root ~ treatment
# But we are also testing if there is a block, genotype and individual effect in rooting as well as interactions between these variables

###~~~
#2.1 ANOVAs with model
###~~~

###
#Root quantitative
###

#Defining model and runing analyses
anovaOUT <- aov(Root ~ Treatment * Block * Genotype * IndGeno, data = MERGE)
#Summary
summary(anovaOUT)
#Save analysis
saveRDS(anovaOUT,"rds/anova.rds")

#####
#RESULTS ANOVA: Print screen
#                     Df Sum Sq Mean Sq F value   Pr(>F)    
# Treatment            4  217.2   54.29  26.066  < 2e-16 ***
# Block                4   23.4    5.84   2.805  0.02538 *  
# Genotype             1    0.0    0.01   0.004  0.94814    
# IndGeno             39  642.4   16.47   7.908  < 2e-16 ***
# Treatment:Block     16  106.3    6.65   3.191 3.18e-05 ***
# Treatment:Genotype   4    6.3    1.57   0.753  0.55665    
# Treatment:IndGeno  156  454.1    2.91   1.397  0.00426 ** 
# Residuals          450  937.3    2.08                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#There are significant (at least **):
# - Treatment (***)
# - Individual (***)
# - Treatment:Block (***) and Treatment:Indvidual (**) --> These are due to the makeup of the blocks or individual variation within those (see Q2) 

###~~~
#2.2 Tukey tests on siginficant variables
###~~~

#Tukey tests
significant_vars <- c("Treatment", "IndGeno", "Treatment:IndGeno", "Treatment:Block")
anovaOUT_tukey <- TukeyHSD(anovaOUT, significant_vars)

#Save analysis
saveRDS(anovaOUT_tukey,"rds/anova_tukey.rds")

#What treatment(s) are significantly increasing rooting? (*: 0.01)
anovaOUT_tukey$Treatment[which(as.numeric(anovaOUT_tukey$Treatment[,4]) <= 0.01),]

####
# RESULTS Tukey tests: treatments
#                     diff        lwr         upr        p adj
# IBA_05-Control  1.1925926  0.7114648  1.67372039 3.441614e-10
# IBA_1-Control   1.6962963  1.2151685  2.17742409 0.000000e+00
# NAA_05-Control  0.6148148  0.1336870  1.09594261 4.635671e-03
# NAA_1-Control   0.8814815  0.4003537  1.36260928 7.452445e-06
# NAA_05-IBA_05  -0.5777778 -1.0589056 -0.09664998 9.526645e-03
# NAA_05-IBA_1   -1.0814815 -1.5626093 -0.60035369 1.644969e-08
# NAA_1-IBA_1    -0.8148148 -1.2959426 -0.33368702 4.518265e-05

#Ranking treatments (worst to best)
# - (c) Control (outperformed by all treatments) (***)
# - (b) NAA_1 and NAA_05
# - (a) IBA_1, IBA_05 (although IBA_1 bettern than IBA_05 with 0.05)

##
#Boxplot: Boxplot comparing rooting per treatment
pdf("Fig_1_Rooting_results.pdf")
boxplot(MERGE$Root ~ MERGE$Treatment, xlab = "Treatment", ylab = "Number of roots per tip", frame = "none", ylim=c(0,12.5))

#Add ranking
text(x=c(1,2,3,4,5), y=rep(12.5,5), c("c","a","a",'b','b'))

#Add ID of outliers and letters (see above)
dev.off()


###
#Root binary
###

#Defining model and runing analyses
anovaOUTbin <- aov(RootBin ~ Treatment * Block * Genotype * IndGeno, data = MERGE)
#Summary
summary(anovaOUTbin)
#Save analysis
saveRDS(anovaOUTbin,"rds/anovabin.rds")

#####
#RESULTS ANOVA: Print screen
#                     Df Sum Sq Mean Sq F value   Pr(>F)    
#Treatment            4  19.19   4.798  32.064  < 2e-16 ***
#Block                4   4.20   1.050   7.015 1.75e-05 ***
#Genotype             1   0.28   0.283   1.894  0.16948    
#IndGeno             39  34.56   0.886   5.923  < 2e-16 ***
#Treatment:Block     16   5.93   0.371   2.479  0.00124 ** 
#Treatment:Genotype   4   1.21   0.302   2.016  0.09127 .  
#Treatment:IndGeno  156  28.47   0.182   1.220  0.06007 .  
#Residuals          450  67.33   0.150                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#There are significant (at least **):
# - Treatment (***)
# - Block (***)
# - Individual (***)
# - Treatment:Block (***) and Treatment:Indvidual (**) --> These are due to the makeup of the blocks or individual variation within those (see Q2) 

###~~~
#2.2 Tukey tests on siginficant variables
###~~~

#Tukey tests
significant_vars <- c("Treatment", "IndGeno", "Treatment:IndGeno", "Treatment:Block")
anovaOUTbin_tukey <- TukeyHSD(anovaOUTbin, significant_vars)

#Save analysis
saveRDS(anovaOUT_tukey,"rds/anovabin_tukey.rds")

#What treatment(s) are significantly increasing rooting? (*: 0.01)
anovaOUTbin_tukey$Treatment[which(as.numeric(anovaOUTbin_tukey$Treatment[,4]) <= 0.01),]

####
# RESULTS Tukey tests: treatments
#                     diff        lwr         upr        p adj
# IBA_05-Control  1.1925926  0.7114648  1.67372039 3.441614e-10
# IBA_1-Control   1.6962963  1.2151685  2.17742409 0.000000e+00
# NAA_05-Control  0.6148148  0.1336870  1.09594261 4.635671e-03
# NAA_1-Control   0.8814815  0.4003537  1.36260928 7.452445e-06
# NAA_05-IBA_05  -0.5777778 -1.0589056 -0.09664998 9.526645e-03
# NAA_05-IBA_1   -1.0814815 -1.5626093 -0.60035369 1.644969e-08
# NAA_1-IBA_1    -0.8148148 -1.2959426 -0.33368702 4.518265e-05

#Ranking treatments (worst to best)
# - (c) Control (outperformed by all treatments) (***)
# - (b) NAA_1 and NAA_05
# - (a) IBA_1, IBA_05 (although IBA_1 bettern than IBA_05 with 0.05)

##
#Boxplot: Boxplot comparing rooting per treatment
pdf("Fig_1_Rooting_results.pdf")
boxplot(MERGE$Root ~ MERGE$Treatment, xlab = "Treatment", ylab = "Number of roots per tip", frame = "none", ylim=c(0,12.5))

#Add ranking
text(x=c(1,2,3,4,5), y=rep(12.5,5), c("c","a","a",'b','b'))

#Add ID of outliers and letters (see above)
dev.off()

###
#Q2: Is there individual variation in rooting rates?
# Yes, as shown in Q1, there is a significant individual effect.
# Individuals are clustering into 3 groups: 
# - Black: No rooting independently of treatment
# - Pink: Some rooting, but majority of tips are not producing any roots.
# - Blue: Significant rooting with IBA_1 being the best treatment.
###

# Ansewring this question is key since we would like to establish and maintain lines of each genotypes
# for genome assembly and GxE experiments.s

# To investigate this question a heatmap was inferred based on rooting data at individual level (from MERGE)
# and samples were sorted into groups based on a clustering analyses (using hclust) 

###~~~
#2.3 Tukey test: Individual effect
###~~~

#Individual effect
IndTukey <- anovaOUT_tukey$IndGeno[which(as.numeric(anovaOUT_tukey$IndGeno[,4]) <= 0.01),]

#What are the individuals that are the most prolific at producing roots?
#Those should be prioritized for establishing lines
IndTukeyOut <- IndTukey[which(IndTukey[,1] > 0),]
#Extract best performing ind (first in rownames of comparison)
P1 <- sapply(strsplit(rownames(IndTukey), split='-'), "[[",1)

#Percentage of individuals that are outperformed by P1
TopInd <- 100*(sort(table(P1), decreasing=T)/length(unique(MERGE$IndGeno)))

#Select only ind that are outperforming at least 20% of individuals 
TopInd <- TopInd[which(TopInd >= 20)]

###
#Results of top individuals
# Only 3 ind are fitting criterion and they are all G2!
# G2_b27_1 G2_b24_1  G2_b4_1 
# 75.55556 33.33333 20.00000 

###~~~
# 3. Heatmap of N of roots per tip (total of tips/ind = 15) 
# where individuals are sorted into groups using a clustering approach (by inferring a euclidian distance using rooting data)
###~~~
#Load requiered packages
require(RColorBrewer)
require(gplots)

# Creates our own color palette from green to red for heatmap
my_palette <- colorRampPalette(c("green", "yellow", "red"))(n = 299)

#Input data for heatmap
mat_data <- table(MERGE$IndGeno, MERGE$Root)

#Load file with color of individuals (which group based on clusters they belong to)
colInd <- read.csv("Individuals.csv")

#Annotate best perfoming individuals (based on Tukey test) and add P1, P2, P3
toAdd <- which(rownames(mat_data) %in% names(TopInd))
rownames(mat_data)[toAdd] <- paste0(rownames(mat_data)[toAdd], c(" (P1)", " (P2)", " (P3)"))

#Draw plot
pdf("Fig_2_Heatmap_compare_ind_rooting.pdf")
out <- heatmap.2(mat_data,
          cellnote = mat_data,  # same data set for cell labels
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          notecex = .8, # cex of numbers in cells
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          xlab = 'N. roots per tip (1 ind. has 15 tips)', #Label of x axis
          cexRow = .8, # cex of rows (= samples ID)
          margins =c(4,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row", # only draw a row dendrogram
          key=F,
          keysize = 0.4,
          RowSideColors = as.vector(colInd$Group),
          Colv="NA") # turn off column clustering
dev.off()

###~~~
#4. Are IBA_1 & IBA_05 still the best treatments for top performers?
###~~~~

###~~~
#4.1 ANOVAs with simplified model
###~~~

#Subset MERGE dataset to only include best ind
Best <- subset(MERGE, MERGE$IndGeno %in% names(TopInd))
#Perform ANOVA with simplified model
anovaOUT2 <- aov(Root ~ Treatment * IndGeno, data = Best)
#Summary
summary(anovaOUT2)
#Save analysis
saveRDS(anovaOUT2,"rds/anovaBest.rds")

####
# Results: There are only significant (*) treatment and Treatment:Ind effects
#                   Df Sum Sq Mean Sq F value Pr(>F)  
# Treatment          4 105.64  26.411   3.580 0.0168 *
# IndGeno            2  40.53  20.267   2.747 0.0803 .
# Treatment:IndGeno  8 137.69  17.211   2.333 0.0443 *
# Residuals         30 221.33   7.378                 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

###~~~
#4.2 Tukey tests on siginficant variables
###~~~

#Tukey tests
significant_vars <- c("Treatment", "Treatment:IndGeno")
anovaOUT_tukey2 <- TukeyHSD(anovaOUT2, significant_vars)
#Save analysis
saveRDS(anovaOUT_tukey2,"rds/TukeyRoot.rds")

#Investigate treatment effect (*)
TreatTukey2 <- anovaOUT_tukey2$Treatment[which(as.numeric(anovaOUT_tukey2$Treatment[,4]) <= 0.05),]
###
# Results
# IBA_1 outperforms control (0.01, **), but all others are not significant 
#                     diff        lwr       upr      p adj
# IBA_1-Control   4.5555556  0.8415243 8.2695868 0.01029258
##

####----
# 5. CALLUS DATA: STATISTICAL ANALYSES: ANOVAs and Tukey tests
####-----

#5.1. What proportion of tips generated calli?
#15*45 = Total number of shoot tips
100*(table(MERGE$Callus)/(15*45))
###
#Results
# 2/3 of tips generated calli!
#       0        1 
#33.62963 66.37037 
# 0: no callus, 1: callus present.

#5.2 Test whether there is a significant treatment or individual effect in explaining callus data
###
#ANOVA
###
#Perform ANOVA using same model as root
anovaCallus <- aov(Callus ~ Treatment * Block * Genotype * IndGeno, data = MERGE)
#Summary
summary(anovaCallus)
#Save analysis
saveRDS(anovaCallus,"rds/anovaCallus.rds")

####
# Results
# Treatment and Individual effects in explaining callus production
#                    Df Sum Sq Mean Sq F value   Pr(>F)    
# Treatment            4  68.88  17.221 200.414  < 2e-16 ***
# Block                4   6.66   1.665  19.379    1e-14 ***
# Genotype             1   0.14   0.141   1.643 0.200574    
# IndGeno             39  15.99   0.410   4.772  < 2e-16 ***
# Treatment:Block     16   3.49   0.218   2.537 0.000934 ***
# Treatment:Genotype   4   1.09   0.273   3.177 0.013641 *  
# Treatment:IndGeno  156  15.74   0.101   1.174 0.104393    
# Residuals          450  38.67   0.086                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

###
#Tukeys
###
significant_vars <- c("Treatment", "IndGeno")
anovaCallus_tukey2 <- TukeyHSD(anovaCallus, significant_vars)
TreatTukeyCallus <- anovaCallus_tukey2$Treatment[which(as.numeric(anovaCallus_tukey2$Treatment[,4]) <= 0.01),]
#Save analysis
saveRDS(TreatTukeyCallus,"rds/TukeyCallus.rds")

###
#Results
# Control tips didn't produce calli!
#                    diff       lwr       upr       p adj
# IBA_05-Control 0.7259259 0.6282063 0.8236455 0.000000000
# IBA_1-Control  0.8148148 0.7170952 0.9125344 0.000000000
# NAA_05-Control 0.7851852 0.6874656 0.8829048 0.000000000
# NAA_1-Control  0.8444444 0.7467248 0.9421641 0.000000000
# NAA_1-IBA_05   0.1185185 0.0207989 0.2162381 0.008543496

#Ranking treatments (worst to best)
# - (c) Control (outperformed by all treatments) (***) 
# - (b) Not significant differences between IBA_1, NAA_05
# - (a) NAA_1 outperforms IBA_05

###
#Heatmap
###

# Create our own color palette from green to red for heatmap
my_palette <- colorRampPalette(c("green", "yellow", "red"))(n = 299)

#Input data for heatmap
mat_data <- table(MERGE$IndGeno, MERGE$Callus)

#Load file with color of individuals (which group based on clusters they belong to)
colInd <- read.csv("Individuals.csv")

#Annotate best perfoming individuals (based on Tukey test) and add P1, P2, P3
toAdd <- which(rownames(mat_data) %in% names(TopInd))
rownames(mat_data)[toAdd] <- paste0(rownames(mat_data)[toAdd], c(" (P1)", " (P2)", " (P3)"))

#Draw plot
pdf("Fig_3_Heatmap_compare_ind_callus.pdf")
heatmap.2(mat_data,
                 cellnote = mat_callus,  # same data set for cell labels
                 main = "", # heat map title
                 notecol="black",      # change font color of cell labels to black
                 notecex = .8, # cex of numbers in cells
                 density.info="none",  # turns off density plot inside color legend
                 trace="none",         # turns off trace lines inside the heat map
                 xlab = 'Presence (1) / absence (0) of callus in tip', #Label of x axis
                 cexRow = .8,
          cexCol = .9,# cex of rows (= samples ID)
                 margins =c(4,9),     # widens margins around plot
                 col=my_palette,       # use on color palette defined earlier
                 #breaks=col_breaks,    # enable color transition at specified limits
                 dendrogram="row", # only draw a row dendrogram
                 key=F,
                 keysize = 0.4,
                 RowSideColors = as.vector(colInd$Group),
                 Colv="NA")
dev.off()


##Additional code

### Mean of measuremments by individuals

The objective is to average data at individual level by accounting for variations in the three replicates. We can do that because there is no replicate effect when ANOVA is run on raw data. The code below create a table displaying averaged rooting and callus data at indvidual level based on the three replicates (Table \@ref(tab:avgInd)):
  
  ```{r avgInd2, echo=T, eval=T, warning=F}
#List of Ind and treatments
IndList <- sort(unique(MERGE$IndGeno))
treat <- sort(unique(MERGE$Treatment))
#Create matrix with Ind and treatment (this will be used to subset MERGE and infer means)
MERGEind <- matrix(ncol=7, nrow=5*length(IndList))
colnames(MERGEind) <- c("Block", "Treatment", "Genotype", "IndGeno","Mean_Callus", "Mean_Root", "Mean_RootBin")
MERGEind <- as.data.frame(MERGEind)

#Populate Ind and treatments
MERGEind$IndGeno <- sort(rep(IndList, length(treat)))
MERGEind$Treatment <- rep(treat, length(IndList))

#Loop to infer means of root and callus
for(i in 1:nrow(MERGEind)){
  #Subset to individual
  tmp <- subset(MERGE, MERGE$IndGeno == MERGEind$IndGeno[i] & MERGE$Treatment == MERGEind$Treatment[i])
  #Block
  MERGEind$Block[i] <- unique(tmp$Block)
  
  #Genotype
  MERGEind$Genotype[i] <- unique(tmp$Genotype)
  #Mean Callus
  MERGEind$Mean_Callus[i] <- round(mean(tmp$Callus),2)
  #Mean root
  MERGEind$Mean_Root[i] <- round(mean(tmp$Root),2)
  #Mean root bin
  MERGEind$Mean_RootBin[i] <- round(mean(as.numeric(tmp$RootBin)),2)
  
}

#Plot table
knitr::kable(MERGEind, caption = "Mean rooting and callus data at individual level.")
```

Plot of mean rooting by treatment (based on binary data `Mean_RootBin`)

```{r plotmeanroot, echo=T, eval=T}
#Table with mean, sd for each treatment
#List treatment
TreatTypes <- sort(unique(MERGEind$Treatment))

RootStats <- matrix(ncol=2, nrow=length(TreatTypes))
colnames(RootStats) <- c("Mean_RootBin_treat", "sd")
rownames(RootStats) <- TreatTypes

for(i in 1:length(TreatTypes)){
  tmp2 <- subset(MERGEind, MERGEind$Treatment == TreatTypes[i])
  
  #Avg
  RootStats[i,1] <- round(mean(as.numeric(tmp2$Mean_RootBin)),2)
  #Sd
  RootStats[i,2] <- round(sd(as.numeric(tmp2$Mean_RootBin)),2)
}

#Plot
plot(as.numeric(RootStats[,1]), names=rownames(RootStats), ylab="Mean rooting at individual level", col='white', ylim = c(min(as.numeric(RootStats[,1])-as.numeric(RootStats[,2])), max(as.numeric(RootStats[,1])+as.numeric(RootStats[,2]))))

segments(y1=as.numeric(RootStats[,1])-as.numeric(RootStats[,2]), y0=as.numeric(RootStats[,1])+as.numeric(RootStats[,2]), x0=c(0.75,1.75,2.5,3.5,4.5), x1=c(0.75,1.75,2.5,3.5,4.5))

```


