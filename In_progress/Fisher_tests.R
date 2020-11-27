###~~~
#FISHER tests
###~~~

###~~~
#Treatment effect
###~~~

#based on `survheigh`

#Clean dataset
survheig <- survheig[-which(is.na(survheig$Cluster) == T),]

#Gather data
Treat <- unique(survheig$TreatConc)[c(4,5,3,2,1)]

###~~~
#Build matrix for Fisher test: Dead, live and # shoot tips per treatment
###~~~
Fisher_treat <- matrix(ncol=3, nrow=length(Treat))
rownames(Fisher_treat) <- Treat
colnames(Fisher_treat) <- c("Dead", "Live", "N. shoot tips")

for(i in 1:nrow(Fisher_treat)){
  #Subset matrix to include target treatment
  tmp <- subset(survheig, survheig$TreatConc == rownames(Fisher_treat)[i])
  
  #Pivot table
  FiveWsurvAll <- table(tmp$X5_w_survival)
  
  #Number of shoot tips
  Fisher_treat[i,3] <- sum(FiveWsurvAll)

  #Dead
  deadCol <- which(names(FiveWsurvAll) == 0)
  if(length(deadCol) > 0){
    Fisher_treat[i,1] <- as.vector(FiveWsurvAll[deadCol])
  }else{
    Fisher_treat[i,1] <- 0
  }
  
  #Live
  if(length(deadCol) > 0){
    Fisher_treat[i,2] <- sum(FiveWsurvAll) - as.vector(FiveWsurvAll[deadCol])
  }else{
    Fisher_treat[i,2] <- sum(FiveWsurvAll)
  }
}

#Perform Fisher test
fisher.test(Fisher_treat[,c(1,2)])


###~~~
#Cluster effect
###~~~

#Gather data
Cluster <- unique(survheigTreat$Cluster)

###~~~
#Build matrix for Fisher test: Dead, live and # shoot tips per treatment
###~~~
Fisher_cluster <- matrix(ncol=3, nrow=length(Cluster))
rownames(Fisher_cluster) <- Cluster
colnames(Fisher_cluster) <- c("Dead", "Live", "N. shoot tips")

for(i in 1:nrow(Fisher_cluster)){
  #Subset matrix to include target treatment
  tmp <- subset(survheig, survheig$Cluster == rownames(Fisher_cluster)[i])
  
  #Pivot table
  FiveWsurvAll <- table(tmp$X5_w_survival)
  
  #Number of shoot tips
  Fisher_cluster[i,3] <- sum(FiveWsurvAll)
  
  #Dead
  deadCol <- which(names(FiveWsurvAll) == 0)
  if(length(deadCol) > 0){
    Fisher_cluster[i,1] <- as.vector(FiveWsurvAll[deadCol])
  }else{
    Fisher_cluster[i,1] <- 0
  }
  
  #Live
  if(length(deadCol) > 0){
    Fisher_cluster[i,2] <- sum(FiveWsurvAll) - as.vector(FiveWsurvAll[deadCol])
  }else{
    Fisher_cluster[i,2] <- sum(FiveWsurvAll)
  }
}

#Perform Fisher test
fisher.test(Fisher_cluster[,c(1,2)])
