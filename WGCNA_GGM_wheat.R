#### Set parameters ####

normalizedData<- read.table(file="Development_tpm_HC.tsv", sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)
subData <- normalizedData[,c("Sample_19A","Sample_20A","Sample_21A","Sample_52A","Sample_53A","Sample_54A")]
subData_1 <-unique(subData[subData>0.5,])
normalizedData<-subData_1
nThreads=4
percThreshold=0.15
fileName="Development_tpm_HC.tsv"
path="./"

### Integrating networks ###

S1N=0.15
S2N=1
dataPath="."
geneNames = NULL

#### WGCNA ####

library(WGCNA)
library(corpcor)

options(stringsAsFactors = FALSE,future.globals.maxSize = 4000 * 1024^7)
normalizedData[1,1] <- as.numeric(normalizedData[1,1])
normalizedDataMatrix <- data.matrix(t(normalizedData))
randomMatrix <- matrix(runif((ncol(normalizedDataMatrix)*nrow(normalizedDataMatrix)),min=1e-20,max=2e-20),nrow=nrow(normalizedDataMatrix))
normalizedDataMatrix <- normalizedDataMatrix + randomMatrix
geneNames <- colnames(normalizedDataMatrix)

allowWGCNAThreads(nThreads = as.integer(nThreads))
powers = seq(from =2,to=20,by=2)
sft <- pickSoftThreshold(normalizedDataMatrix, powerVector = powers, verbose = 5)
R2value <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
maxIndex <- which.max(R2value)
optimalPower <- sft$fitIndices[maxIndex,1]
adjacencyMatrix <- adjacency(normalizedDataMatrix,power=optimalPower)
TOMMatrix <- TOMsimilarity(adjacencyMatrix)
TOMVector <- sm2vec(TOMMatrix)
indexes = sm.index(TOMMatrix)
edgeList = cbind(indexes,TOMVector)
sortIndex <- order(edgeList[,3],decreasing=TRUE)
edgeList <- edgeList[sortIndex,]

# permutation
randomRowIndexes <- sample(nrow(normalizedData),nrow(normalizedData),replace=FALSE)
randomMatrix <- normalizedData[randomRowIndexes,]
for(i in seq(1,length(randomRowIndexes))){
  randomColIndexes <- sample(ncol(normalizedData),ncol(normalizedData),replace=FALSE)
  randomMatrix[i,] <- randomMatrix[i,randomColIndexes] 
}
randomMatrix[1,1] <- as.numeric(randomMatrix[1,1])
randomMatrix <- data.matrix(t(randomMatrix))
randomsft <- pickSoftThreshold(randomMatrix, powerVector = powers, verbose = 5)
R2value <- -sign(randomsft$fitIndices[,3])*randomsft$fitIndices[,2]
maxIndex <- which.max(R2value)
optimalPower <- randomsft$fitIndices[maxIndex,1]
randomAdjacencyMatrix <- adjacency(randomMatrix,power=optimalPower)
randomTOMMatrix <- TOMsimilarity(randomAdjacencyMatrix)
randomTOMVector <- sm2vec(randomTOMMatrix)
sortIndex <- order(randomTOMVector,decreasing=TRUE)
randomTOMVector <- randomTOMVector[sortIndex]
threValue <- randomTOMVector[round(as.numeric(percThreshold)*length(randomTOMVector))]
aboveIndexes <- which(edgeList[,3] > threValue)
edgeList <- edgeList[aboveIndexes,]
fileName <- gsub("\\.txt$","",fileName)
fileName <- gsub(".*\\/", "", fileName)
save(edgeList,file=paste(path,paste(fileName,"_wgcna_network_final_edgeIndex.RData", sep=""),sep="/"))

library(corpcor)
wgcnaAdjacencyMatrix <- matrix(0,nrow=length(geneNames),ncol=length(geneNames))
wgcnaBoolMatrix <- matrix(1,nrow=length(geneNames),ncol=length(geneNames))
wgcnaNumMatrix <- matrix(0,nrow=length(geneNames),ncol=length(geneNames))
wgcnaNets <- list.files(dataPath,pattern="wgcna_network_final_edgeIndex\\.RData$",full.name=TRUE)

for(wgcnaNetIndex in seq(1,length(wgcnaNets))){
  tempMatrix <- matrix(0,nrow=length(geneNames),ncol=length(geneNames))
  edgeList[,3] <- (edgeList[,3]-min(edgeList[,3]))/(max(edgeList[,3])-min(edgeList[,3]))
  tempMatrix[as.matrix(edgeList[,1:2])] <- edgeList[,3]
  wgcnaAdjacencyMatrix <- wgcnaAdjacencyMatrix + tempMatrix
  wgcnaNumMatrix[as.matrix(edgeList[,1:2])] <- wgcnaNumMatrix[as.matrix(edgeList[,1:2])] + 1
}
wgcnaAdjacencyMatrix <- wgcnaAdjacencyMatrix/(length(wgcnaNets)-1)

#### ggm ####

library(GeneNet)
library(corpcor)

normalizedData<- read.table(file="Development_tpm_HC.tsv", sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)
subData <- normalizedData[,c("Sample_19A","Sample_20A","Sample_21A","Sample_52A","Sample_53A","Sample_54A")]
subData_1 <-unique(subData[subData>0.5,])
normalizedData<-subData_1


### Estimate partial correlation matrix
normalizedDataMatrix <- data.matrix(normalizedData)
normalizedDataMatrix <- t(normalizedDataMatrix)
randomMatrix <- matrix(runif((ncol(normalizedDataMatrix)*nrow(normalizedDataMatrix)),min=1e-20,max=2e-20),nrow=nrow(normalizedDataMatrix))
normalizedDataMatrix <- normalizedDataMatrix + randomMatrix
normalizedDataMatrix <- normalizedDataMatrix[,colSums(is.na(normalizedDataMatrix))<nrow(normalizedDataMatrix)]
partialCorMatrix <- ggm.estimate.pcor(normalizedDataMatrix)
pcor = sm2vec(partialCorMatrix)
indexes = sm.index(partialCorMatrix)
edgeList_ggm = cbind(indexes,pcor)
sortIndex <- order(edgeList_ggm[,3],decreasing=TRUE)
edgeList_ggm <- edgeList_ggm[sortIndex,]
fdr.out = fdrtool(edgeList_ggm[,3],statistic="correlation",plot=FALSE)
pval = fdr.out$pval
qval = fdr.out$qval
prob = 1 - fdr.out$lfdr
pcor <- edgeList_ggm[,3]

### permutation
randomRowIndexes <- sample(nrow(normalizedData),nrow(normalizedData),replace=FALSE)
randomMatrix <- normalizedData[randomRowIndexes,]
for(i in seq(1,length(randomRowIndexes))){
  randomColIndexes <- sample(ncol(normalizedData),ncol(normalizedData),replace=FALSE)
  randomMatrix[i,] <- randomMatrix[i,randomColIndexes] 
}
randomMatrix[1,1] <- as.numeric(randomMatrix[1,1])
randomMatrix <- data.matrix(t(randomMatrix))
randomMatrix <- randomMatrix[,colSums(is.na(randomMatrix))<nrow(randomMatrix)]
randomCorMatrix <- ggm.estimate.pcor(randomMatrix)

randomcor = sm2vec(randomCorMatrix)
posrandomcor <- randomcor[which(randomcor >= 0)]
negrandomcor <- abs(randomcor[which(randomcor < 0)])
sortIndex1 <- order(posrandomcor,decreasing=TRUE)
posrandomcor <- posrandomcor[sortIndex1]
sortIndex2 <- order(negrandomcor,decreasing=TRUE)
negrandomcor <- negrandomcor[sortIndex2]

threValue1 <- posrandomcor[round(as.numeric(percThreshold)*length(posrandomcor))]
threValue2 <- negrandomcor[round(as.numeric(percThreshold)*length(negrandomcor))]
cat(fileName,threValue1,threValue2,"\n")

aboveIndexes <- which(edgeList_ggm[,3] > threValue1 | edgeList_ggm[,3] < -threValue2)
edgeList_ggm <- edgeList_ggm[aboveIndexes,]
fileName <- gsub("\\.txt$","",fileName)
fileName <- gsub(".*\\/", "", fileName)
save(edgeList_ggm,file=paste(path,paste(fileName,"_ggm_network_final_edgeIndex.RData", sep=""),sep="/"))

## ggm ##
ggmNets <- list.files(dataPath,pattern="Count_ggm_network_final_edgeIndex\\.RData$",full.name=TRUE)
ggmAdjacencyMatrix <- matrix(0,nrow=length(geneNames),ncol=length(geneNames))
ggmBoolMatrix <- matrix(1,nrow=length(geneNames),ncol=length(geneNames))
ggmNumMatrix <- matrix(0,nrow=length(geneNames),ncol=length(geneNames))
for(ggmNetIndex in seq(1,length(ggmNets))){
  tempMatrix <- matrix(0,nrow=length(geneNames),ncol=length(geneNames))
  edgeList_ggm[,3] <- (edgeList_ggm[,3]-min(edgeList_ggm[,3]))/(max(edgeList_ggm[,3])-min(edgeList_ggm[,3]))
  tempMatrix[as.matrix(edgeList_ggm[,1:2])] <- edgeList_ggm[,3]
  ggmAdjacencyMatrix <- ggmAdjacencyMatrix + tempMatrix
  ggmNumMatrix[as.matrix(edgeList_ggm[,1:2])] <- ggmNumMatrix[as.matrix(edgeList_ggm[,1:2])] + 1
}

ggmAdjacencyMatrix <- ggmAdjacencyMatrix/(length(ggmNets)-1)

## Combine and write ##

adjacencyMatrix <- (wgcnaAdjacencyMatrix + ggmAdjacencyMatrix)/2
boolMatrix <- wgcnaBoolMatrix  + ggmBoolMatrix
candiIndexes <- which(boolMatrix >= S2N,arr.ind=TRUE)
cat("candiIndexes:",candiIndexes[1,],"\n")
edgeList <- cbind(candiIndexes,adjacencyMatrix[candiIndexes])
cat("edgeList:",edgeList[1:5,3],"\n")
sortIndex <- order(edgeList[,3],decreasing=TRUE)
edgeList <- edgeList[sortIndex,]
edgeList[,3] <- (edgeList[,3]-min(edgeList[,3]))/(max(edgeList[,3])-min(edgeList[,3]))
# save(edgeList,file=paste(dataPath,"final_geneNet_edgeIndex.RData",sep="/"))
netMatrix <- cbind(geneNames[edgeList[,1]],geneNames[edgeList[,2]],edgeList[,3])
write.table(netMatrix,file=paste(dataPath,"final_geneNet.tab",sep="/"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")