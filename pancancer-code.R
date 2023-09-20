# CCCN-CFN-PCN
# September 19, 2023
# Yu Lu

library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
setwd("F:/luyu/PTM-metabolic/pancancer")
#raw data process
#download url :
pancancer <- read.table("acetylome_X.tsv",sep = "\t",header = T)
#ref_peptide <- read.table("LRG_RefSeqGene.txt", sep = "\t")
#ref_peptide <- ref_peptide[,c(3,8)]
l <- strsplit(pancancer[,1],"_")
t <- as.data.frame(l)
max(as.numeric(t[4,]))  #max = 2

#substr(l[[944]][3],1,nchar(l[[944]][3]) / 2 -1)
#substr(l[[944]][3],nchar(l[[944]][3]) / 2 + 1,nchar(l[[944]][3])-1)
#paste(l[[1]][1],l[[1]][2],sep="_") %in% ref_peptide[,2]


refseq_point <- data.frame()
refseq <- c()
for (i in 1:length(l)){
  tem <- paste(l[[i]][1],l[[i]][2],sep="_")
  tem_1 <- substr(tem,1,nchar(tem)-2)
  refseq <- c(refseq,tem_1)
  tem_frame <- data.frame(t(data.frame(c(tem,tem_1))))
  refseq_point <- rbind(refseq_point,tem_frame)
}

eg <- bitr(refseq,fromType = 'REFSEQ',
           toType = c('SYMBOL'),
           OrgDb='org.Hs.eg.db',
)

#substr(l[[944]][3],1,nchar(l[[944]][3]) / 2 -1)
#substr(l[[944]][3],nchar(l[[944]][3]) / 2 + 1,nchar(l[[944]][3])-1)
data_row_name <- data.frame()
rows_no_symbol <- c()  #record the rows that don't have symbol
for (i in 1:length(l)){
  tem <- paste(l[[i]][1],l[[i]][2],sep="_")
  tem_1 <- substr(tem,1,nchar(tem)-2)
  if (tem_1 %in% eg[,1]){
    sym <- eg[which(eg$REFSEQ==tem_1),2]
  }
  else {
    sym <- tem
    rows_no_symbol <- c(rows_no_symbol,i)
  }
  if (l[[i]][4] == 1){
    residue <- substr(l[[i]][3],1,nchar(l[[i]][3])-1)
    result_rowname <- paste(sym, "ack", residue, sep = " ")
  }
  if (l[[i]][4] == 2){#max = 2
    residue1 <- substr(l[[i]][3],1,nchar(l[[i]][3]) / 2 - 1)
    residue2 <- substr(l[[i]][3],nchar(l[[i]][3]) / 2 + 1,nchar(l[[i]][3])-1)
    result_rowname <- paste(sym, "ack", residue1, ";", sym, "ack", residue2, sep = " ")
  }
  data_row_name <- rbind(data_row_name,result_rowname)
}

pancancer1 <- cbind(data_row_name,pancancer)
pancancer1 <- pancancer1[,-2]
pancancer1 <- pancancer1[-rows_no_symbol,]
#duplicated
pancancer1 <- pancancer1[!duplicated(pancancer1[,1]),]
rownames(pancancer1) <- pancancer1[,1]
pancancer1 <- pancancer1[,-1]
write.csv(pancancer1,"acetylome_data.csv")

###################################################################
##Calculate Spearman dissimilarity, Euclidean distance, and SED
# Getting the raw data
result1 <- 2^pancancer1
gzdata.allt  <- result1

# correlation calculation
gzallt.cor <- cor(t(gzdata.allt), use = "pairwise.complete.obs", method = "spearman")
# Pearson Works, though there were many zero std. deviations
diag(gzallt.cor) <- NA
dissimilarity.gzallt <- 1 - abs(gzallt.cor)
diss.gzallt.noabs <- 1 - gzallt.cor
# set NA to two orders of magnitude higher than max distance
dissimilarity.gzallt[is.na(dissimilarity.gzallt)] <- 100*max(dissimilarity.gzallt, na.rm=T)
diss.gzallt.noabs[is.na(diss.gzallt.noabs)] <- 50*max(diss.gzallt.noabs, na.rm=T) 


# Euclid	
gzallt.dist = as.matrix (dist (gzdata.allt), method = "euclidean")  
# check 


### !!!!!!!! load Data_Input_Formatting.R   AND Dissimilarity_Calculations.R
max.na(gzallt.dist)	# Check for Inf values here
gzallt.dist[is.na(gzallt.dist)] <- 100*max(gzallt.dist, na.rm=T)
gzallt.dist.1 <- 100*gzallt.dist/max(gzallt.dist, na.rm=T)  # now max=100
# SED: combinde Euclid and Spearman w/o taking absolute value.	
gzallt.sed <- (gzallt.dist.1 + diss.gzallt.noabs)/2
#



## Calculate t-SNE embeddings from Spearman dissimilarity, Euclidean distance, and SED

##{r t-SNE embeddings}
# Rtsne: Barnes-Hut implementation of t-Distributed Stochastic Neighbor Embedding
library(Rtsne)
# Do with perplexity = 15 for slightly more resolution; theta = 0.25 for slightly more accuracy, max_iter=5000, for stabiliztion of groups

eu.gzallt.tsne.list <- Rtsne(as.matrix(gzallt.dist), dims = 3, perplexity = 15, theta = 0.25, max_iter=5000, check_duplicates = FALSE, pca=FALSE)
eu.gzallt.tsne <- eu.gzallt.tsne.list$Y
sp.gzallt.tsne.list <- Rtsne(as.matrix(dissimilarity.gzallt), dims = 3, perplexity = 15, theta = 0.25, max_iter=5000, check_duplicates = FALSE, pca=FALSE)
sp.gzallt.tsne <- sp.gzallt.tsne.list$Y
# SED without abs(cor)
sed.gzallt.tsne.list <- Rtsne(as.matrix(gzallt.sed), dims = 3, perplexity = 15, theta = 0.25, check_duplicates = FALSE, pca=FALSE)
sed.gzallt.tsne <- sed.gzallt.tsne.list$Y
#######################################################################################
# A view of the t-SNE embeddings can be obtained by plotting in three dimensions:
plot3d(eu.gzallt.tsne, type="s", radius=1.5, col="forestgreen")
plot3d(sp.gzallt.tsne, type="s", radius=1.5, col="green")
plot3d(sed.gzallt.tsne,  type="s", radius=1.5, col="red")

##########################################################################################
# Function make.clusterlist() is used to generate a list of PTM clusters from the t-SNE embedding
# Note that some tweaking of the toolong variable is necessary to obtain resonable cluster sizes.

eu.gzallt.list <- make.clusterlist(eu.gzallt.tsne, 3.8, gzdata.allt) 
esizes.gzallt <- sapply(eu.gzallt.list[[1]], function(x) dim(x)[1])

sp.gzallt.list <- make.clusterlist(sp.gzallt.tsne, 3.8, gzdata.allt) 
spsizes.gzallt <- sapply(sp.gzallt.list, function(x) dim(x)[1])

sed.gzallt.list <- make.clusterlist(sed.gzallt.tsne, 3.0, gzdata.allt) 
sedsizes.gzallt <- sapply(sed.gzallt.list, function(x) dim(x)[1])



# if there an error
##  esizes.gzallt <- sapply(eu.gzallt.list[[1]], function(x) dim(x)[1])
##  spsizes.gzallt <- sapply(sp.gzallt.list[[1]], function(x) dim(x)[1])
##  sedsizes.gzallt <- sapply(sed.gzallt.list[[1]], function(x) dim(x)[1])

hist(esizes.gzallt, breaks=100, col="red")
hist(spsizes.gzallt, breaks=100, col="blue")
hist(sedsizes.gzallt, breaks=100, col="green")

# if there an error




# Focus on intersect of all clusters from Euclid, Spearman, and SED using list.common()
#  Concatenate group list files
#		make group names unique
eu.gzallt.df <- ldply(eu.gzallt.list)[,2:3]  ##4449 2
sp.gzallt.df <- ldply(sp.gzallt.list)[,2:3] #4449 2
sed.gzallt.df <- ldply(sed.gzallt.list)[,2:3]	

### maybe there a bug(above),run this 
##  eu.gzallt.df <- ldply(eu.gzallt.list[[1]])[,2:3]  ##4449 2
##  sp.gzallt.df <- ldply(sp.gzallt.list[[1]])[,2:3] #4449 2
##  sed.gzallt.df <- ldply(sed.gzallt.list[[1]])[,2:3]	

# Further partition large groups
eu.gzallt.df $group <- paste(noquote(eu.gzallt.df $group), noquote("e"), sep="", collapse=NULL)
sp.gzallt.df $group <- paste(noquote(sp.gzallt.df $group), noquote("s"), sep="", collapse=NULL)
sed.gzallt.df $group <- paste(noquote(sed.gzallt.df $group), noquote("sed"), sep="", collapse=NULL)
gzalltgroups.df <- rbind(eu.gzallt.df, sed.gzallt.df, sp.gzallt.df)##13347 2
# For checking:
gzalltgroups.tab <- table(gzalltgroups.df)	
# Get the gene nemase from PTM clusters
eu.gzallt.genes <- lapply(eu.gzallt.list, extract.genes.from.clist)  

##if eu.gzallt.genes is null
##eu.gzallt.genes <- lapply(eu.gzallt.list[[1]], extract.genes.from.clist) 

sp.gzallt.genes <- lapply(sp.gzallt.list, extract.genes.from.clist)
## if sp.gzallt.genes is null ,run this
## sp.gzallt.genes <- lapply(sp.gzallt.list[[1]], extract.genes.from.clist)
sed.gzallt.genes <- lapply(sed.gzallt.list, extract.genes.from.clist)
##if sed.gzallt.genes is null ,run this
## sed.gzallt.genes <- lapply(sed.gzallt.list[[1]], extract.genes.from.clist)
# Extract the PTM names
eu.gzallt.peps <- lapply(eu.gzallt.list, extract.peps.from.clist)
sp.gzallt.peps <- lapply(sp.gzallt.list, extract.peps.from.clist)
sed.gzallt.peps <- lapply(sed.gzallt.list, extract.peps.from.clist)

#if null run this
eu.gzallt.peps <- lapply(eu.gzallt.list[[1]], extract.peps.from.clist)
sp.gzallt.peps <- lapply(sp.gzallt.list[[1]], extract.peps.from.clist)
sed.gzallt.peps <- lapply(sed.gzallt.list[[1]], extract.peps.from.clist)

eu.sp.gene <- list.common(eu.gzallt.genes, sp.gzallt.genes, keeplength=2)
eu.sp.gene.sizes <- sapply(eu.sp.gene, length)
eu.sp.sed.gene <- list.common(eu.sp.gene, sed.gzallt.genes, keeplength=2)
eu.sp.sed.gene.sizes <- sapply(eu.sp.sed.gene, length) 



eu.sp.gzallt <- list.common(eu.gzallt.peps, sp.gzallt.peps, keeplength=2)
eu.sp.gzallt.sizes <- sapply(eu.sp.gzallt, length)
eu.sp.sed.gzallt <- list.common(eu.sp.gzallt, sed.gzallt.peps, keeplength=2)
eu.sp.sed.gzallt.sizes <- sapply(eu.sp.sed.gzallt, length) 



