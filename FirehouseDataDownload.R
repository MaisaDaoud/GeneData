#  For more about the data and survival analysis https://www.biostars.org/p/153013/
setwd("/User/myUser/tutorial")
library(survival)

# read RNA file 
rna <- read.table('RNA/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt',nrows=20533, header=T,row.names=1,sep='\t')
# and take off first row cause we don't need it
rna <- rna[-1,]

# first remove genes whose expression is == 0 in more than 50% of the samples:
rem <- function(x){
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  r <- as.numeric(apply(x,1,function(i) sum(i == 0)))
  remove <- which(r > dim(x)[2]*0.5)
  return(remove)
}
remove <- rem(rna)
rna <- rna[-remove,]
# see the values
table(substr(colnames(rna),14,14))
# 0      1 
# 534   72
# get the index of the normal/control samples
n_index <- which(substr(colnames(rna),14,14) == '1')
t_index <- which(substr(colnames(rna),14,14) == '0')


# apply voom function from limma package to normalize the data
vm <- function(x){
  cond <- factor(ifelse(seq(1,dim(x)[2],1) %in% t_index, 1,  0))
  d <- model.matrix(~1+cond)
  x <- t(apply(x,1,as.numeric))
  ex <- voom(x,d,plot=F)
  return(ex$E)
}

rna_vm  <- vm(rna)
colnames(rna_vm) <- gsub('\\.','-',substr(colnames(rna),1,12))


# and check how data look, they should look normally-ish distributed
hist(rna_vm)
# we can remove the old "rna" cause we don't need it anymor
rm(rna)
#z = [(value gene X in tumor Y)-(mean gene X in normal)]/(standard deviation X in normal)


# calculate z-scores
scal <- function(x,y){
  mean_n <- rowMeans(y)  # mean of normal
  sd_n <- apply(y,1,sd)  # SD of normal
  # z score as (value - mean normal)/SD normal
  res <- matrix(nrow=nrow(x), ncol=ncol(x))
  colnames(res) <- colnames(x)
  rownames(res) <- rownames(x)
  for(i in 1:dim(x)[1]){
    for(j in 1:dim(x)[2]){
      res[i,j] <- (x[i,j]-mean_n[i])/sd_n[i]
    }
  }
  return(res)
}
z_rna <- scal(rna_vm[,t_index],rna_vm[,n_index])
# set the rownames keeping only gene name
rownames(z_rna) <- sapply(rownames(z_rna), function(x) unlist(strsplit(x,'\\|'))[[1]])

rm(rna_vm) #we don't need it anymore
