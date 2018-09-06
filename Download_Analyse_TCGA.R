
#devtools::install_github('BioinformaticsFMRP/TCGAbiolinks') install TCGAbiolinks
library(TCGAbiolinks)
library(dplyr)
#https://github.com/BioinformaticsFMRP/TCGAbiolinks/issues/156
#create a folder for the project and name it "GeneExpressionProject"
setwd('path /to/GeneExpressionProject')

# Query platform Illumina HiSeq with a list of barcode 
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  experimental.strategy = "RNA-Seq",
                  platform = "Illumina HiSeq",
                  file.type = "results",
                  legacy = TRUE)
# Download IF NOT ALREADY EXISTS, this downloads 1215 samples (16/10/2017)
GDCdownload(query,method='client')

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
dataBRCA <- GDCprepare(query)
# normalization of genes
dataNorm <- TCGAanalyze_Normalization(tabDF = dataBRCA, geneInfo =  geneInfo)
# ncol(dataNorm) # 1215 samples
# nroq(dataNorm) # 19858 Genes
# quantile filter of genes
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)
#get the tumor subtypes 
dataSubt <- TCGAquery_subtype(tumor = "BRCA")
#lumA <- dataSubt[which(dataSubt$PAM50.mRNA == "Luminal A"),1]

subtypes <- dataSubt$BRCA_Subtype_PAM50
head(subtypes)
count(dataSubt,BRCA_Subtype_PAM50)
######################
 #   x freq
#1  Basal  192
#2   Her2   82
#3   LumA  562
#4   LumB  209
#5     NA    2
#6 Normal   40
#####################

#####################GBM project ###########
''' 
  Original.Subtype     n
  <fct>            <int>
1 Classical          146
2 G-CIMP              39
3 Mesenchymal        157
4 Neural              83
5 Proneural           99
6 <NA>                82
'''
#############################################

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("NT"))

# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), 
                                   typesample = c("TP"))
#01A = Tumor sample,11A healthy sample some times 11A and 01A are from the same patient, 01B 
# means samples of the same aliquots from the same biopsy from the same patient,
samplesTP_with_known_subtype = dataSubt$patient[dataSubt$patient %in% substr(samplesTP,1,12)]
pat_id_label = dataSubt$patient[dataSubt$patient %in%samplesTP_with_known_subtype]
class_labels = dataSubt$BRCA_Subtype_PAM50[dataSubt$patient %in%samplesTP_with_known_subtype]
#dataFilt indecis
indexLst<-which(substr(colnames(dataFilt),1,15)%in%paste(samplesTP_with_known_subtype,"-01",sep=""))

dataFilt_update = dataFilt[,indexLst]
patient_ids = colnames(dataFilt)[indexLst]

                                
 nrow(dataFilt_update)
#[1] 14893
ncol(dataFilt_update)
#[1] 1083

match_class_data= class_labels[which(substr(patient_ids,1,12)%in%pat_id_label)]
# Diff.expr.analysis (DEA)  '#mat2 = dataFilt[,samplesTP],
dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                            mat2 = dataFilt_update,
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")

# DEGs table with expression values in normal and tumor samples TCGAanalyze_LevelTab(dataDEGs,"Tumor","Normal",
                                         # dataFilt[,samplesTP],dataFilt[,samplesNT])
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"Tumor","Normal",
                                          dataFilt_update,dataFilt[,samplesNT])
#nrow(dataDEGsFiltLevel) #3509 def. expressed genes

# take the indecies of the difretially expressed genes from (brca.exp) 
indicies<-which(row.names(dataNorm) %in% row.names(dataDEGsFiltLevel))
"
indecies<-matrix(ncol=1,nrow= nrow(dataDEGsFiltLevel))
i<-1

for(idx1 in 1:nrow(dataDEGsFiltLevel)){

 name1=row.names(dataDEGsFiltLevel)[idx1]
 Flag<-TRUE
 idx2<-1
 while(Flag==TRUE){
   
   name2<- row.names(dataNorm)[idx2]
   
   if(is.na(pmatch(name1,name2))==FALSE){
      indecies[i]<-idx2
      i<-i+1
      Flag<-FALSE
      }
   
   idx2<-idx2+1
   }}"
 #Z_score across NT genes 
 #NT_mat<-dataNorm[,samplesNT]
 #scaling <- function(gene_for_samples){
  #mean_n <- mean(gene_for_samples)  # mean of normal
  #sd_n <- sd(gene_for_samples)  # SD of normal
  # z score as (value - mean normal)/SD normal
 # for(i in 1:length(gene_for_samples)){
     # res[i] <- scale(gene_for_samples)#((gene_for_samples[i]-mean_n)/sd_n)}
  #return(res)}}
  samplesNT_mat<-dataNorm[indecies,samplesNT] #the normalized defrintially expressed normal samples

  samplesNT_mat2<-matrix(,nrow(samplesNT_mat),ncol(samplesNT_mat))
 for(idx in 1:nrow(samplesNT_mat)){
 samplesNT_mat2[idx,]<-scale(samplesNT_mat[idx,])} # zero mean, one unit variance
write.table(samplesNT_mat2, "samplesNT_scaled.csv",
            na = "",
            row.names = TRUE,
            col.names = TRUE,
            append = TRUE,
            sep = ",")
 #NT_mat2 : z_scores for normalized normal-samples
 #Z_score across TP samples genes
samplesTP_mat<-dataNorm[indicies,colnames(dataFilt_update)]
 #samplesTP_mat<-dataNorm[indecies,samplesTP]
 samplesTP_mat2<-matrix(,nrow(samplesTP_mat),ncol(samplesTP_mat))
 
  for(idx in 1:nrow(samplesTP_mat)){
 samplesTP_mat2[idx,]<-scale(samplesTP_mat[idx,])} # zero mean, one unit variance
write.table(samplesTP_mat2, "samplesTP_scaled.csv",
            na = "",
            row.names = FALSE,
            col.names = FALSE,
            append = TRUE,
            sep = ",")
 # if you want both the normal and tumor samples in the same file
 samples_mat<-dataNorm[indecies,colnames(dataFilt)]
 samples_mat2<-matrix(,nrow(samples_mat),ncol(samples_mat))
 for(idx in 1:nrow(samples_mat)){
 samples_mat2[idx,]<-scale(samples_mat[idx,])} 
 write.table(samples_mat2, "All_samples_scaled.csv",
            na = "",
            row.names = TRUE,
            col.names = TRUE,
            append = TRUE,
            sep = ",")
