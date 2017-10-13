library(TCGAbiolinks)
GDCdownload(query)

# Query platform Illumina HiSeq for all available barcodes (samples) 
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  experimental.strategy = "RNA-Seq",
                  platform = "Illumina HiSeq",
                  file.type = "results",
                  
                  legacy = TRUE)
# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
BRCARnaseqSE <- GDCprepare(query)
# normalization of genes
dataNorm <- TCGAanalyze_Normalization(tabDF = dataBRCA, geneInfo =  geneInfo)

# quantile filter of genes
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("NT"))

# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), 
                                   typesample = c("TP"))

# Diff.expr.analysis (DEA)
dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                            mat2 = dataFilt[,samplesTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")

# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"Tumor","Normal",
                                          dataFilt[,samplesTP],dataFilt[,samplesNT])

# take the indecies of the difretially expressed genes from (brca.exp) 
indecies<-matrix(ncol=1,nrow = 4)#nrow(dataDEGsFiltLevel))
i<-0
for(idx1 in 1:nrow(dataDEGsFiltLevel)){
 i<-i+1
 name1=row.names(dataDEGsFiltLevel)[i]
 
 for(idx2 in 1:nrow(brca.exp)){
   
   name2<- row.names(brca.exp)[idx2]
   
   if(is.na(pmatch(name1,name2))==FALSE){
   indecies[i]<-idx2}}}
 #Z_score across NT genes and Z_score across PT samples genes
 