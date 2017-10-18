# GeneData
This project introduces some basics about downloading gene expression datasets from GDC and firebrowse repositories. This project uses R for downloading and pre-processing the data.
 ###### To download Data from FireBrowse
-  http://gdac.broadinstitute.org/ 
-  Select "Breast invasive carcinoma" -> "data" -> "Browse".
-  From "mRNASeq" select "illuminahiseq_rnaseqv2-RSEM_genes_normalized" and save it on your computer.
-  Unzip the folder and rename it to RNA.
###### To Download and process data from GDC/TCGA
- TCGA repository  https://portal.gdc.cancer.gov/ 
- use  Download_Analyze_TCGA.R to download and analyse the data
-The code normalize the data to extract 19k+ genes from the originally 20k+, then only the def. expresessed genes (3k+) were selected. 
- we scaled [0 mean and one unit variance] the def. expressed genes across genes.
- the code generates 3 different files: (1) scaled Normal: samplessamplesNT_scaled.csv (2) scaled tumor samples: samplesTP_scaled.csv (3) a scaled data file combining both normal and tumor samples: All_samples_scaled.csv

 
