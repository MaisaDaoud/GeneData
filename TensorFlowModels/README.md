# TensorFlowModels
  -This project implements stacked denoising autoencoder (SDA) that reads data as (*.csv) file located in (data) folder and save the        trained model in (checkpoint) folder
###### Use the following command to train  the model:
      python main.py --train --dataset="TCGA_BRCA_NT_scaled_train"  --epochs=1000 --batch_size=25
      
###### Use the following command to test the model:
      python main.py  --dataset="TCGA_BRCA_NT_scaled_test"  

