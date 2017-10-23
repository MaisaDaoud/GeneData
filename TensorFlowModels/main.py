import os
import numpy as np
import tensorflow as tf
import pprint
from Autoencoder import Autoencoder

flags= tf.app.flags
flags.DEFINE_integer("epochs", 1, "number of training epochs")
flags.DEFINE_float("learning_rate",0.001,"learning rate")
flags.DEFINE_integer("batch_size",100,"number of samples/batch")
flags.DEFINE_integer("training_size",100,"number of training samples") 
#flags.DEFINE_integer("input_length",1000," gene expression (input) length") 
#flags.DEFINE_integer("output_length",1000,"output length")
flags.DEFINE_string("dataset","TCGA_BRCA","dataset name")
flags.DEFINE_string("extension",".csv", "image extension")
flags.DEFINE_string("checkpoint_dir","checkpoint","directory to save the models")
flags.DEFINE_string("sample_dir","samples","directory to save a sample of test images")
flags.DEFINE_boolean("train", False, "true to train , false to test the model")
FLAGS= flags.FLAGS

def main(_):
 #... printing the flags ....
 pp=pprint.PrettyPrinter()
 pp.pprint(flags.FLAGS.__flags)
 
 
 # Definning the model based on the dataset
 with tf.Session() as sess:
  
      auto=Autoencoder(
        sess,
        epochs=FLAGS.epochs,
        learning_rate=FLAGS.learning_rate,
        batch_size=FLAGS.batch_size,
        training_size=FLAGS.training_size,
        #input_height= FLAGS.input_length,
        #output_height=FLAGS.output_length,
        dataset=FLAGS.dataset,
        checkpoint_dir=FLAGS.checkpoint_dir,
        samples_dir=FLAGS.sample_dir,
        train=FLAGS.train)
   
 
 
if __name__=="__main__":
 tf.app.run()