from __future__ import division
import tensorflow as tf
import time
import os
from glob import glob
from six.moves import xrange
import numpy as np
import csv

class Autoencoder(object):
  def __init__(self, sess,epochs=1, learning_rate=0.001, batch_size=100, training_size=100,
   dataset="TCGA_BRCA", checkpoint_dir="checkpoint", samples_dir="samples", train=False):
   
    self.sess=sess
    self.epochs=epochs
    self.learning_rate=learning_rate
    self.batch_size=batch_size
    self.dataset=dataset
    self.checkpoint_dir=checkpoint_dir
    self.samples=samples_dir
    self.train=train
    self.file= csv.reader(open(os.getcwd()+"\\data\\"+self.dataset+".csv",'r'))
    self.data1= []
    for row in self.file:
          try:
             self.data1.append(np.float32(row))
          except ValueError:
             continue
    
    self.data=np.array(self.data1)
    print( "[*] Sequence length:",len(self.data))
    self.data=self.data.transpose()
    print("[*] Number of Samples: ",len(self.data) )
    self.training_size=len(self.data) #training_size
    self.input_length = len(self.data[1])
    self.output_length = self.input_length
    self.n_input=len(self.data[1])
    print(self.n_input)
    self.n_hidden_1=1000
    self.n_hidden_2=500
    self.n_hidden_3=250
    self.n_hidden_4=30
    #self.test_num= .2 * len(self.data)#20% of the # of samples
    # Initialization Method 
  
     
    self.weights = {

    'encoder_h1': tf.Variable(tf.random_normal([self.n_input, self.n_hidden_1])),
    'encoder_h2': tf.Variable(tf.random_normal([self.n_hidden_1, self.n_hidden_2])),
    'encoder_h3': tf.Variable(tf.random_normal([self.n_hidden_2, self.n_hidden_3])),
    'encoder_h4': tf.Variable(tf.random_normal([self.n_hidden_3, self.n_hidden_4])),
    'decoder_h1': tf.Variable(tf.random_normal([self.n_hidden_4, self.n_hidden_3])), 
    'decoder_h2': tf.Variable(tf.random_normal([self.n_hidden_3, self.n_hidden_2])),
    'decoder_h3': tf.Variable(tf.random_normal([self.n_hidden_2, self.n_hidden_1])),
    'decoder_h4': tf.Variable(tf.random_normal([self.n_hidden_1, self.n_input])),

    }
    self.biases = {
    'encoder_b1': tf.Variable(tf.random_normal([self.n_hidden_1])),
    
    
    'encoder_b2': tf.Variable(tf.random_normal([self.n_hidden_2])),
    'encoder_b3': tf.Variable(tf.random_normal([self.n_hidden_3])),
    'encoder_b4': tf.Variable(tf.random_normal([self.n_hidden_4])),
    'decoder_b1': tf.Variable(tf.random_normal([self.n_hidden_3])),
    'decoder_b2': tf.Variable(tf.random_normal([self.n_hidden_2])),
    'decoder_b3': tf.Variable(tf.random_normal([self.n_hidden_1])),
    'decoder_b4': tf.Variable(tf.random_normal([self.n_input])),
    }
    
    
    self.build_model()
  def build_model(self):

        self.training_samples = tf.placeholder(tf.float32, [self.batch_size], name='training_samples')
        #self.testing_samples = tf.placeholder(tf.float32, [self.test_num], name='test_outputs')

        training_samples=self.training_samples
        #testing_samples= self.testing_samples
        self.saver= tf.train.Saver()
        if (self.train):
            
            self.training(training_samples)
            
        else:
            self.test()
            '''
            self.test(testing_samples)
            '''
            '''
            self.e_test=self.encoder(testing_samples)
            self.d_test=self.decoder(self.e_test)
            '''
  def training(self, x):
        X=tf.placeholder('float', [None,self.n_input])
        self.e= self.encoder(X)
        self.d= self.decoder(self.e)
            
        cost=tf.reduce_mean(tf.pow(self.d-X,2))
        optimizer = tf.train.RMSPropOptimizer(self.learning_rate).minimize(cost)
        # Initializing the variables
        try:
            tf.global_variables_initializer().run()
        except:
            tf.initialize_all_variables().run()
        counter=1
        display_step =2
        start_time= time.time()
    
        could_load = self.load(self.checkpoint_dir)
        if could_load:
           
           print(" [*] Load SUCCESS")
        else:
           print(" [!] Load failed...")
        
        for epoch in xrange(self.epochs):
            
            batch_idxs = min(len(self.data),self.training_size) // self.batch_size
            
            #start trainning
            for idx in xrange(batch_idxs):
                
                batch_samples = self.data[idx*self.batch_size:(idx+1)*self.batch_size]#np.array(batch).astype(np.float32)
               
                _, c = self.sess.run([optimizer, cost], feed_dict={X: self.mask_noise(np.reshape(batch_samples,[self.batch_size,self.n_input]),20)})

            if epoch % display_step == 0:
                   print("Epoch:", '%04d' % (epoch+1),"cost=", "{:.9f}".format(c))
            if np.mod(counter, 2) == 0:
                   self.save(self.checkpoint_dir, counter)
            counter+=1
  def encoder(self,x):
        '''
        with tf.variable_scope("encoder") as scope:
             #encoder weigths
             encoder_h1= tf.Variable(tf.random_normal([n_input, n_hidden_1]))
             encoder_h2= tf.Variable(tf.random_normal([n_hidden_1, n_hidden_2]))
             encoder_h3= tf.Variable(tf.random_normal([n_hidden_2, n_hidden_3]))
             encoder_h4= tf.Variable(tf.random_normal([n_hidden_3, n_hidden_4]))
             #encoder biases
             encoder_b1= tf.Variable(tf.random_normal([n_hidden_1]))
             encoder_b2= tf.Variable(tf.random_normal([n_hidden_2]))
             encoder_b3= tf.Variable(tf.random_normal([n_hidden_3]))
             encoder_b4= tf.Variable(tf.random_normal([n_hidden_4]))
        '''
        layer_1 = tf.nn.sigmoid(tf.add(tf.matmul(x, self.weights['encoder_h1']),
                                   self.biases['encoder_b1']))
              
        layer_2 = tf.nn.sigmoid(tf.add(tf.matmul(layer_1, self.weights['encoder_h2']),
                                   self.biases['encoder_b2']))
        layer_3 = tf.nn.sigmoid(tf.add(tf.matmul(layer_2, self.weights['encoder_h3']),
                                   self.biases['encoder_b3']))
        layer_4 = tf.nn.sigmoid(tf.add(tf.matmul(layer_3, self.weights['encoder_h4']),
                                   self.biases['encoder_b4']))
    
        return layer_4
  def decoder(self,x):
        '''
        with variable_scope("decoder") as scope:
             #decoder weights
             decoder_w1=tf.Variable(tf.random_normal([n_hidden_4, n_hidden_3]))
             decoder_w2= tf.Variable(tf.random_normal([n_hidden_3, n_hidden_2]))
             decoder_w3= tf.Variable(tf.random_normal([n_hidden_2, n_hidden_1]))
             decoder_w4= tf.Variable(tf.random_normal([n_hidden_1, n_input]))
             #decoder biases
             decoder_b1= tf.Variable(tf.random_normal([n_hidden_3]))
             decoder_b2= tf.Variable(tf.random_normal([n_hidden_2]))
             decoder_b3= tf.Variable(tf.random_normal([n_hidden_1]))
             decoder_b4= tf.Variable(tf.random_normal([n_input]))
        '''
        layer_1 = tf.nn.sigmoid(tf.add(tf.matmul(x, self.weights['decoder_h1']),
                                   self.biases['decoder_b1']))
    
        layer_2 = tf.nn.sigmoid(tf.add(tf.matmul(layer_1, self.weights['decoder_h2']),
                                   self.biases['decoder_b2']))
        layer_3 = tf.nn.sigmoid(tf.add(tf.matmul(layer_2, self.weights['decoder_h3']),
                                   self.biases['decoder_b3']))
        layer_4 = tf.nn.sigmoid(tf.add(tf.matmul(layer_3, self.weights['decoder_h4']),
                                   self.biases['decoder_b4']))
    
        return layer_4
  def test(self):
        #saver.restore(sess,'./modelMultiDenoisedBad.ckpt')
        X=tf.placeholder('float', [None,self.n_input])
        self.e= self.encoder(X)
        self.d= self.decoder(self.e)
            
        cost=tf.reduce_mean(tf.pow(self.d-X,2),1)
        could_load = self.load(self.checkpoint_dir)
        if could_load:
           
           print(" [*] Load SUCCESS")
        else:
           print(" [!] Load failed...")
        
        test_samples = self.data#np.array(batch).astype(np.float32)
                #print("[***] patch samples shape = ",batch_samples.shape)
        encode_decode = self.sess.run(cost, feed_dict={X: np.reshape(self.data,[len(self.data),self.n_input])})
        print(encode_decode)
        np.savetxt("{}_{}".format(self.dataset,"cost.csv"), encode_decode, delimiter=",")
        reps = self.sess.run(self.e, feed_dict={X: np.reshape(self.data,[len(self.data),self.n_input])})
        np.savetxt("{}_{}".format(self.dataset,"reps.csv"), np.transpose(reps, (0,1)), delimiter=",")

        '''
        images_list=mnist.test.images[:500]
        x=np.reshape(x,(1,784))
        # Applying encode and decode over test set
        encode_decode = sess.run(
        y_pred, feed_dict={X: x}) #  examples_to_show
        np.savetxt("NonoptimizedImagesMulti10Examples.csv", encode_decode, delimiter=",")
        tl.visualize.W(W=np.transpose(collection, (1,0)), second=10, saveable=True, shape=[28,28], name='OptimizeMulti10Examples', fig_idx=2396512)
        '''
  def xavier_init(self,nin,nout,const=1):
     low = -const * np.sqrt(1/ (nin + nout))
     high = const * np.sqrt(1/ (nin + nout))
    
     return tf.random_uniform((nin, nout), minval=low, maxval=high)
    # Noising Method
  def mask_noise(self,x,v):
    x_noise = x.copy()

    n_samples = x.shape[0]
    n_features = x.shape[1]

    for i in range(n_samples):
        mask = np.random.randint(0, n_features, v)

        for m in mask:
            x_noise[i][m] = 0.

    return x_noise
  @property
  def model_dir(self):
        
           return "{}_{}_{}_{}".format(
           self.dataset, self.batch_size, self.input_length, self.output_length)
        
  def save(self, checkpoint_dir, step):
        model_name = "Autoencoder.model"
        checkpoint_dir = os.path.join(checkpoint_dir, self.model_dir)

        if not os.path.exists(checkpoint_dir):
          os.makedirs(checkpoint_dir)

        self.saver.save(self.sess,
            os.path.join(checkpoint_dir, model_name))

  def load(self, checkpoint_dir):
        import re
        print(" [*] Reading checkpoints...")
        if self.train: 
           checkpoint_dir = os.path.join(checkpoint_dir, self.model_dir)
        else:
           
            i=self.dataset.find("test")
            if (i>0):
                name= self.dataset[0:i-1] + "_train" + self.dataset[i+4:]
                model_dir= "{}_{}_{}_{}".format(
                name,self.batch_size,self.input_length,self.output_length)
                            
                print("[************] Model_dir = "+ model_dir)
                checkpoint_dir = os.path.join(checkpoint_dir, model_dir)
            else:
                checkpoint_dir = os.path.join(checkpoint_dir, self.model_dir)
        ckpt = tf.train.get_checkpoint_state(checkpoint_dir)
        if ckpt and ckpt.model_checkpoint_path:
          ckpt_name = os.path.basename(ckpt.model_checkpoint_path)
          self.saver.restore(self.sess, os.path.join(checkpoint_dir, ckpt_name))
          
          print(" [*] Success to read {}".format(ckpt_name))
          return True
        else:
          print(" [*] Failed to find a checkpoint")
          return False

