#!/usr/bin/env python
# coding: utf-8

# In[1]:



#get_ipython().run_line_magic('load_ext', 'autoreload')
#get_ipython().run_line_magic('autoreload', '2')
import numpy as np
import bardensr
import tensorflow as tf
import time
import numpy.random as npr
import matplotlib.pyplot as plt
import bardensr.plotting
import tqdm.notebook
import h5py
import scipy as sp
import scipy.io
import time
import imageio
import pandas as pd

import glob


#get_ipython().run_line_magic('matplotlib', 'inline')
basefn=''


# In[2]:


# control slice indices and target fdr threshold
controlidx=range(25,35) #control slice numbers
fdrthresh=0.05 #set target fdr threshold
trim=160


# 
# # get code for loading data

# In[3]:


# gene names
matfile = scipy.io.loadmat(basefn+'codebookLC_base.mat')
gene_names = np.array([str(a[0]) for a in (matfile['codebook'][:, 0])])

# codebook
codebook=scipy.io.loadmat(basefn+'codebookforbardensr.mat')['codebookbin1'].astype(float)
codebook=codebook.transpose([2,1,0])
R,C,J=codebook.shape

# flattened version, flattening over all 28 frames
codeflat=np.reshape(codebook,(R*C,J))

# get list of codes which are for error correction (ERror Correction Codes :) )
ercc_codes=np.array([x.startswith('unused') for x in gene_names])


# In[4]:


def grab_xiaoyin_data(fn,trim):
#use this for basecalling the actual image,but also allows trimming. The trimming should be fairly small (usually starting at 5%) 
    D=[]
    for i in range(1,R+1):
        with imageio.get_reader(fn+'%02d.tif'%i) as f:
            for c in range(C):
                D.append(f.get_data(index=c)[None])
    D1=np.array(D)
    ydim=np.shape(D1)[2]
    xdim=np.shape(D1)[3]
    return D1[:,:,trim:-trim,trim:-trim]  # <-- 28 x 1 x 3200 x 3200


# In[5]:


def grab_xiaoyin_data_center(fn):
#use this for calculating threshold from the center of the images
    cropfactor2=0.4 #this should be <0.5
    D=[]
    for i in range(1,R+1):
        with imageio.get_reader(fn+'%02d.tif'%i) as f:
            for c in range(C):
                D.append(f.get_data(index=c)[None])
    D1=np.array(D)
    ydim=np.shape(D1)[2]
    xdim=np.shape(D1)[3]
    return D1[:,:,round(ydim*cropfactor2):round(ydim*(1-cropfactor2)),round(xdim*cropfactor2):round(xdim*(1-cropfactor2))]  # <-- 28 x 1 x 3200 x 3200


# In[6]:


folderlist=glob.glob('processed/MAX*/aligned/')
folderlistcontrol=[folderlist[i] for i in controlidx]


# # get normalizing constants (using all fovs) and perform normalization
# 

# In[7]:


maxs=np.array([grab_xiaoyin_data_center(fn+'alignedfixedn2vgeneseq').max(axis=(1,2,3)) for fn in folderlistcontrol])
maxmax=np.median(maxs,axis=0) #max of each channel/cycle across fovs. Used to normalize all images. 


# read predefined threshold.

with open('thresh_refined.txt') as f:
    thresh_refined=float(f.readline())

# # spot call each fov, using the thresholds we decided on, and the normalization we decided on

# In[14]:

noisefloor=0.05

for fn in folderlist:
    et=bardensr.spot_calling.estimate_density_singleshot(
        grab_xiaoyin_data(fn+'alignedfixedn2vgeneseq',trim)/maxmax[:,None,None,None],
        codeflat,
        noisefloor
    )

    spots=bardensr.spot_calling.find_peaks(et,thresh_refined,use_tqdm_notebook=True)
    spots.loc[:,'m1']=spots.loc[:,'m1']+trim
    spots.loc[:,'m2']=spots.loc[:,'m2']+trim
    spots.to_csv(fn+'bardensrresult.csv',index=False)
    print(f"in {fn} found {len(spots)} spots")
    

