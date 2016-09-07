#!/usr/bin/env python
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Check force and velocity directions in CSV files extracted from data

# Set filename (Replace with a dir and for loop later)
folder = '/home/raeed/Dropbox/Research/ForceKin/ForceKin Paper/Data/'
for fname in os.listdir(folder):
    if fname.endswith('.csv'):
        # Open CSV file and read first two columns
        thv = pd.read_csv(folder+fname,usecols=[0])
        thf = pd.read_csv(folder+fname,usecols=[1])
        
        # scatter plot of points
        plt.figure(figsize=(12,9))
        plt.xlim(-np.pi,np.pi)
        plt.ylim(-np.pi,np.pi)
        plt.scatter(thf,thv)
        plt.gca().set_aspect('equal', adjustable='box')
        
        # Save figure
        plt.savefig(folder + 'domains/default/' + os.path.splitext(fname)[0] + '.png',bbox_inches='tight')
        continue
    else:
        continue


