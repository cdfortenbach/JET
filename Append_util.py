# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 17:47:50 2019

@author: Charles Fortenbach
"""

import numpy as np


###############################################################################
# Merge hi number Pdxo_Output file with low number Pdxo_Output file
###############################################################################

# Must run this within the Linux file structure

# From terminal within low number JET directory, 
# enter cmd: python Append_util.py

dataset_low = np.loadtxt('Pdxo_Output.txt', ndmin=2, skiprows=0)
targets_dataset_low = np.max(dataset_low, axis=0)
last_target_low = targets_dataset_low[2]

PATH = input(' Enter the directory path to the hi number target Pdxo_Output file: ')
# for example:  '/media/Win7share/JET_wkg_2/'

dataset_hi = np.loadtxt(PATH + 'Pdxo_Output.txt', ndmin=2, skiprows=0)
targets_dataset_hi = np.min(dataset_hi, axis=0)
first_target_hi = targets_dataset_hi[2]

# check compatibility then append data
if first_target_hi == last_target_low +1:
    with open('Pdxo_Output.txt', 'ba') as f: 
        np.savetxt(f, dataset_hi)
        print('Pdxo_Output file updated')
        
else: print('mismatch in Pdxo_Output files')        