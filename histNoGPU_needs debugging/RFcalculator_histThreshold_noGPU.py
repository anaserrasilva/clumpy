# -*- coding: utf-8 -*-
"""
Helper script for ...
"""

from sys import argv
import ete3 as ete
import numpy as np
#import os
#import matplotlib.pyplot as plt
import plotille as tille
import glob
import math


#set arguments
Script = argv[0]
Distance = argv[1]
binWidth = argv[2]

binWidth = int(binWidth)

######load trees and prep them analysis
tree_collection = []
for newick in open('compareTrees.tre'):
    tree_collection.append(ete.Tree(newick))

for t in tree_collection:
    t.unroot()
	
    
# #collapse unsupported branches in ST
# for node in tree_collection[0].get_descendants():
    # if not node.is_leaf() and node.support == '?':
        # node.delete()
    

###RF-matrix from ete

RFmatrixUncorrected = np.zeros((len(tree_collection),), dtype=object)
#RFmatrixNormalised = np.zeros((1,len(tree_collection)), dtype=object)
RFmatrixWeighted = np.zeros((len(tree_collection),), dtype=object)
counter = 0
for t in tree_collection:
    
    if t == tree_collection[0]:
    
        for p in tree_collection:
            
            if t == p:
                pass
            
            elif tree_collection.index(t) > tree_collection.index(p):
                pass
            
            else:
                #calculate RF distance
                RF = tree_collection[0].robinson_foulds(p, unrooted_trees=True)
                
                #RFnorm = abs(RF[0]/(len(tree_collection[0].get_edges())))
                            
                wRF = RF[0] * (len(tree_collection[0].get_edges())/((len(RF[3]) - len(RF[5]))))
                
                RFmatrixUncorrected[tree_collection.index(p),] = RF[0]
                #RFmatrixNormalised[0, tree_collection.index(p)] = RFnorm
                RFmatrixWeighted[tree_collection.index(p),] = wRF
            
        counter =+ 1

if Distance == 'uRF':
    
    temp = glob.glob('uncorrectedRFMatrix*')
    tempName = ["uncorrectedRFMatrix_", str(len(temp)+1), ".txt"]
    fileName = ''.join(tempName)        
    np.savetxt(fileName, RFmatrixUncorrected, fmt = '%f', delimiter = ',')
    
    RFmax = math.ceil(np.amax(RFmatrixUncorrected))
    if (RFmax % binWidth) != 0:
        RFmax = (RFmax + (binWidth - (RFmax % binWidth)))
        
    if RFmax == 0:
        RFmax =+ binWidth
        
    RFmatrix = np.copy(RFmatrixUncorrected)
    
if Distance == 'wRF':
    temp = glob.glob('weightedRFMatrix*')
    tempName = ["weightedRFMatrix_", str(len(temp)+1), ".txt"]
    fileName = ''.join(tempName)        
    np.savetxt(fileName, RFmatrixWeighted, fmt = '%.2f', delimiter = ',')
    
    RFmax = math.ceil(np.amax(RFmatrixWeighted))
    if (RFmax % binWidth) != 0:
        RFmax = (RFmax + (binWidth - (RFmax % binWidth)))
        
    if RFmax == 0:
        RFmax =+ binWidth
        
    RFmatrix = np.copy(RFmatrixWeighted)

np.savetxt('RFmatrix.txt', RFmatrix, fmt = '%.2f', delimiter = ',')

# binList = []
# for i in range(0, (RFmax+binWidth), binWidth):
    # binList.append(i)    

print(np.delete(RFmatrix,0).tolist(), width=int(binWidth))

threshold = input("Enter threshold value: ")
f = open("threshold.txt", "w")
f.write(threshold)
f.close()

