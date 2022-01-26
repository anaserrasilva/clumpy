# -*- coding: utf-8 -*-
"""
Helper script for ...
"""

import ete3 as ete
import numpy as np
#import os
import matplotlib.pyplot as plt
#import plotille as tille
import glob


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

RFmatrixUncorrected = np.zeros((1,len(tree_collection)), dtype=object)
RFmatrixNormalised = np.zeros((1,len(tree_collection)), dtype=object)
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
                RF = tree_collection[0].robinson_foulds(p, unrooted_trees=True)[0]
                
                RFnorm = abs(RF/(tree_collection[0].robinson_foulds(t, unrooted_trees=True, polytomy_size_limit=len(t.get_leaves()))[1]))
                            
                RFmatrixUncorrected[0, tree_collection.index(p)] = RF
                RFmatrixNormalised[0, tree_collection.index(p)] = RFnorm
            
        counter =+ 1
        
np.savetxt('uncorrectedRFMatrix.txt', RFmatrixUncorrected, fmt = '%d', delimiter = ',')

RFmax = np.amax(RFmatrixUncorrected)
if (RFmax % 2) != 0:
    RFmax = (RFmax + 1)

if RFmax == 0:
    RFmax =+ 2

#plt.hist(RFmatrixUncorrected.tolist(), bins=100)
print("Close plot window to enter threshold value.")
plt.hist(RFmatrixUncorrected.tolist(), histtype = 'stepfilled', bins=int(RFmax/2))
temp = glob.glob('histogram*')
tempName = ["histogram_", str(len(temp)+1), ".svg"]
fileName = ''.join(tempName)
plt.savefig(fileName)
plt.show()

#print(tille.hist(RFmatrixUncorrected[0,].tolist(), bins=110))

threshold = input("Enter threshold value: ")
f = open("threshold.txt", "w")
f.write(threshold)
f.close()

