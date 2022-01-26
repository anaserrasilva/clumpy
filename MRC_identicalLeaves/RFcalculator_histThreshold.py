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

    
# #collapse unsupported branches in ST
# for node in tree_collection[0].get_descendants():
    # if not node.is_leaf() and node.support == '?':
        # node.delete()
    
# dictTrees = {}
# for t in tree_collection:
    # key = tree_collection.index(t)
    # dictTrees[key] = lTrees[key]
    
    
#make dictionary equivalent of translation table, set taxon names to numbers
taxonDict = {}
counter = 1
for t in tree_collection:
    if t == tree_collection[0]:
        for leaf in t.get_leaves():
            #print(leaf.name)
            key = leaf.name
            taxonDict[key] = str(counter)
            leaf.name = str(taxonDict[leaf.name])
            counter += 1
            
    else:
        for leaf in t.get_leaves():
            if leaf.name in taxonDict.keys():
                leaf.name = taxonDict[leaf.name]
            
            else:
                #print(leaf.name)
                key = leaf.name
                taxonDict[key] = str(counter)
                leaf.name = str(taxonDict[leaf.name])
                counter += 1


####get splits
#using structured array
# create empty array, will be better than 0-filled array in multrees script

splits = np.empty((0, (len(tree_collection)+1)), dtype = object)
counter = 1
for t in tree_collection:
    
    m = []
    if len(t.get_leaves()) != len(taxonDict):
        
        names = []
        for leaf in t.get_leaves():
            names.append(leaf.name)
            
        m = list(set(list(taxonDict.values())) - set(names))
        #print(m)
                        
    for node in t.traverse():
        if node in t.get_leaves():
            pass
        
        else:
            split = []
            descend = node.get_descendants()
            for i in descend:
                if i in t.get_leaves():
                    split.append(i.name)
                
            temp = [ '.' for d in range(0, len(taxonDict))]
            
            if len(m) != 0:
                for i in m:
                    temp[int(i)-1] = '?' #partial splits
                    #need to write section to compare them
                
            for n in split:
                temp[int(n)-1] = '*'
                        
                s = str(temp).replace("'", "")
                s = s.replace(",", "")
                s = s.replace(" ", "")
                #print(s)
            
            if all(temp) == '.':
                pass
            
            elif s.count('*') == len(t.get_leaves()):
                pass
            
            elif s.count('*') == 1:
                pass
            
            elif any(splits[:,0]) == s:
                pass
            
            
            else:
                if splits.size == 0:
                    splits1 = np.zeros((1, (len(tree_collection)+1)), dtype= object)
                    #print(s)
                    splits1[0,0] = s #adds split
                    splits1[0,counter] = 1 #adds support value
                    splits = np.append(splits, splits1, axis = 0)
                
                elif s in splits[:,0]:
                    #print(s)
                    ind = str(np.where(splits[:,0] == s)) #finds index of matching split
                    ind = ind.partition('[')[2].partition(']')[0] #strips number from array
                    #print(ind)
                    splits[int(ind),counter] = 1 #adds support value
                    
                else:
                    splits1 = np.zeros((1, (len(tree_collection)+1)), dtype= object)
                    #print(s)
                    splits1[0,0] = s
                    splits1[0,counter] = 1
                    splits = np.append(splits, splits1, axis = 0)
                    
    counter += 1


#print(splits) 

###Compute shared leaves
sharedLeaves = np.zeros((1,len(tree_collection)), dtype=int)
nonsharedLeaves = [] #need for RFmatrix step
for t in tree_collection:
    
    if t == tree_collection[0]:
    
        tnames = []
        for leaf in t.get_leaves():
                #print(leaf.name)
                tnames.append(leaf.name)
        
        for p in tree_collection:
                if p == t:
                    pass
                
                elif sharedLeaves[0, tree_collection.index(p)] != 0:
                    pass
                
                else:
                    
                    pnames = []
                    for leaf in p.get_leaves():
                        #print(leaf.name)
                        pnames.append(leaf.name)
                        
                    tnamesSet = set(tnames)
                    l = len(tnamesSet.intersection(pnames))
                    sharedLeaves[0, tree_collection.index(p)] = l
                    
                    d = tnamesSet.symmetric_difference(pnames)
                    nonsharedLeaves.append(list(d))


###RF-matrix from splits

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
            
            #prevent trees without shared tips from having RF=0
            elif sharedLeaves[0, tree_collection.index(p)] == 0:
                
                RFmatrixUncorrected[0, tree_collection.index(p)] = -1
                RFmatrixNormalised[0, tree_collection.index(p)] = -1
                
            elif RFmatrixUncorrected[0, tree_collection.index(p)] != 0:
                pass
            
            else:
                
                #extract only columns of interest from splits array
                prune = splits[:, [0, (tree_collection.index(t) + 1), (tree_collection.index(p) + 1)]]
                prune = prune[np.any((prune == 1), axis=1)]
                sTable = prune[:,0]
                
                #prune splits to shared leaves
                ind = nonsharedLeaves[counter]
                ind = [ int(x)-1 for x in ind]
                irr = []
                for i, j in enumerate(sTable):
                    
                    j = j.replace('[', '')
                    j = j.replace(']', '')
                    j = j.replace(" ", "")
                    k = list(j)
                    
                    for u in sorted(ind, reverse=True):
                        del k[u]
                    
       
                    if k[0] == '*':
                        
                        if k.count('*') == len(k):
                            irr.append(i)
                        
                        k = str(k).replace("*", "#")
                        k = k.replace(".", "*")
                        k = k.replace("#", ".")
                        k = k.replace("'", "")
                        k = k.replace(",", "")
                        k = k.replace(" ", "")
                        
                        sTable[i] = k
                        
                    else: 
                        
                        if k.count('.') == len(k):
                            irr.append(i)
                        
                        k = str(k).replace("'", "")
                        k = k.replace(",", "")
                        k = k.replace(" ", "")
                        #print(k)
                
                        sTable[i] = k
                
                ##compare pruned splits and find shared pruned splits
                prune[:,0] = sTable
                prune = np.delete(prune, irr, axis=0)
                
                #remove duplicate splits from t
                pruneT = prune[prune[:,1] == 1, 0:2]
                
                pruneT = np.unique(pruneT[:, 0])
                pruneT = np.reshape(pruneT,(pruneT.size, 1))
                t1 = np.full((len(pruneT), 1), 1)
                pruneT = np.append(pruneT, t1, axis = 1)
                lT = len(pruneT)
    
                #remove duplicate splits from p
                pruneP = prune[prune[:,2] == 1, 0]
                pruneP = np.unique(pruneP)
                pruneP = np.reshape(pruneP,(pruneP.size, 1))
                lP = len(pruneP)
           
                #check for shared splits between t and p
                for r in pruneP:
                    if r in pruneT[:,0]:
                        ind = str(np.where(pruneT[:,0] == r)) #finds index of matching split
                        ind = ind.partition('[')[2].partition(']')[0] #strips number from array
                        pruneT[int(ind),1] = 2
                        
                    else:
                        prune1 = np.zeros((1, 2), dtype= object)
                        #print(s)
                        prune1[0,0] = r
                        prune1[0,counter] = 1
                        pruneT = np.append(pruneT, prune1, axis = 0)
                
                #calculate RF distance
                RF = lT + lP - (2*np.count_nonzero(pruneT[:,1] == 2))
                
                RFnorm = abs(RF/(lT+lP))
                            
                RFmatrixUncorrected[0, tree_collection.index(p)] = RF
                RFmatrixNormalised[0, tree_collection.index(p)] = RFnorm
            
        counter =+ 1
        
np.savetxt('uncorrectedRFMatrix.txt', RFmatrixUncorrected, fmt = '%d', delimiter = ',')

RFmax = np.amax(RFmatrixUncorrected)
if (RFmax % 2) != 0:
    RFmax = (RFmax + 1)

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

