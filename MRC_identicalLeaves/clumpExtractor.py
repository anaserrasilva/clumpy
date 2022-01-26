# -*- coding: utf-8 -*-
"""
Created on Thu May 20 10:37:22 2021

@author: anad4
"""

import dendropy as dendro
import numpy as np
#import os
import glob
#import shutil
import re

#read in trees and distance matrix
trees = dendro.TreeList.get(path='compareTrees.tre', schema='newick')
#tempST = glob.glob('supertree*')
#tempSTname = ["supertree_", str(len(tempST)), ".tre"]
#fileNameST = ''.join(tempSTname)
#trees[0].write(path = fileNameST, schema = 'newick',)


RFmatrix = np.loadtxt('uncorrectedRFMatrix.txt', delimiter=',', dtype=int)

threshold = int(np.loadtxt('threshold.txt', dtype=int))
temp = glob.glob('threshold*')
tempName = ["threshold_", str(len(temp)), ".txt"]
fileName = ''.join(tempName)
#save threshold file matching clump number
f = open(fileName, "w")
f.write(str(threshold))
f.close()


#clump extraction
x = threshold
l = len(trees)
mat = RFmatrix

#adding property, equivalent of colour in graph-based clustering approaches
p = np.full((1, l), 'a')
p[0,0] = 'b'
s = np.where(RFmatrix <= x)
p[0,s] = 'b'

#find indexes of trees to extract
r = np.where(p[0,:] == 'b')
r = r[0]
r = r.tolist()

#extract clump
clump = dendro.TreeList()
for index in r:
    if index != 0:
        clump.append(trees[index])

tempClump = glob.glob('clump_*')
tempClumpName = ["clump_", str(len(tempClump)+1), ".tre"]
fileNameClump = ''.join(tempClumpName)
clump.write(path = fileNameClump, schema = 'newick',)

#edit inputTrees file
inputTrees = dendro.TreeList.get(path='inputTrees.tre', schema='newick')

r.remove(r[0])
r = [x - 1 for x in r]
if len(inputTrees) >= len(clump):
    
    # if len(inputTrees) == (len(trees)-1):
    #     print('No more trees to analyse.')
    #     f = open("inputTrees.tre", "w")
    #     f.close()
        
    #     exit()
        
    # else:
    counter = 0
    for index in r:
        
        inputTrees.remove(inputTrees[int(index - counter)])
        counter += 1    


if len(inputTrees) != 0:
    inputTrees.write(path = 'inputTrees.tre', schema = 'newick')
    
if len(inputTrees) == 0:
    print('No more trees to analyse.')
    f = open("inputTrees.tre", "w")
    f.close()


