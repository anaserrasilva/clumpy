# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 17:55:45 2021

@author: anad4
"""

import ete3 as ete
import numpy as np
import os
import glob
import re


##Set working directory
os.chdir("C:\\Users\\anad4\\Dropbox\\clumpExtractors\\resultsWeightedRF\\Chan\\legacy\\wRF\\legacy_break")

##Locate and read in input supertrees
temp_iST = glob.glob('supertrees/supertree_*')

#sort list of files into proper numerical order of trees
sorted_iST = [None] * len(temp_iST)
for i in temp_iST:
    ind = int(re.findall(r'\d+', i)[0])
    sorted_iST[ind-1] = ete.Tree(i)

for t in sorted_iST:
    t.unroot()

##Locate and read in clump supertrees
temp_cST = glob.glob('clumpSupertrees/clump_*')

#sort list of files into proper numerical order of trees
sorted_cST = [None] * len(temp_cST)
for c in temp_cST:
    ind = int(re.findall(r'\d+', c)[0])
    sorted_cST[ind-1] = ete.Tree(c)

for t in sorted_cST:
    t.unroot()

##Calculate RF between iST and cST
#distance between iST and cST pair
RFdist_iST_cST = [None] * len(sorted_iST)
for index_i, phylo_i in enumerate(sorted_iST):
    for index_c, phylo_c in enumerate(sorted_cST):
        if index_i == index_c:
            RFdist_iST_cST[index_i] = phylo_i.robinson_foulds(phylo_c, unrooted_trees=True)[0]
            

#distance between iST[0] and each cST
RFdist_iSTfirst_cST = [None] * len(sorted_iST)
for index_c, phylo_c in enumerate(sorted_cST):
    RFdist_iSTfirst_cST[index_c] = sorted_iST[0].robinson_foulds(phylo_c, unrooted_trees=True)[0]

#save tables
with open('RFdist_inputST_clumpSTPair.txt', 'w') as file1:
    for item in RFdist_iST_cST:
        file1.write("%s\n" % item)

with open('RFdist_firstST_clumpSTs.txt', 'w') as file2:
    for item in RFdist_iSTfirst_cST:
        file2.write("%s\n" % item)

