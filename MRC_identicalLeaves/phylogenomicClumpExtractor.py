'''
clumpExtractor_histogramGPU

Usage: python3 phylogenomicClumpExtractor.py *.tre

    
     splitSupport:   Reads in a collection of trees in Newick format and extracts the tree clumps identifiable under the desired threshold.
                     Nexus files are not supported at this time.
                     
                     To use Astral as the supertree calculator, edit the path to Astral in the Astral_islands.sh file.
                     To ensure the script will run in your system use the command "chmod +x Astral_islands.sh" before the first use.
                     
                     The first argument corresponds to the file of Newick formatted trees.
                     
                     After the Robinson-Foulds distance of all trees to the supertree is calculated, a window showing the histogram of the distances will open.
                     After selecting a threshold type it into the console.
                     
------------------------------------
For Python 3.8
Python modules required:
        - ete3
        -Numpy
        -subprocess
        -shutil
        -matplotlib
		
Dependencies: Astral

R packages required:
		- ape
------------------------------------
TO BE NAMED
https://github.com/anaserrasilva/TOBENAMED
Written by Ana Serra Silva
Ideas by Ana Serra Silva and Mark Wilkinson
a.da-silva@nhm.ac.uk; m.wilkinson@nhm.ac.uk
ADD DATE
Distributed under the 
MIT License
'''
#import modules
from sys import argv
import os
import subprocess
import shutil
import glob
import numpy as np
import ete3 as ete
import dendropy as dendro


#set arguments
Script = argv[0]
Infile = argv[1]


#copy input file
shutil.copy(Infile, 'inputTrees.tre')

#checks length of input file
#lt = [l.strip() for l in open("inputTrees.tre")]

#island extractor
while len([l.strip() for l in open("inputTrees.tre")]) != 0: #check input file is not empty
    
    #run Astral and concatenate ST and input file
    #subprocess.check_output(["./Astral_islands.sh"])
    #compute M50 consensus
    trees = dendro.TreeList.get(path='inputTrees.tre', schema='newick')
    conTree = trees.consensus(min_freq=0.5)
    conString = conTree.as_string(schema='newick').replace('[&U] ', '')
    
    #save M50    
    temp = glob.glob('consensus*')
    tempName = ["consensus_", str(len(temp)+1), ".tre"]
    fileName = ''.join(tempName)
    
    f = open(fileName, 'w')
    f.write(conString)
    f.close()
    
    f1 = open('consensusTree.tre', 'w')
    f1.write(conString)
    f1.close()
    
    #conTree.write(path = fileName, schema = 'newick',)
    #conTree.write(path = 'consensusTree.tre', schema = 'newick',)
    
    #concatenate tree files
    #subprocess.run(["cat consensusTree.tre inputTrees.tre > compareTrees.tre", shell=True])
    filenames = ['consensusTree.tre', 'inputTrees.tre']
    with open('compareTrees.tre', 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
   
    os.remove('consensusTree.tre')
    
    
    #run RFcalculator.py, includes pre-analysis collapse of unsupported branches (localPP=?)
    #also IDs minimum threshold
    subprocess.call(["python3", "-u", "RFcalculator_histThreshold_noSplits.py"])
    
    #run islandExtractor.R
    subprocess.call(["python3","-u", "clumpExtractor.py"])


#clean temp files
os.remove('inputTrees.tre')
os.remove('compareTrees.tre')
os.remove('uncorrectedRFMatrix.txt')
os.remove('threshold.txt')
