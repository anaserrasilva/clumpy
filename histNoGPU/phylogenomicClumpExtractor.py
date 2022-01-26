'''
clumpExtractor_histogramGPU

Usage: python3 phylogenomicClumpExtractor.py *.tre Distance binWidth

    
     splitSupport:   Reads in a collection of trees in Newick format and extracts the tree clumps identifiable under the desired threshold.
                     Nexus files are not supported at this time.
                     
                     To use Astral as the supertree calculator, edit the path to Astral in the Astral_islands.sh file.
                     To ensure the script will run in your system use the command "chmod +x Astral_islands.sh" before the first use.
                     
                     The first argument corresponds to the file of Newick formatted trees.
                     
                     The second argument, Distance, corresponds to your Robinson-Foulds distance of choice, 'uncorrectedRF' corresponds to the raw RF distance between two partially overlapping trees and 'weightedRF' corresponds to the RF modification described in Serra Silva and Wilkinson (202X).
                     
                     The third argument, binWidth, sets the histogram bin width.
                     
                     After the Robinson-Foulds distance of all trees to the supertree is calculated, a window showing the histogram of the distances will open.
                     After selecting a threshold type it into the console.
                     
------------------------------------
For Python 3.8
Python modules required:
        -ete3
        -Numpy
        -subprocess
        -shutil
        -plotille
		-glob
		-math
		
Dependencies: Astral

R packages required:
		- ape
------------------------------------
Clumpy
https://github.com/anaserrasilva/clumpy
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
import numpy as np
import ete3 as ete


#set arguments
Script = argv[0]
Infile = argv[1]
Distance = argv[2]
binWidth = argv[3]


#copy input file
shutil.copy(Infile, 'inputTrees.tre')

#checks length of input file
#lt = [l.strip() for l in open("inputTrees.tre")]

#island extractor
while len([l.strip() for l in open("inputTrees.tre")]) != 0: #check input file is not empty
    
    #run Astral and concatenate ST and input file
    subprocess.call(["./Astral_islands.sh"])
    
    #run RFcalculator.py, includes optional pre-analysis collapse of unsupported branches (localPP=?), commented out
    #also IDs minimum threshold
    subprocess.call(["python3", "-u", "RFcalculator_histThreshold_noGPU.py", Distance, binWidth])
    
    #run clumpExtractor.py
    subprocess.call(["python3", "-u", "clumpExtractor.py"])
	
	if len(glob.glob("supertree_*")) != len(glob.glob("histogram_*")):
        print("The RF calculation failed. Please check all arguments present. If so, allocate more memory to the supertree building step.")
        exit()


#clean temp files
os.remove('supertree.tre')
os.remove('inputTrees.tre')
os.remove('compareTrees.tre')
os.remove('uncorrectedRFMatrix.txt')
os.remove('threshold.txt')