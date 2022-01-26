#!/bin/bash

#What this script does...

#Set the next line to Astral's path in your machine
ASTRAL=~/Astral/astral.5.7.7.jar

#make ST
java -jar $ASTRAL -i inputTrees.tre -o supertree.tre

#combine ST and input trees into single file
cat supertree.tre inputTrees.tre > compareTrees.tre

