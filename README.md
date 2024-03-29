# clumpy
## Python pipeline to identify and extract clumps of trees from sets of phylogenetic trees with partially overlapping leaf sets

The pipeline was originally described in Serra Silva (2022) and Serra Silva and Wilkinson (202#).

### What does *clumpy* do?

*clumpy* is a Python3 pipeline that allows you to identify and extract clumps of trees from sets of trees with completely or partially overlapping leaf sets. At the moment, supertrees infereed during the pipeline, using Astral-III [(Zhang et al. 2018)](https://github.com/smirarab/ASTRAL) are not checked for polytomies (hard or otherwise), nor for ambiguously placed taxa (branch score = '?').

### "Installing" *clumpy*
You can download a zipped file or clone the repository to get the *clumpy* scripts, but you will have to install Astral-III separately.

If you want to clone the repository do:

```sh
git clone --recursive  https://github.com/anaserrasilva/clumpy
```

To run *clumpy* ensure that you have Python3 installed, as well as the following modules: 
  * [numpy](https://numpy.org/)
  * [ete3](http://etetoolkit.org/)
  * [dendropy](https://dendropy.org/)
  * [matplotlib](https://matplotlib.org/)
  
If your machine does not have graphical  capacities, you will need to download the following module to print histograms directly on the terminal:
  * [plotille](https://github.com/tammoippen/plotille)

Before trying to run *clumpy*, edit the Astral_islands.sh file to include YOUR path to Astral.

```sh
#Set the next line to Astral's path in your machine
ASTRAL=~/Astral/astral.5.7.7.jar
```

You may also need to change its permissions with:

```sh
chmod +x Astral_islands.sh
```

You should now be ready to use *clumpy*.

NOTE: If you would like to use a different software to obtain your supertrees, you need only edit the sh file to call your preferred method.

### Using *clumpy*

Example files will be provided at a later date, but the command to run *clumpy* is 

```sh
python3 phylogenomicClumpExtractor.py *.tre Distance binWidth
```

The argument *Distance* can be set to uRF (the uncorrected Robinson-Foulds distance) or wRF, which is a correction for tree size described in Serra Silva (202#). The argument *binWidth* sets the histogram bin widths, for a discussion of this setting see Serra Silva (202#). Also, please ensure your tree file is Newick formatted.

Example
```sh
python3 phylogenomicClumpExtractor.py *.tre wRF 2
```

### Potential additions/modifications

This pipeline is still *under construction* with functionalities being added.

Near future changes include a more organized output file structure, incorporation of the accessory scripts (clumpSmallTreeLocator.r and deltaInputSTclumpST.py) into the pipeline and inference/comparison of clump supertrees.

Planned modifications also include support for polytomous trees, as well as rooted trees.


### Citation

To be added


### Author contributions

To be added


### Useful readings

To be added
