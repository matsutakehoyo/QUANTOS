# QUANTOS
Scripts for QUANTOS (QUalitative and quantitative ANalysis using Bayes Theorem Optimized for Synapse evaluation)

## Contents
QUANTOS scripts and sample IHC images are conteined in the "Scripts_SAmpleData" folder. Once you successfully run the scripts following the procedure described below, you should get exactly the same files as in "Results" folder.

### System Requirements
Scripts were confirmed to run under following environment:
* macOS Hight Sierra (version 10.13.3)
* ImageJ (version 1.51)
* ImageJ Fiji (version)
* R (version 3.4.3)
* RStudio (Version 1.1.414)

### Software preparation
Please install R, RStudio, ImageJ, ImageJ Fiji softwares before running the scripts.
Also, please make sure that following plug-ins and packages are installed:

* ImageJ Fiji plug-ins: [Adjustable watershed](http://imagejdocu.tudor.lu/doku.php?id=plugin:segmentation:adjustable_watershed:start), 
[RATS](https://imagej.net/RATS:_Robust_Automatic_Threshold_Selection)
* R Studio packages: MASS, grid,gridExtra,ks, tidyverse, wesanderson, propagate, bde

### Sample preparation
Image file for analysis should be TIFF format with 1024 by 1024 pixel resolution, and have 3 channels with DAPI, pre-synaptic marker, post-synaptic marker respectively. This analysis requires 5 sequential z-stack image.
A sample IHC Image is in "Scripts_SampleData" folder.


### Procedure
1. Open ImageJ Fiji
2. run "1_ImageProcessing.ijm"
3. Select the directory of image file to analyze. Once the image processing is done, csv files with graphical information of ROIs will be created. 
4. Open Studio
5. Run "2_data.prep.fiji.R". This script gathers csv file information and export "data2D.RDa" that contains necessary data for R analysis.
6. Run "3_find_pair". This script identifys pre- and post-synaptic markers pairs within 1.2 um, and export "data_pair.rds".
7. Run "4_estimate_bg.R" This script gathers background intensity information and exports "data_pair_bg.rds" and "data_bg.rds"
8. Run "5_find_synapse.R" This script references probability density functions (PDFs) of Ideal Synapse and Ideal Noise, to estimate the posterior probability of synapse candidates. Once executed, following files will be exported:
	log : number of detected synapses are exported  in log subdirectory. "Number of synapses (Bayes)" is the final result of the analysis
	synapse.csv : center coordinates of detected synapses
	synapse_coord.pdf : plot of synapse coordinates.
	marker_coord.pdf : plot of all marker coordinates.
	ROI sats probability distribution (parameters).pdf : Histogram of each parameters of sample image overlaid with Ideal Synapse and Ideal Noise distribution.
	Naive Bayes scatter plot.pdf : scatter plot of likelihood being synapse and noise.
9. Run "6_find_synapse_maturation.R" This script analyzes the maturation of detected synapses. Mature synapse numbers will be exported in "log_mature" directory. Similar files are exported as find_synapse.R script, except they have "maturation" at the bottom of the file names.
made by Ryutaro Akiba and Take Matsuyama
