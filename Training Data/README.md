Contents
QUANTOS scripts for generating training data are contained in this "Scripts_Trainingdata" folder.
In "Scripts_Sampledata" folder, pdf used in QUANTOS manuscript is already there.
You only need to run this script when you want to generate PDFs from different training data.

System Requirements
Scripts were confirmed to run under following environment:
macOS Hight Sierra (version 10.13.3)
ImageJ (version 1.51)
ImageJ Fiji (version)
R (version 3.4.3)
RStudio (Version 1.1.414)

Software preparation
Please install R, RStudio, ImageJ, ImageJ Fiji softwares before running the scripts.
Also, please make sure that following plug-ins and packages are installed:

[ImageJ Fiji plug-ins]
	1, Adjustable watershed
	http://imagejdocu.tudor.lu/doku.php?id=plugin:segmentation:adjustable_watershed:start
	2, RATS
	https://imagej.net/RATS:_Robust_Automatic_Threshold_Selection

[R Studio packages]
	1, MASS
	2, grid
	3, gridExtra
	4, ks
	5, plyr
	6, tidyverse
	7, wesanderson
	8, propagate
	9, bde
	10, HDInterval


Sample preparation
Image file for analysis should be TIFF format with 1024 by 1024 pixel resolution, and have 3 channels with DAPI, pre-synaptic marker, post-synaptic marker respectively. This analysis requires 5 sequential z-stack image.
A sample IHC Image is in "Scripts_SampleData" folder.


Procedure
[Image process of Ideal Synapse]
1, Make "Ideal Synapse" folder, and make folder inside with sample ID name. Then place sample IHC image in the sample ID name folder.
1, Open ImageJ Fiji
2, run "ImageProcessing_IdealSynapse.ijm"
3, Select the directory of image file to analyze. 
4, During macro, dialog box appears and requires you to select OPL area with rectangle selection. Select OPL area, and press "OK". Macro will clear outside the ROI.
5, Once the image processing is done, csv files with graphical information of ROIs will be created. 
6, Open R Studio
7, Run "2_data.prep.fiji.R". This script gathers csv file information and export "data2D.RDa" that contains necessary data for R analysis.
8, Run "3_find_pair". This script identifiess pre- and post-synaptic markers pairs within 1.2 um, and export "data_pair.rds".
9, Run "4_estimate_bg.R" This script gathers background intensity information and exports "data_pair_bg.rds" and "data_bg.rds"
11, If you have several sample images to use as training data, iterate 1-10.

[Image process of Ideal Noise]
Same as Ideal Synapse, except that you need to create "Ideal Noise" folder at the beginning, and ImageJ macro you need to run is "ImageProcessing_IdealNoise.ijm"

[Generating PDFs]
This step generates probability density functions (PDFs) which contain characteristics of graphical parameters extracted from training data.
1, Place "Ideal Synapse" and "Ideal Noise" folder under same directory.
2, Run "MakeIdealSynapseDistribution.R" in the directory.
3, If you succeeded in running the script ,folder named "pdf" will appear in the directory.
4, Copy the "pdf" folder into the directory where you have sample IHC image to process with QUANTOS. Make sure you do this before you use QUANTOS.
5, Run QUANTOS as instructed in "Scripts_Sampledata".

