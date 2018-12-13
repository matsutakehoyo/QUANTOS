setBatchMode(true);
roiManager("Reset");
run("Clear Results");

//Getting directory path
path = getDirectory("Select a Directory");

//setting working directory of Fiji
open(path + File.separator + "image.tif");
name = getTitle();
print(name);

//Making OPL ROI
selectWindow("image.tif");
Stack.setDisplayMode("composite");
setTool("polygon");
waitForUser("Select OPL area then press OK");
roiManager("Add");
roiManager("Save",  path + File.separator + "OPL.zip");
roiManager("Reset");
roiManager("Show All");
roiManager("Show None");

//delete all ROIs | avoid errors in excution with preexisting ROIs
roiManager("Reset");
//this is the original imag
selectWindow("image.tif");
ori = getImageID();
run("Duplicate...", "duplicate");
dup = getImageID();
//Close original image
selectImage(ori);
run("Close");
//Split channels
selectImage(dup);
run("Split Channels");


//tag DAPI image
selectImage(1);
dapi=getImageID();
run("Blue");
//tag postynaptic marker image
selectImage(2);
post=getImageID();
run("Red");
saveAs("Tiff", path + File.separator +  "post.tif");
run("Close");
//tag presynaptic marker image
selectImage(2);
pre=getImageID();
run("Green");
saveAs("Tiff", path + File.separator +  "pre.tif");
run("Close");

//select DAPI image
selectImage(dapi);
run("Z Project...", "start=2 stop=6 projection=[Average Intensity]");
saveAs("Tiff", path + File.separator +  "DAPIZproject.tif");
dapiZ=getImageID();
selectImage(dapi);
run("Close");
run("Duplicate...", "duplicate");
//run("Subtact background...")
run("Subtract Background...", "rolling=18");
//run("Bandpass filter...)
run("Bandpass Filter...", "filter_large=32 filter_small=6 suppress=None tolerance=5 process");
saveAs("Tiff", path + File.separator + "DAPIZproject_bandpass.tif");
dapiZbandpass=getImageID();
//RATS
run("Robust Automatic Threshold Selection", "noise=3 lambda=1 min=205");
saveAs("Tiff", path + File.separator +  "DAPIZproject_threshold.tif");
run("Options...", "iterations=1 count=1 black do=Dilate");
//run("Invert LUT");
run("Adjustable Watershed", "tolerance=0.5");
saveAs("Tiff", path + File.separator +  "DAPIZproject_watershed.tif");
run("Close All");
//Analyze particle
setBatchMode(false);
open(path + File.separator + "DAPIZproject_watershed.tif");
run("Analyze Particles...", "clear add slice");
roiManager("Save", path + File.separator + "ROIset_dapi.zip");
run("Close");
//Measure
open(path + File.separator + "DAPIZproject.tif");
roiManager("Show None");
roiManager("Show All");
run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis area_fraction limit display scientific add nan redirect=None decimal=9");
roiManager("Measure");
saveAs("Results", path + File.separator +  "DAPI.csv");
run("Clear Results");
roiManager("Reset");
run("Close All");


//Select Post
setBatchMode(true);

open(path + File.separator + "post.tif");
post=getImageID();
//Make PostZ
run("Z Project...", "start=2 stop=6 projection=[Average Intensity]");
saveAs("Tiff", path + File.separator +  "PostZproject.tif");
postZ=getImageID();
selectImage(post);
run("Close");
//delete all ROIs | avoid errors in excution with preexisting ROIs
roiManager("Reset");
//Get rid of noise 
selectImage(postZ);
run("Duplicate...", "duplicate");
run("Bandpass Filter...", "filter_large=20 filter_small=1 suppress=None tolerance=5 process");
run("Maximum...", "radius=0.5");
saveAs("Tiff", path + File.separator +  "PostZproject_bandpass.tif");
postZbandpass=getImageID();

//create empty image
selectImage(postZbandpass);
run("Duplicate...", "duplicate");
run("Select All");
run("Clear", "slice");
postEmpty=getImageID();
wait(800);
//Threshold individual ROIs and Interporate into 1 image
selectImage(postZbandpass);
iniT = 8;
run("Find Maxima...", "noise=8 output=[Segmented Particles]");
run("Analyze Particles...", "size=0-infinity clear add");
roiManager("Save", path + File.separator + "PreMaxROI.zip");
nROI = roiManager("count");
print(nROI);

//Loop for Every ROI and Threshold Skeletonize
for(i=0; i<nROI; i++){
	print(i);
	selectImage(postZbandpass);
	run("Duplicate...", "duplicate");
	temp=getImageID();
	roiManager("select", i);
	run("Clear Outside");
	getStatistics(area, mean, min, max, std, histogram);
	//print(max);
	//print(min);
	tol_i=(max - min)/2.6;
	if(tol_i<iniT){
		tio_i=iniT;
	}
	print(tol_i);
	run("Select None");
	//run("Robust Automatic Threshold Selection", "noise=5 lambda=3 min=205");
	//run("Sharpen");
	run("Find Maxima...", "noise=&tol_i output=[Maxima Within Tolerance]");
	temp2=getImageID();
	run("Select All");
	getStatistics(area, mean);
	//print(mean);
	run("Select None");
	if(mean>254){
	run("Invert");
	}
	imageCalculator("Add add", postEmpty, temp2);
	selectImage(temp);
	run("Close");
	selectImage(temp2);
	run("Close");
	//Clear memory
	//call("java.lang.System.gc");
}
selectImage(postEmpty);
saveAs("Tiff", path + File.separator +  "PostZproject_MaximaThreshold.tif");
saveAs("Tiff", path + File.separator +  "postZproject_Threshold_trim.tif");
//Clear  OPL
roiManager("Reset");
roiManager("Open", path + File.separator + "OPL.zip");
roiManager("select", 0);
run("Clear");
roiManager("Reset");
roiManager("Show None");
roiManager("Show All");
saveAs("Tiff", path + File.separator +  "postZproject_OPL.tif");
run("Analyze Particles...", "size=0-2 clear add");
roiManager("Save",  path + File.separator + "ROIset_post.zip");
//Use Post ROI for postZ
open(path + File.separator + "postZproject.tif");
postZ=getImageID();
roiManager("Show None");
roiManager("Show All");
run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis area_fraction limit display scientific add nan redirect=None decimal=9");
roiManager("Measure");
saveAs("Results", path + File.separator +  "Post.csv");
roiManager("Reset");
run("Clear Results");
run("Close All");

//Clear memory
call("java.lang.System.gc");



roiManager("Reset");

setBatchMode(true);

//Open Pre
open(path + File.separator + "pre.tif");
pre=getImageID();
//Make PreZ
run("Z Project...", "start=2 stop=6 projection=[Average Intensity]");
saveAs("Tiff", path + File.separator +  "PreZproject.tif");
open(path + File.separator + "preZproject.tif");
preZ=getImageID();
run("Duplicate...", "duplicate");
preZdup=getImageID();
//delete all ROIs | avoid errors in excution with preexisting ROIs
roiManager("Reset");
run("Subtract Background...", "rolling=10");
//run(Bandpass Filter...)
run("Bandpass Filter...", "filter_large=20 filter_small=1 suppress=None tolerance=5 process");
run("Smooth");
saveAs("Tiff", path + File.separator +  "PreZproject_bandpass.tif");
preZbandpass=getImageID();
//Get rid of noise 
run("Find Maxima...", "noise=90 output=[Segmented Particles]");
wait(800);

//create empty image
selectImage(preZbandpass);
run("Duplicate...", "duplicate");
run("Select All");
run("Clear", "slice");
preEmpty=getImageID();
wait(800);
//Threshold individual ROIs and Interporate into 1 image
selectImage(preZbandpass);
run("Find Maxima...", "noise=10 output=[Segmented Particles]");
run("Analyze Particles...", "size=0-infinity clear add");
roiManager("Save", path + File.separator + "PreMaxROI.zip");
nROI = roiManager("count");
print(nROI);

//Loop for Every ROI and Threshold Skeletonize
for(i=0; i<nROI; i++){
	print(i);
	selectImage(preZbandpass);
	run("Duplicate...", "duplicate");
	temp=getImageID();
	roiManager("select", i);
	run("Clear Outside");
	getStatistics(area, mean, min, max, std, histogram);
	print(max);
	print(min);
	tol_i=(max - min)/1.5;
	print(tol_i);
	run("Select None");
	//run("Robust Automatic Threshold Selection", "noise=5 lambda=3 min=205");
	//run("Sharpen");
	run("Find Maxima...", "noise=&tol_i output=[Maxima Within Tolerance]");
	temp2=getImageID();
	run("Select All");
	getStatistics(area, mean);
	print(mean);
	run("Select None");
	if(mean>254){
	run("Invert");
	}
	imageCalculator("Add add", preEmpty, temp2);
	selectImage(temp);
	run("Close");
	selectImage(temp2);
	run("Close");
}
selectImage(preEmpty);
saveAs("Tiff", path + File.separator +  "PreZproject_skeletonize.tif");
saveAs("Tiff", path + File.separator +  "PreZproject_skeletonize_trim.tif");
//Clear  OPL
roiManager("Reset");
roiManager("Open", path + File.separator + "OPL.zip");
roiManager("select", 0);
run("Clear");
roiManager("Reset");
roiManager("Show None");
roiManager("Show All");
saveAs("Tiff", path + File.separator +  "PreZproject_skeletonize_OPL.tif");
run("Close All");
run("Close All");
setBatchMode(false);
//Analyze Particlew
open(path + File.separator + "PreZproject_skeletonize_trim.tif");
run("Analyze Particles...", "size=0.2-4.0 clear add");
roiManager("Save",  path + File.separator + "ROIset_pre.zip");
//Measure
open(path + File.separator + "PreZproject.tif");
preZ=getImageID();
roiManager("Show None");
roiManager("Show All");
run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis area_fraction limit display scientific add nan redirect=None decimal=9");
roiManager("Measure");
saveAs("Results", path + File.separator +  "Pre.csv");
roiManager("Reset");
run("Clear Results");
run("Close All");

//Measure DAPI again
open(path + File.separator + "DAPIZproject.tif");
roiManager("Open", path + File.separator + "ROIset_dapi.zip");
roiManager("Show None");
roiManager("Show All");
run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis area_fraction limit display scientific add nan redirect=None decimal=9");
roiManager("Measure");
saveAs("Results", path + File.separator +  "DAPI.csv");
run("Clear Results");
roiManager("Reset");
run("Close All");


//Getting Text image for Backgroundestimation
	//DAPI
	open(path + File.separator + "DAPIZproject.tif");
	saveAs("Text Image",  path + File.separator + "DAPIZproject.txt");
	
	
	//Pre
	open(path + File.separator + "PreZproject.tif");
	saveAs("Text Image",  path + File.separator + "PreZproject.txt");
	
	//Post
	open(path + File.separator + "PostZproject.tif");
	saveAs("Text Image",  path + File.separator + "PostZproject.txt");
	
	run("Close All");

//Clear memory
call("java.lang.System.gc");

//print("quitting...");
//run("Quit");
print("quitting...");
eval("script", "System.exit(0);"); 