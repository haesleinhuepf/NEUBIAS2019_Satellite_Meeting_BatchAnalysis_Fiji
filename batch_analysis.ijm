// Vesice colocalisation batch processing macro
// ============================================
// 
// Expected input data: 
//  * Folder of lif files
//  * Four channels each
//    * 1: DAPI / nuclei channel
//    * 3, 4: vesicles
//
// Computes the Jaccard index / overlap between channel 3 and 4
// 
// This is a academic code example prepared for the NEUBIAS Satellite
// meeting BioImage Analyst Best Practices seminar 2019.   
//
// Author: Robert Haase, Myers lab, MPICBG, rhaase@mpi-cbg.de
// January 2019
// ----------------------------------------
inputFolder = "C:/structure/data/Tischi_neubias-2019/neubias-2018/";
resultCSVFilename = "C:/structure/data/Tischi_neubias-2019/result.tsv";

// threshold algorithm candidates are "Otsu", "Moments", "Yen", "Intermodes", "MaxEntropy"
vesicleThresholdingAlgorithm = "Otsu";

// zero-based channel indices to compute overlap
channelA = 1;
channelB = 3;
nucleiChannel = 0;

binaryOpeningDistance = 3;

// main: run batch processing of a folder
function main() {
	File.append("Image/Cell\tJaccard Index", resultCSVFilename);
	
	processFolder(inputFolder);
}

// go through a folder and open all lif images
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);

	for (i = 0; i < list.length; i++) {
		if(endsWith(list[i], ".lif")) {
			analyseFile(input + list[i]);
			break;
		}
	}
}

// open one particular image
function analyseFile(filename) {
	
	numberOfImagesInFile = 1000;
	firstTitle = "";
	for (i = 1 ; i < numberOfImagesInFile; i ++) {
		run("Close All");
		run("Bio-Formats", "open=[" + filename + "] " + 
		"autoscale color_mode=Default rois_import=[ROI manager] " +
		"view=Hyperstack stack_order=XYCZT series_" + i);
		
		if (getTitle() == firstTitle) {
			// no image open
			break;
		}
		getDimensions(width, height, channels, slices, frames);
		if (slices < 2 || channels != 4) {
			continue;
		}
		if (firstTitle == "") {
			firstTitle = getTitle();
		}
	
		rename("input"); 
		analyseImage(filename, filename + "_" + i);
	}
}

// analyse one particulare image by segmenting cells, vesicles and measuring
// overlap
function analyseImage(inputFilename, outputFilename) {
	pickSlice();
	if (!File.exists(outputFilename + ".slice.tif")) {
		saveAs("Tiff", outputFilename + ".slice.tif");
	}
	segmentAllChannels();
	
	// analysis of overlap in the whole image
	jaccardIndexWholeImage = computeOverlapBetweenChannels(channelA, channelB);
	
	print(filename + " Jaccard Index: " + jaccardIndexWholeImage);
	
	count = 0;
	sumJaccardIndex = 0;
	while(true) {
		count = count + 1;
		cellFound = segmentCell(outputFilename, count);
		if (!cellFound) {
			break;
		}

		// analysis of overlap in individual cell
		jaccardIndexIndividualCell = computeOverlapBetweenChannels(channelA, channelB);

		// summarize over whole image
		sumJaccardIndex	= sumJaccardIndex + jaccardIndexIndividualCell;
		
		print(outputFilename + " cell " + count + " Jaccard Index: " + jaccardIndexIndividualCell);
		
		File.append(outputFilename + " cell " + count + "\t" + jaccardIndexIndividualCell, resultCSVFilename);
		
		close();
	}

	// determine average of cells in one image and write to result log file
	File.append(inputFilename + " cell average\t" + (sumJaccardIndex / count), resultCSVFilename);
	File.append(inputFilename + " whole image\t" + jaccardIndexWholeImage, resultCSVFilename);
	
}

// cell segmentation using nuclei segmentation + voronoi maps
// vesicle segmentation using Otsu method
function segmentCell(outputFilename, cellNumber) {
	run("Duplicate...", "duplicate ");
	run("Duplicate...", "duplicate channels=" + (nucleiChannel + 1));
	run("Fill Holes");

	
	for (o = 0; o <binaryOpeningDistance; o++) {
		run("Dilate", "slice");
	}
	for (o = 0; o <binaryOpeningDistance; o++) {
		run("Erode", "slice");
	}

	// segment cells using a voronoi diagram
	setOption("BlackBackground", true);
	run("Ultimate Points");
	run("Make Binary");
	run("Voronoi");
	setAutoThreshold("Otsu dark");
	setThreshold(0, 1);
	run("Convert to Mask");

	// connected components labeling
	run("Analyze Particles...", "  show=[Count Masks]");
	if (!File.exists(outputFilename + ".cells.tif")) {
		saveAs("Tiff", outputFilename + ".cells.tif");
	}
	
	// select only one cell
	setThreshold(cellNumber, cellNumber);
	run("Convert to Mask");

	// check if something was selected
	run("Set Measurements...", "mean redirect=None decimal=3");
	run("Measure");
	mean = getResult("Mean", nResults() - 1);
	if (mean ==0) { 
		return false;
	}

	// if yes, set all other pixels to zero
	run("Create Selection");
	getDimensions(width, height, channels, slices, frames);
	
	run("Make Inverse");
	close();
	close();
	run("Restore Selection");
	run("Multiply...", "value=0 stack");
	run("Select None");
	return true;
}

// Use channel 1 (DAPI) to find the slice where the most information is in. Concentrate on this plane.
function pickSlice() {
	// an alternative could be maximum projection; for colocalization use with care
	//run("Z Project...", "projection=[Max Intensity]");
	
	getDimensions(width, height, channels, slices, frames);
	
	run("Set Measurements...", "mean redirect=None decimal=3");
	maximumIntensity = -1;
	argMax = 0;
	for (z = 0; z < slices; z++) {
		Stack.setSlice(z + 1); // slices are 1-indiced
		run("Measure");
		intensity = getResult("Mean", nResults() - 1);
		if (intensity > maximumIntensity) {
			maximumIntensity = intensity;
			argMax = z;
		}
	}
	run("Duplicate...", "duplicate slices=" + (argMax + 1));
}

function segmentAllChannels() {
	setOption("BlackBackground", true);
	run("Convert to Mask", "method=" + vesicleThresholdingAlgorithm + " background=Dark calculate black");
}

function computeOverlapBetweenChannels(channel1, channel2) {
	inputImageTitle = getTitle();
	
	run("Duplicate...", "duplicate channels=" + (channel1 + 1));
	rename("channelA");
	channel1ImageTitle = getTitle();

	selectWindow(inputImageTitle);
	run("Duplicate...", "duplicate channels=" + (channel2 + 1));
	rename("channelB");
	channel2ImageTitle = getTitle();

	run("Set Measurements...", "area redirect=None decimal=3");

	// determine union area
	imageCalculator("OR create", channel1ImageTitle, channel2ImageTitle);
	rename("Union");
	setThreshold(127, 255);
	run("Create Selection");
	run("Measure");
	areaUnion = getResult("Area", nResults() - 1);
	
	// determine intersection area
	imageCalculator("AND create", channel1ImageTitle, channel2ImageTitle);
	rename("Intersection");
	setThreshold(127, 255);
	run("Create Selection");
	run("Measure");
	areaIntersection = getResult("Area", nResults() - 1);

	// close single channel images
	selectWindow("channelA");
	close();
	selectWindow("channelB");
	close();
	selectWindow("Union");
	close();
	selectWindow("Intersection");
	close();
	
	selectWindow(inputImageTitle);
	
	// return Jaccard index
	return areaIntersection / areaUnion;
}

run("Close All");

main();

