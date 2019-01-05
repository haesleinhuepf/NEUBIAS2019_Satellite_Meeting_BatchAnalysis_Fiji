filename = "C:/structure/data/Tischi_neubias-2019/neubias-2018/180222 sp5a.lif";

// threshold algorithm candidates are "Otsu", "Moments", "Yen", "Intermodes", "MaxEntropy"
vesicleThresholdingAlgorithm = "Otsu";

// zero-based channel indices to compute overlap
channelA = 1;
channelB = 3;

invertSelection = true;

function analyseImage(filename) {
	numberOfImagesInFile = 2;	
	for (i = 1 ; i < 2; i ++) {
		run("Close All");
		if ( endsWith(filename, ".tif") ) {
			numberOfImagesInFile = 1; // stop for loop after processing this image
			open(filename);
		} else {
			run("Bio-Formats", "open=[" + filename + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_" + i);
		}
		rename("input");
		
		pickSlice();

		segmentAllChannels();
		
		jaccardIndexWholeImage = computeOverlapBetweenChannels(channelA, channelB);

		print(filename + " Jaccard Index: " + jaccardIndexWholeImage);
	}
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

	imageCalculator("AND create", channel1ImageTitle, channel2ImageTitle);
	rename("Union");
	run("Create Selection");
	if (invertSelection) {
		run("Make Inverse");
	}
	run("Measure");
	areaUnion = getResult("Area", nResults() - 1);
	//close();

	imageCalculator("OR create", channel1ImageTitle, channel2ImageTitle);
	rename("Intersection");
	run("Create Selection");
	if (invertSelection) {
		run("Make Inverse");
	}
	run("Measure");
	areaIntersection = getResult("Area", nResults() - 1);
	//close();

	// close single channel images
	//close();
	//close();

	selectWindow(inputImageTitle);
	
	// return Jaccard index
	return areaUnion / areaIntersection;
}

run("Close All");

analyseImage(filename);



