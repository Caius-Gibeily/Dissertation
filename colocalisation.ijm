// @String(label="Macro Mode",choices={"Current Image","Current Folder"}) modeChoice

// @String(value="JACOP Parameters", visibility="MESSAGE") JACOPparamtxt
// @Integer(label="Channel A",value=1) chA
// @Integer(label="Channel A",value=2) chB
// @String(label="What Threshold to use?",choices={ "Default","Huang ", "Intermodes ", "IsoData ", "IJ_IsoData ", "Li ", "MaxEntropy ", "Mean ", "MinError ", "Minimum ", "Moments ", "Otsu ", "Percentile ", "RenyiEntropy ", "Shanbhag ", "Triangle ", "Yen ","Use Manual Threshold Below"}) thresholdChoiceChA
// @String(label="What Threshold to use?",choices={ "Default","Huang ", "Intermodes ", "IsoData ", "IJ_IsoData ", "Li ", "MaxEntropy ", "Mean ", "MinError ", "Minimum ", "Moments ", "Otsu ", "Percentile ", "RenyiEntropy ", "Shanbhag ", "Triangle ", "Yen ","Use Manual Threshold Below"}) thresholdChoiceChB
// @Integer(label="Use Manual Threshold, Channel A",value=0) manualThChA
// @Integer(label="Use Manual Threshold, Channel B",value=0) manualThChB

// @Boolean(label="Crop ROIs", value=true) crop_rois_Bool
// @Boolean(label="Consider Z slices Separately", value=true) consider_z_slices_separately_Bool
// @Boolean(label="Set Auto Thresholds Stack Histogram", value=true) stackHisto_Bool
// @Boolean(label="Report as Vertical Montage", value=true) verticalReport_Bool
// @Boolean(label="All Coloc parameters", value=true) getMeasure_Bool
// @Boolean(label="Fluorogram", value=true) get_fluorogram_Bool


/*
 * romain dot guiet (at) epfl dot ch
 * olivier dot burri (at) epfl dot ch
 */
 
if (modeChoice == "Current Image"){
	testMode = true ;
} else {
	testMode = false ;
}


// Install the BIOP Library (from PTBIOP update site)
call("BIOP_LibInstaller.installLibrary", "BIOP"+File.separator+"BIOPLib.ijm");
//run("Close All");
//run("Clear Results");
//roiManager("Reset");
setForegroundColor(255, 255, 255);
setBackgroundColor(0, 0, 0);

run("Set Measurements...", "mean display redirect=None decimal=3");

if (testMode){
	if (nImages>0){
		processImage();	
	} else{
		showMessage("Please, open an image first.");
		return;
	}
	
} else { // process a folder
	setBatchMode(!testMode);
	
	//toolName();
	
	//get directory 
	imageFolder = getImageFolder();
	saveFolder = getSaveFolder();
	imageNumber = getNumberImages();

	for (imageIndex = 0 ; imageIndex < imageNumber ; imageIndex++){
		roiManager("Reset");
		openImage(imageIndex);					// open the image (an its assiciated roiset.zip)
		processImage();							// process the image
		//saveRois("Open");//saveRois("Save")		// save the ROIset 
		//saveCurrentImage();						// save the current image (within "\Processed")
		run("Close All");						// close all
	}
	
	// save the results
	if( isOpen("Results") ){
		selectWindow("Results");
		saveAs("Results", imageFolder+"results.txt");
	}
	

	setBatchMode(false);
}
showMessage("Jobs DONE!");


// required functions 

function toolName() {
	return "BIOP Basics Parameters";
}

function processImage(){
	imageName = getTitle();
	getDimensions(width, height, channels, slices, frames);
	getVoxelSize(widthPixel, heightPixel, depthPixel, unitPixel);
	rowIndexBefore = nResults;
	roiCount = roiManager("Count");

	
	reportOption = "";
	
	if (crop_rois_Bool){
		reportOption = reportOption+"crop_rois ";
	} 
	if (consider_z_slices_separately_Bool) {
		reportOption = reportOption+"consider_z_slices_separately ";
	}
	if (stackHisto_Bool) {
		reportOption = reportOption+"set_auto_thresholds_on_stack_histogram ";
	}
	if (getMeasure_Bool){
			reportOption = reportOption+"get_pearsons get_manders get_overlap get_li_ica ";
	}
	if (get_fluorogram_Bool){
			reportOption = reportOption+"get_fluorogram ";
	}
	if (verticalReport_Bool) {
		reportOption = reportOption+"report_as_vertical_montage ";
	}

	print(reportOption);
	
	imageNumberBefore = nImages;
					
	run("BIOP JACoP", 	"channel_a="+chA+" "+
						"channel_b="+chB+" "+
						"threshold_for_channel_a=["+thresholdChoiceChA+"] "+
						"threshold_for_channel_b=["+thresholdChoiceChB+"] "+
						"manual_threshold_a=["+manualThChA+"] "+
						"manual_threshold_b=["+manualThChB+"] "+
						reportOption);
	
	imageNumberAfter = nImages;
	
	
	imageToSaveNumber = imageNumberAfter - imageNumberBefore;
	
	for (selectedImageIndex = imageNumberBefore+1; selectedImageIndex <= imageNumberAfter ; selectedImageIndex++ ){
		//print(selectedImageIndex);
		//selectImage(selectedImageIndex);
		//saveCurrentImage();sss
	}
	
	//selectWindow("Results");
	updateResults;
}

