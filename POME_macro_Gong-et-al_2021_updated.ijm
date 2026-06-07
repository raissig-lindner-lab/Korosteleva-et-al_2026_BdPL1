// Part I: Initialization
Dialog.create("POME 1.0.0");
Dialog.addString("Membrane Width: ", 1);
Dialog.addString("Minimum Cell Size: ", 10); 
Dialog.addString("Cell Radius: ", 100);
Dialog.addString("Minimum Circularity: ", 0.4);
Dialog.addString("Cell Outline Channel: ", 2);
Dialog.addString("Polarity Protein Channel: ", 1);


Dialog.addMessage("Input Details:");

Dialog.addMessage("Membrane Width: Width of the membrane\ndetermines the number of pixels to be\nsegmented as the membrane and measured\nin POME. Unit: physical units (e.g. microns,\nin images with scale) or pixels (in images\nwithout scale).");

Dialog.addMessage("Minimal Cell Size: This parameter helps\nto exclude interference (e.g. small particles,\nnucleus) to be segmented as a cell.\nAnything below this size will be excluded.\nUnit: same as the width of the membrane.");

Dialog.addMessage("Cell Radius: This parameter determines\nthe number of pixels for the measurement\nline originated from the cell centroid. Unit:\npixels.");

Dialog.addMessage("Minimum Circularity: In some cases,\nminimal circularity of the cell needs\nto be lowered (e.g. to 0.3) here to segment\ncells with extreme shape or aspect ratio.");

Dialog.addMessage("Cell Outline Channel: The channel\ncontaining the cell outline marker.");

Dialog.addMessage("Polarity Protein Channel: The channel\ncontaining the marker of interest.");

//Dialog.addHelp("file:///C:/Raghav/Stuff/Important%20Things/College/Research/CheungLab/Thing/POME_paper_eLIFE_submission/Appendix_1_Brief_user_guide_for_POME_analysis.pdf")

Dialog.show();

MembraneSize=Dialog.getString();// Set the membrane width in scale unit (e.g. micron) in images with scale or number of pixel in images without scale
MinCellSize=Dialog.getString();// Set the minimal cell size in scale unit (e.g. square micron) in images with scale or number of pixel in images without scale
CellRadius=Dialog.getString();// Set the length of the rotating measurement line in number of pixel
MinCirc = Dialog.getString();
mCh = Dialog.getString();
gfpCh = Dialog.getString();



// Part II: Segmentation
waitForUser("Please open an image and select the region of interest, then click 'OK'.");

numb = nSlices();
bool = numb > 2;

while (bool) {
	waitForUser("Please convert your z-stack to a single image for analysis. Open a new image and select a different region of interest, then click 'OK' to retry or click 'Esc' on keyboard to abort.");
	numb = nSlices();
	bool = numb > 2;
}



run("Set Measurements...", "area mean min centroid center redirect=None decimal=0");
run("Duplicate...", "duplicate");
run("Make Composite"); 
Stack.setChannel(mCh);
rename("Sample");
selectWindow("Sample");
run("Duplicate...", "title=C" + mCh + "-mask channels=0");
run("Despeckle");
run("Invert");
setAutoThreshold("Huang dark");
run("Convert to Mask");
run("Analyze Particles...", "size=MinCellSize-Infinity circularity=MinCirc-1.00 exclude add");// Adjust the lower boundary of the circularity here for cells with extrene shape and aspect ratio
rCount = roiManager("Count");




while (rCount < 1) {
	selectWindow("Sample");
	close();
	selectWindow("C" + mCh + "-mask");
	close();
	waitForUser("Sorry, ROI doesn't work. Please open an image and select a different region of interest, then click 'OK' to retry or click 'Esc' on keyboard to abort.");
	run("Set Measurements...", "area mean min centroid center redirect=None decimal=0");
	run("Duplicate...", "duplicate");
	rename("Sample");
	selectWindow("Sample");
	run("Duplicate...", "title=C" + mCh + "-mask channels=0");
	run("Despeckle");
	run("Invert");
	setAutoThreshold("Huang dark");
	run("Convert to Mask");
	run("Analyze Particles...", "size=MinCellSize-Infinity circularity=MinCirc-1.00 exclude add");// Adjust the lower boundary of the circularity here for cells with extrene shape and aspect ratio
	rCount = roiManager("Count");
}

roiManager("Select", 0); //ERRORS HERE
run("Enlarge...", "enlarge=MembraneSize");
roiManager("Add");
roiManager("Select", newArray(0,1));
roiManager("XOR");
roiManager("Add");


selectWindow("Sample");
run("Split Channels");

selectWindow("C" + gfpCh + "-Sample");
run("Green");
run("Duplicate...", "title=C" + gfpCh + "-pixels");
run("32-bit");
run("Duplicate...", "title=C" + gfpCh + "-pixelsmask");

selectWindow("C" + mCh + "-Sample");
run("Magenta");
roiManager("Deselect");
roiManager("Select", 1);
roiManager("Measure");
AreaX = getResult("X");
AreaY = getResult("Y");
roiManager("Select", 2);
roiManager("Measure");
MembraneXM = getResult("XM");
MembraneYM = getResult("YM");
run("Color Picker...");
setBackgroundColor(0, 0, 0);
run("Clear Outside");
selectWindow("C" + gfpCh + "-Sample");
roiManager("Select", 2);
roiManager("Measure");
PolarizedXM = getResult("XM");
PolarizedYM = getResult("YM");
run("Color Picker...");
setBackgroundColor(0, 0, 0);
run("Clear Outside");
IJ.renameResults("CenterOfMass_Results");
saCyt = getResult("Area",0);
memReg = getResult("Area",1);
centX = getResult("X",1);
centY = getResult("Y",1);
comX = getResult("XM",2);
comY = getResult("YM",2);
print("Cytosol Surface Area: " + toString(saCyt));
print("Membrane Region Area: " + toString(memReg));
print("Cell Centroid Location: (" + toString(centX) + ", " + toString(centY) + ")");
print("Center of Mass of Polarity Protein Channel: (" + toString(comX) + ", " + toString(comY) + ")");

selectWindow("C" + gfpCh + "-pixels");
roiManager("Select", 2);
run("Clear Outside");
selectWindow("C" + gfpCh + "-pixelsmask");
roiManager("Select", 2);
setForegroundColor(0, 255, 0);
run("Fill", "slice");
run("Clear Outside");
run("Convert to Mask");
run("Divide...", "value=255");
imageCalculator("Divide create 32-bit", "C" + gfpCh + "-pixels","C" + gfpCh + "-pixelsmask");
selectWindow("C" + gfpCh + "-pixels");
run("Close");
selectWindow("C" + gfpCh + "-pixelsmask");
run("Close");

selectWindow("C" + mCh + "-mask");
run("Cyan");
run("Select All");
setBackgroundColor(0, 0, 0);
run("Clear", "slice");
roiManager("Select", 2);
setForegroundColor(0, 0, 255);
run("Draw", "slice");
roiManager("Select", newArray(0,1,2));
run("Select All");
roiManager("Delete");
selectWindow("CP");
run("Close");

// Part III: Measure fluorescence intensity in radius
Angle = 0;

//Use constant angle (GW)
/*Angle = atan((PolarizedYM-MembraneYM)/(PolarizedXM-MembraneXM));
(if (PolarizedXM-MembraneXM<0) {
	Angle = Angle+PI;
	}
else if (PolarizedXM-MembraneXM>0 && PolarizedYM-MembraneYM<0) {
		Angle = 2*PI+Angle;
	}
*/

	
getPixelSize(unit, pixelWidth, pixelHeight);
X = AreaX/pixelWidth;
Y = AreaY/pixelHeight;
selectWindow("Result of C" + gfpCh + "-pixels");

row = 0;
	for (n=Angle-5; n<=Angle+5; n += 0.1) {
	selectWindow("Result of C"+ gfpCh + "-pixels");
	makeLine(X, Y, X+CellRadius*cos(n), Y+CellRadius*sin(n));
	setResult("Angle", row, n*180/PI);
 	run("Plot Profile");
 	Plot.getValues(x, y);
	 	for (i=1; i<x.length; i++) {
 			if (isNaN(y[i-1]) && y[i]==0) {
 				y[i]= NaN;
 			}
  		setResult(x[i-1], row, y[i-1]);
	 	}
	Truey = newArray();
		sum = 0;
 		for (j=0; j<y.length; j++) {
  			if (isNaN(y[j])==0) {
			Truey = Array.concat(Truey, y[j]);
			sum = sum + y[j];
  			}
 		}
 	selectWindow("Plot of Result of C" + gfpCh + "-pixels");
 	run("Close");
 	Array.getStatistics(Truey, Min, Max, Mean, Std);
 	setResult("Mean", row, Mean);
 	setResult("Max", row, Max);
 	setResult("Std", row, Std);
 	setResult("Sum", row, sum);
    row++;
	}

IJ.renameResults("Polarity_Results");

// Part IV: Clean up and output
print("Angle of Polarization = ", Angle*180/PI, "\n");	
selectWindow("C"+mCh+"-mask");
close();
run("Merge Channels...", "c6=C" + mCh + "-Sample c" + gfpCh + "=C" + gfpCh + "-Sample");
rename("POME 1.0.0 Segmentation Visualization");
selectWindow("POME 1.0.0 Segmentation Visualization");
makeLine(X, Y, X+100*cos(Angle), Y+100*sin(Angle));

setForegroundColor(255, 255, 255);
run("Draw", "slice");
selectWindow("Result of C" + gfpCh + "-pixels");
close();
selectWindow("CenterOfMass_Results");
run("Close");
selectWindow("ROI Manager");
run("Close");
IJ.renameResults("Polarity_Results","POME 1.0.0 Polarity Results");