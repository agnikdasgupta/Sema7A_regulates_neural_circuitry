run("ROI Manager...");
run("Select All");
roiManager("Deselect");
roiManager("Delete");
roiManager("Open", getInfo("image.directory")+getInfo("image.filename")+".zip")
run("ROI Manager...");
run("Select All");
roiManager("Save", getInfo("image.directory")+getInfo("image.filename")+".zip")
numROIs=roiManager("count");
nr=0;
run("Set Measurements...", "stack redirect=None decimal=3");

for (i=0; i<numROIs; i++) {
	roiManager("Select", i);
	b = selectionType();
	if (b==10) {
		roiManager("measure");
		setResult("Name", i, Roi.getName);
		saveAs("Results", getInfo("image.directory")+getInfo("image.filename")+"-"+Roi.getName+".csv");
		run("Clear Results");
	}
	if (b==7) {
		getSelectionCoordinates(x, y);
		z = getSliceNumber();
		for (j=0; j<x.length; j++) {
			setResult("X", j+nr, x[j]);
			setResult("Y", j+nr, y[j]);
			setResult("Z", j+nr, z);
			setResult("Name", j+nr, Roi.getName);
		}
		nr+=x.length;
		updateResults();
	}
}
saveAs("Results", getInfo("image.directory")+getInfo("image.filename")+".csv");
run("ROI Manager...");