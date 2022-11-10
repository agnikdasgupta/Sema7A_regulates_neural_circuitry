

roiManager("select",0);
Roi.getContainedPoints(x_top,y_top);
makeLine(x_top[0],y_top[0],x_top[0]+0.001,y_top[0]);
roiManager("Add");

roiManager("open", getInfo("image.directory")+replace(getInfo("image.filename"),"_w1iSIM488-525","_w2iSIM561-630")+".zip");

roiManager("deselect");
roiManager("select",newArray(0,2,3));
roiManager("delete");

ROIcount_int=roiManager("count");
flag = 0;
for(i=0;i<ROIcount_int;i++){
	roiManager("Deselect");
	roiManager("Select", i);
	name = Roi.getName;
//	if(matches(name, "Path .3.-XY-0001")==1){
//		roiManager("Select", i-1);
//		break;
//	}
	if(matches(name, "Path .2.-XY-0001")==1){
		flag = 1;
	}
	if(flag==1&&matches(name, "Path .2.-XY-0001")==0){
		flag = 0;
		roiManager("Select", i-1);
		break;
	}
}
position_bot = Roi.getPosition(channel, slice_bot, frame);
Roi.getContainedPoints(x_bot,y_bot);
makeLine(x_bot[0],y_bot[0],x_bot[0]+0.001,y_bot[0]);
roiManager("Add");

ROIcount_int = roiManager("count");
for(k=0;k<ROIcount_int;k++){
	roiManager("Deselect");
	roiManager("Select", k);
	label = Roi.getName;
	if(matches(label, "Path .1.-.*")||matches(label, "Path .2.-.*")==1){
		roiManager("Delete");
		ROIcount_int--;
		k--;
	}
}



run("Clear Results");
ROIcount_int=roiManager("count");
roiManager("Deselect");
n=0;
for(i=0;i<ROIcount_int;i++){
	roiManager("Select", i);
	name=Roi.getName;
//	Roi.getPosition(channel, z, frame);
	z=getSliceNumber();
	run("Line to Area");
	Roi.getContainedPoints(x,y);
	for(j=0;j<x.length;j++){	
		setResult("Name", j+n, name);
		setResult("X", j+n, x[j]);
		setResult("Y", j+n, y[j]);
		setResult("Z", j+n, z);
	}
	n+=x.length;
	updateResults();
}


roiManager("Save", getInfo("image.directory")+replace(getInfo("image.filename"),"_w1iSIM488-525","_w2iSIM561-630")+"-top_and_bot.zip");
updateResults();
saveAs("Results", getInfo("image.directory")+replace(getInfo("image.filename"),"_w1iSIM488-525","_w2iSIM561-630")+"-top_and_bot.csv");

//pos = eval("js","IJ.getImage().getRoi().getPosition()");
//print("position="+pos);

//getDimensions(width, height, channels, slices, frames);

roiManager("deselect");
roiManager("delete");
close("*");