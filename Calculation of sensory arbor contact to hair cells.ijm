folder_path = getDirectory("Choose folder to analyse");
subfolder_path = folder_path+"analysis/";  
filenames_list = getFileList(folder_path);


run("Bio-Formats Macro Extensions");
for(i=0;i<filenames_list.length;i++){
	if (matches(filenames_list[i], ".*\.nd") == 1){
		print(i);
		print(filenames_list[i]);
		namesplits_string = split(replace(filenames_list[i],".nd",""),"_");

		namesplitssplits_string = split(namesplits_string[2],"f");
		imagenumber_string = namesplitssplits_string[1];
		print(imagenumber_string);
		Ext.openImagePlus(folder_path+filenames_list[i]);
		origtitle_string = replace(getTitle(), ".TIF", "");
		print(origtitle_string);

		newtitle_string = replace(origtitle_string, "dpf11_", "dpf"+imagenumber_string+"_");
		print(newtitle_string);
		selectWindow(origtitle_string+".TIF");

		run("Duplicate...", "duplicate channels=1");
		close(origtitle_string+".TIF");
		selectWindow(origtitle_string+"-1.TIF");
		run("Z Project...", "projection=[Max Intensity]");
		close(origtitle_string+"-1.TIF");
		selectWindow("MAX_"+origtitle_string+"-1.TIF");
		saveAs("Tiff", subfolder_path+"MAX_"+newtitle_string+".tif");
		selectWindow("MAX_"+newtitle_string+".tif");
		saveAs("PNG", subfolder_path+"MAX_"+newtitle_string+".png");
		close("MAX_"+newtitle_string+".png");		
	}    
}
waitForUser("Make selection in Photoshop...");


subfilenames_list = getFileList(subfolder_path);
counter = 0;
counter2 = 0;
flag = false;
for(i=0;i<subfilenames_list.length;i++){
	if (matches(subfilenames_list[i], ".*\.png") == 1){
		title_string = replace(subfilenames_list[i],".png","");
		print(title_string);
		open(subfolder_path+title_string+".png");
		selectWindow(title_string+".png");
		run("8-bit");
		setAutoThreshold("Default dark");
		setOption("BlackBackground", true);
		run("Convert to Mask");
		run("Create Selection");
		roiManager("Add");
		roiManager("Deselect");
		roiManager("Delete");
		roiManager("Add");
		close(title_string+".png");
		open(subfolder_path+title_string+".tif");
		selectWindow(title_string+".tif");
		roiManager("Select", 0);
		run("Enlarge...", "enlarge=0.5");
		run("Make Inverse");
		roiManager("Add");
		counter2 = counter;
		for(j=0;j<filenames_list.length;j++){
			print(counter);
			print(counter2);
			if(matches(filenames_list[j], ".*\.zip")==1){
				if(flag==false){
					if(counter2==0){
						print(j);
						roiManager("Open", folder_path+filenames_list[j]);
						run("ROI Manager...");
						run("Select All");
						roiManager("Save", subfolder_path+filenames_list[j]);
						flag = true;

						roiManager("Deselect");
						roiManager("Select", newArray(0,2,3));
						roiManager("Delete");
						
						ROIcount_int = roiManager("count");
						ROIcount_int = roiManager("count");
						for(k=0;k<ROIcount_int;k++){
							roiManager("Deselect");
							roiManager("Select", k);
							label = Roi.getName;
							print(label);
							if(matches(label, "Path .1.-.*")||matches(label, "Path .2.-.*")==1){
								print("hit");
								roiManager("Delete");
								ROIcount_int--;
								k--;
							}
						}
						
						
						ROIcount_int=roiManager("count");
						print(ROIcount_int);
						selection_array = newArray(ROIcount_int);
						for (k=0;k<ROIcount_int;k++){
							selection_array[k] = k;
						}
						selection_array[0]=1;
						roiManager("Select", selection_array);
						roiManager("Combine");
						roiManager("Add");
						roiManager("Delete");
						roiManager("Select", newArray(0,1));
						roiManager("AND");
						roiManager("Add");
						roiManager("Deselect")
						roiManager("Select", 0);
						roiManager("Delete");

						run("ROI Manager...");
						run("Select All");
						roiManager("Save", getInfo("image.directory")+getInfo("image.filename")+".zip");

//						run("Set Measurements...", "stack redirect=None decimal=0");
						roiManager("Select", 0);
//						run("Interpolate", "interval=1");
						Roi.getContainedPoints(x, y);
					    for (k=0; k<x.length; k++) {
							setResult("Label", k, 0);
							setResult("X", k, x[k]);
							setResult("Y", k, y[k]);
							setResult("Name", k, Roi.getName);
						}
						updateResults();
						saveAs("Results", getInfo("image.directory")+getInfo("image.filename")+"_tree.csv");
						run("Clear Results");
						roiManager("Select", 1);
//						run("Interpolate", "interval=1");
						Roi.getContainedPoints(x, y);
					    for (k=0; k<x.length; k++) {
							setResult("Label", k, 1);
							setResult("X", k, x[k]);
							setResult("Y", k, y[k]);
							setResult("Name", k, Roi.getName);
						}
						updateResults();   
						saveAs("Results", getInfo("image.directory")+getInfo("image.filename")+"_branch.csv");
						run("Clear Results");
					}
					else{
						counter2 = counter2-1;
					}
				}
			}
		}
		counter = counter+1;
		flag = false;
	}    
}