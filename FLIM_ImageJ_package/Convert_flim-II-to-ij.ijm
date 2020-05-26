// Macro to convert FLIM files from FLIM-II for imageJ version 1
//===========================================================================
// Version history:
// Version 1: 28-2-2019 Dorus Gadella, created from Convert_flim-I-to_ij_v2.
//===========================================================================
// It needs an opened reference image stack. This reference is used for the entire directory (if that option is selected).
// Reads in metadata from FLIM-II FLIM stacks (.fli files) and converts them for input in ImageJ for lifetime display.
// Reference phase and modulation and frequency are imported from a reference stack
// Optional background subtraction per phase image. This can be a signed 16-bit correction or a rolling ball correction.
// Includes calculation of average image lifetime statistics based on a fixed background threshold.
//===========================================================================


Dialog.create("Convert FLIM-II ics stack to input for ImageJ Lifetimes15 macro");
Dialog.addChoice("What to convert ",newArray("Convert file","Convert entire directory"),"Convert entire directory");
Dialog.addChoice("Background subtraction in phase images ", newArray("no background subtraction","rolling ball with 500 pixels diameter","rolling ball with 250 pixels diameter","rolling ball with 125 pixels diameter", "sliding parabola with 500 pixels diameter","sliding parabola with 250 pixels diameter","sliding parabola with 125 pixels diameter"),"no background subtraction");
Dialog.addNumber("Gaussian blur sample data, 0=no, >0 is blur size", 1);
Dialog.addNumber("Reference lifetime (ns):", 4.05);
Dialog.addNumber("Gaussian blur reference data, 0=no, >0 is blur size", 4);
Dialog.addCheckbox("Calculate average image statistics:", true);
Dialog.addNumber("Background threshold for image stats", 1000);
Dialog.addCheckbox("Write logfile with metadata and statistics:", true);
Dialog.addString("Output extension of filename", "_ij.tif");
Dialog.show();
choice_open= Dialog.getChoice();
bg_choice= Dialog.getChoice();
blur_choice = Dialog.getNumber();
tauref= Dialog.getNumber();
blur_ref= Dialog.getNumber();
stat= Dialog.getCheckbox();
bg_stat= Dialog.getNumber();
log_choice= Dialog.getCheckbox();
extension = Dialog.getString();
//end input


hfreq=40000000;
phasecorr0=0;
modcorr0=1.0;
logje=0;
v=newArray(160);

timepoints=newArray(10000);
timepoints[0]=0;

suffix=".fli";
getDateAndTime(year,month,dw,dm,hr,mi,sec,msec);
month=month+1;

//==============
setBatchMode(true);
waitForUser("Please select a reference phase stack")
run("Bio-Formats (Windowless)");
files=1;
refdir = getDirectory("image"); name=getTitle(); 
dotIndex = indexOf(name, "."); 
name = substring(name, 0, dotIndex); 
reffilein=refdir+name+suffix;

rename("ref_source");
close();
run("Bio-Formats Importer", "open=["+reffilein+"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");
rename("ref_source");
run("Bio-Formats Importer", "open=["+reffilein+"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_2");
rename("ref_source_bg");
imageCalculator("Subtract stack" ,"ref_source","ref_source_bg");
selectWindow("ref_source_bg");
close();
rename("ref_source");
metadataref=getMetadata("Info");
metadatar=split(metadataref,"\n");
fre=0;
for (j=0;j<metadatar.length;j++) {
	if(startsWith(metadatar[j],"PARAMETERS: ACQUISITION SETTINGS - Frequency")==true) {
		sinfo=split(metadatar[j],"=");
		fre=sinfo[1];
		sinfo=split(fre," ");
		fre=sinfo[0];
		fre=parseFloat(fre);
	}
}
if (fre!=0) {
		hfreq=fre*1000000;
}
selectWindow("ref_source");
run("32-bit");
if (blur_ref>0) {
	run("Gaussian Blur...", "sigma="+blur_ref+" stack");
}
z=nSlices();
y=getHeight();
x=getWidth();
nt=z;
w=hfreq*2*3.1415927/1000000000;
wt=w*tauref;
phicor=atan(wt);
modcor=sqrt(wt*wt+1);
phasecorr=phasecorr0*3.1415927/180-phicor;
modcorr=modcorr0*modcor;
rename("input-bg-j");
run("Duplicate...", "duplicate");
rename("input-bg-jcos");
run("Duplicate...", "duplicate");
rename("input-bg-jsin");
for (i=1;i<=nt;i++) {
	s=sin(2*3.1415927*(i-1)/nt);
	c=cos(2*3.1415927*(i-1)/nt);
	selectWindow("input-bg-jsin");
	setSlice(i);
	run("Multiply...", "value="+s+" slice");		
	selectWindow("input-bg-jcos");
	setSlice(i);
	run("Multiply...", "value="+c+" slice");	
}
selectWindow("input-bg-jsin");
run("Z Project...", "projection=[Sum Slices]");
rename("sinsum");
run("Divide...", "value="+nt);	
selectWindow("input-bg-jcos");
run("Z Project...", "projection=[Sum Slices]");
rename("cossum");
run("Divide...", "value="+nt);	
selectWindow("input-bg-j");
run("Z Project...", "projection=[Sum Slices]");
rename("sum");
run("Divide...", "value="+nt);	
selectWindow("input-bg-jsin");
close();
selectWindow("input-bg-jcos");
close();
selectWindow("input-bg-j");
close();

imageCalculator("Divide create 32-bit", "cossum","sinsum");
rename("tau_phi");
selectWindow("tau_phi");
run("Macro...", "code=v=atan(v)");
run("Add...", "value="+phasecorr);	
selectWindow("tau_phi");
rename("phi_ref");
if (blur_ref>0) {
	run("Gaussian Blur...", "sigma="+blur_ref);
}
imageCalculator("Multiply create 32-bit", "sinsum","sinsum");
rename("sin2");
imageCalculator("Multiply create 32-bit", "cossum","cossum");
rename("cos2");
imageCalculator("Add create 32-bit", "sin2","cos2");
rename("sin2cos2");
run("Square Root");
imageCalculator("Divide create 32-bit" ,"sin2cos2","sum");
rename("tau_mod");
selectWindow("tau_mod");
m=2*modcorr;  
run("Multiply...", "value="+m);	
selectWindow("tau_mod");
rename("mod_ref");
if (blur_ref>0) {
	run("Gaussian Blur...", "sigma="+blur_ref);
}

selectWindow("sum");
close();
selectWindow("sinsum");
close();
selectWindow("cossum");
close();
selectWindow("sin2");
close();
selectWindow("cos2");
close();
selectWindow("sin2cos2");
close();

if (choice_open=="Convert entire directory") {
//	open();

	waitForUser("Please select one (sample) FLIM stack in the directory that you want to convert");
	run("Bio-Formats (Windowless)");
	filedir = getDirectory("image"); 
	close();
	list2=getFileList(filedir);
	list=list2;
	nfiles=list.length;
	files=0;
	for (ifile=0;ifile<nfiles;ifile++) {
		test2=0;
		string=list[ifile];
		test2=endsWith(string, suffix);
		if (test2==1) {
			list[files]=list2[ifile];
			files=files+1;		
		}
	}
}

if (choice_open=="Convert file") {
//	open();
	waitForUser("Please select the (sample) FLIM stack that you want to convert");
	run("Bio-Formats (Windowless)");
	filedir = getDirectory("image"); name=getTitle(); 
	dotIndex = indexOf(name, "."); 
	name = substring(name, 0, dotIndex); 
	fname=filedir+name+suffix;
	close();
	run("Bio-Formats Importer", "open=["+fname+"]  autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");
	filein=filedir+name+suffix;fileout=filedir+name+extension;

	files=1;
	

	filedate=File.dateLastModified(filein);
	filename=name+suffix;
//	print(filedate);
	rename("input");
	run("Bio-Formats Importer", "open=["+fname+"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_2");
	rename("sample_source_bg");
	imageCalculator("Subtract stack" ,"input","sample_source_bg");
	selectWindow("sample_source_bg");
	close();

}




w=w/1000;
for (ifile=0;ifile<files;ifile++) {
	if (files>1) {
//		open(filedir+list[ifile]);
//		run("Bio-Formats (Windowless)", "open=["+filedir+list[ifile]+"]");
		run("Bio-Formats Importer", "open=["+filedir+list[ifile]+"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");
		fname=getTitle;
		dotIndex = indexOf(fname, "."); 
           	 	fname = substring(fname, 0, dotIndex); 
		fileout=filedir+fname+extension;
		filein=filedir+fname+suffix;
		


		filedate=File.dateLastModified(filein);
		filename=fname+suffix;
//		print(filedate);
		rename("input");
		run("Bio-Formats Importer", "open=["+filedir+list[ifile]+"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_2");
		rename("sample_source_bg");
		imageCalculator("Subtract stack" ,"input","sample_source_bg");
		selectWindow("sample_source_bg");
		close();
	
	}
	


	selectWindow("input");
	metadata2=getMetadata("Info");
	metadata3="";
	metadata=split(metadata2,"\n");
	z=nSlices();
	y=getHeight();
	x=getWidth();
	for (j=0;j<metadata.length;j++) {
		if(startsWith(metadata[j],"PARAMETERS: ACQUISITION SETTINGS - NumPhases")==true) {
			sinfo=split(metadata[j],"=");
			nt=sinfo[1];
			nt=parseFloat(nt);
			numstacks=z/nt;
			numstacks_1=numstacks-1;
		}
	}
	if (numstacks>1) {
		for (i=0;i<numstacks;i++){
			for (j=0;j<metadata.length;j++) {
				if(startsWith(metadata[j],"FLIMIMAGE: TIMESTAMPS - t"+i+" =")==true) {
					sinfo=split(metadata[j],"=");
					sinfo1=sinfo[1];
					sinfo=split(sinfo1," ");
					sinfo00=sinfo[0];
					sinfo01=sinfo[1];
					sinfo00=parseInt(sinfo00);
					sinfo01=parseInt(sinfo01);
					if (i==0){ 
						s0=sinfo00;
						s1=sinfo01;
					}
				}
			}
			time0=sinfo00-s0;
			time1=sinfo01-s1;
			time=time0*4294967296+time1;
			time=time/10000;
			time=round(time);
			time=time/1000;
//			print("Time for image "+i+" is "+time+" s");
			timepoints[i]=time;
		}	
	}
	
	test=1;
	if(test==0) {
		close();
	}else{
		run("32-bit");
		if (blur_choice>0) {
			run("Gaussian Blur...", "sigma="+blur_choice+" stack");
		}
		if (bg_choice=="rolling ball with 500 pixels diameter") run("Subtract Background...", "rolling=500 stack");
		if (bg_choice=="rolling ball with 250 pixels diameter") run("Subtract Background...", "rolling=250 stack");
		if (bg_choice=="rolling ball with 125 pixels diameter") run("Subtract Background...", "rolling=125 stack");
		if (bg_choice=="sliding parabola with 500 pixels diameter") run("Subtract Background...", "rolling=500 sliding stack");
		if (bg_choice=="sliding parabola with 250 pixels diameter") run("Subtract Background...", "rolling=250 sliding stack");
		if (bg_choice=="sliding parabola with 125 pixels diameter") run("Subtract Background...", "rolling=125 sliding stack");
		rename("input-bg");
//		numstacks=1;
//		nt=z/numstacks;
		ot=3*numstacks;
		//print(nt,z,numstacks,ot);
		f=hfreq/7510000000*3.1415927/180;
		newImage("Output", "16-bit", x, y, ot);
		if (numstacks>1) {
			selectWindow("input-bg");
			run("Stack to Hyperstack...", "order=xyczt(default) channels=1 slices="+nt+" frames="+numstacks+" display=Color");
			rename("input-bg");
		}
		for (j=1;j<=numstacks;j++) {
			slice=(j-1)*nt+1;
			selectWindow("input-bg");
			if (numstacks>1)  {	
				setSlice(slice);	
				run("Reduce Dimensionality...", "slices keep");
			}
			rename("input-bg-j");
			run("Duplicate...", "duplicate");
			rename("input-bg-jcos");
			run("Duplicate...", "duplicate");
			rename("input-bg-jsin");
			for (i=1;i<=nt;i++) {
				s=sin(2*3.1415927*(i-1)/nt);
				c=cos(2*3.1415927*(i-1)/nt);
				selectWindow("input-bg-jsin");
				setSlice(i);
				run("Multiply...", "value="+s+" slice");		
				selectWindow("input-bg-jcos");
				setSlice(i);
				run("Multiply...", "value="+c+" slice");	
				}
			selectWindow("input-bg-jsin");
			run("Z Project...", "projection=[Sum Slices]");
			rename("sinsum");
			selectWindow("input-bg-jcos");
			run("Z Project...", "projection=[Sum Slices]");
			rename("cossum");
			selectWindow("input-bg-j");
			run("Z Project...", "projection=[Sum Slices]");
			rename("sum");
			selectWindow("input-bg-j");
			close();
			selectWindow("input-bg-jsin");
			close();
			selectWindow("input-bg-jcos");
			close();
			oslice=(j-1)*3+1;
			selectWindow("Output");
			setSlice(oslice);	
			selectWindow("sum");
			if (stat==1){
				// create threshold image
				run("Duplicate...", " ");
				rename("blai");
				run("Divide...", "value="+nt);	
				getStatistics(area,mean,min,max);
				upthres=max+1;
				imin=bg_stat;
				setMinAndMax(imin, upthres);
				run("8-bit");
				setThreshold(1, 254);
				run("Threshold", "thresholded remaining black slice");
				run("Invert");
				run("Divide...", "value=255");
				run("16-bit");
				getStatistics(npi,tau_div);
				tau_area=npi*tau_div;
				rename("threshold2");
			}
			selectWindow("sum");
			run("Duplicate...", " ");
			rename("sumcopy");
			run("Divide...", "value="+nt);	
			setMinAndMax(0,65535);
			run("16-bit");

			run("Copy");
			selectWindow("Output");
			run("Paste");
			selectWindow("sumcopy");
			close();
	
			oslice=oslice+1;
			setSlice(oslice);	
			imageCalculator("Divide create 32-bit", "cossum","sinsum");
			rename("tau_phi");
			selectWindow("tau_phi");
			run("Macro...", "code=v=atan(v)");
			imageCalculator("Subtract 32-bit", "tau_phi", "phi_ref");
			run("Macro...", "code=v=tan(v)");
			run("Divide...", "value="+w);	
			run("Duplicate...", " ");
			rename("tau_phi2");
			setMinAndMax(0,65535);
			run("16-bit");		

			if (stat==1){
				imageCalculator("Multiply create 32-bit", "tau_phi2","threshold2");
				rename("blap");
				setMinAndMax(0,65535);
				run("16-bit");		
				getStatistics(npix,taup_mean,pllp,qllq,taup_sd);
				x2=npix*(taup_sd*taup_sd+taup_mean*taup_mean);
				taup_mean=taup_mean/tau_div;
				taup_sd=sqrt((x2-tau_area*taup_mean*taup_mean)/tau_area);
				taup_mean=round(taup_mean);
				taup_sd=round(taup_sd);
				taup_mean=taup_mean/1000; taup_sd=taup_sd/1000;
				selectWindow("blap");
				close();
			}
			selectWindow("tau_phi2");
			run("Copy");
			selectWindow("Output");
			run("Paste");
			oslice=oslice+1;
			setSlice(oslice);	
			imageCalculator("Multiply create 32-bit", "sinsum","sinsum");
			rename("sin2");
			imageCalculator("Multiply create 32-bit", "cossum","cossum");
			rename("cos2");
			imageCalculator("Add create 32-bit", "sin2","cos2");
			rename("sin2cos2");
			run("Square Root");
			imageCalculator("Divide create 32-bit","sum", "sin2cos2");
			rename("tau_mod");
			selectWindow("tau_mod");
			m=2;
			run("Divide...", "value="+m);	
			imageCalculator("multiply create 32-bit", "mod_ref","tau_mod");
			run("Square");
			run("Subtract...", "value=1");	
			run("Divide...", "value="+m);	

			rename("modje1");
			run("Duplicate...", " ");
			run("Abs");
			rename("modje2");
			selectWindow("tau_mod");
			close();
			imageCalculator("add create 32-bit", "modje1","modje2");
			rename("tau_mod");
			selectWindow("modje1");
			close();
			selectWindow("modje2");
			close();
			run("Square Root");
			run("Divide...", "value="+w);	
			run("Duplicate...", " ");
			rename("tau_mod2");
			setMinAndMax(0,65535);
			run("16-bit");		

			if (stat==1){
				imageCalculator("Multiply create 32-bit", "tau_mod2","threshold2");
				rename("blam");	
				setMinAndMax(0,65535);
				run("16-bit");		
				getStatistics(npix2,taum_mean,pl,ql,taum_sd);
				x2=npix2*(taum_sd*taum_sd+taum_mean*taum_mean);
				taum_mean=taum_mean/tau_div;
				taum_sd=sqrt((x2-tau_area*taum_mean*taum_mean)/tau_area);
				taum_mean=round(taum_mean);
				taum_sd=round(taum_sd);
				taum_mean=taum_mean/1000; taum_sd=taum_sd/1000;
				kk=j-1;
				time=timepoints[kk];
				metadata3=metadata3+"Stack # "+j+"Time_after_start "+time+":Average tau_phi= "+taup_mean+" +- "+taup_sd+" ns, Average tau_mod= "+taum_mean+" +- "+taum_sd+" ns \n";
				selectWindow("blam");
				close();
				selectWindow("threshold2");
				close();
			}
			selectWindow("tau_mod2");
			run("Copy");
			selectWindow("Output");
			run("Paste");
			setMinAndMax(0,65535);
			run("16-bit");
			selectWindow("sum");
			close();
			selectWindow("sinsum");
			close();
			selectWindow("cossum");
			close();
			selectWindow("tau_phi");
			close();
			selectWindow("tau_phi2");
			close();
			selectWindow("tau_mod");
			close();			
			selectWindow("tau_mod2");
			close();
			selectWindow("sin2");
			close();
			selectWindow("cos2");
			close();
			selectWindow("sin2cos2");
			close();
			if(log_choice==1){
				if(logje==0){
  					print("  ");
					print("Conversion start at: "+dm+"-"+month+"-"+year+"  "+hr+":"+mi+":"+sec);
					print("Directory for sample files: "+filedir);
					print("Reference file: "+reffilein);
					print("--------------Sample statistics--------------");
					print("filename \tfile creation date\tstack#\ttime (s)\ttau(phi) (ns)\tsd_tau(phi) (ns)\ttau(mod) (ns)\tsd_tau(mod) (ns)\tIntensity threshold\tnumber of phase images\tH_freq (Hz)\tReference tau (ns)\tPhasecorr (deg)\tModcorr (fold)"); 
					logje=1;
				}
				kk=j-1;
				time=timepoints[kk];	
				print(filename+"\t"+filedate+"\t"+j+"\t"+time+"\t"+taup_mean+"\t"+taup_sd+"\t"+taum_mean+"\t"+ taum_sd+"\t"+ bg_stat+"\t"+nt+"\t"+ hfreq+"\t"+ tauref+"\t"+ phasecorr0+"\t"+ modcorr0);
			}

		}

		selectWindow("Output");
		fre=hfreq/1000000;
		run("Rotate 90 Degrees Right");
		metadata2=metadata2+"------------------------Original FLIM-II metadata-------------------------\n";
		metadata2=metadata2+"Original FLIM stack name: "+fname+"\n";
		metadata2=metadata2+"Original FLIM stack creation date: "+filedate+"\n";
		metadata2=metadata2+"Number of stacks or time points: "+numstacks+" \n";
		metadata2=metadata2+"Number of phase images per FLIM recording: "+nt+" \n";
		metadata2=metadata2+"Homodyne modulation frequency (MHz) ="+fre+" \n";
		metadata2=metadata2+"----------------Metadata used for conversion and statistics-----------------\n";
		metadata2=metadata2+"File converted using Convert_flim-II-to-ij ImageJ macro \n";
		metadata2=metadata2+"Reference file: "+reffilein+" \n";
		metadata2=metadata2+"Tauref: "+tauref+" \n";
		metadata2=metadata2+"Additional phase correction (degrees): "+phasecorr0+" \n";
		metadata2=metadata2+"Additional modulation correction (fold factor): "+modcorr0+" \n";
		metadata2=metadata2+"Gaussian blur size of reference data: "+blur_ref+" \n";	
		metadata2=metadata2+"Gaussian blur size of sample data: "+blur_choice+" \n";
		metadata2=metadata2+"Background_choice: "+bg_choice+" \n";
		metadata2=metadata2+"Statistics_choice: "+stat+" \n";
		metadata2=metadata2+"Background int for stats: "+bg_stat+" \n";
		if (numstacks>1){
			metadata2=metadata2+"--------------------------Timing information----------------------------\n";
			for (j=1;j<=numstacks;j++) {
				kk=j-1;
				metadata2=metadata2+"Time_after_start "+kk+" ="+timepoints[kk]+" \n";
			}
		}
		metadata2=metadata2+"------------------------Average image statistics--------------------------\n";
		metadata2=metadata2+metadata3;
		metadata2=metadata2+"-----------------------------------------------------------------------\n";
		


		setMetadata("Info",metadata2);
		saveAs("Tiff", fileout);
		close();
		if (numstacks>1)  {
			selectWindow("input-bg");	
			close();	
	
		}
	}
}
selectWindow("mod_ref");
close();
selectWindow("phi_ref");
close();
setBatchMode("exit and display");




  
   
  


