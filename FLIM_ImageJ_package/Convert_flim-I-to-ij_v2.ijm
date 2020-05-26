// Macro to convert FLIM files from FLIM-I for imageJ version 2
//===========================================================================
// Version history:
// Version 2   : 28-02-2019 Dorus Gadella, Adjustments and synchronization with flim-II converter
// Version 1.2: 29-10-2018 Dorus Gadella, Adjusted file import settings to be compatible with Bioformats & Fiji
// Version 1.1: 29-10-2018 Dorus Gadella, Adjusted file import settings to be compatible with Bioformats & Fiji
// Version 1   : 28-10-2018 Dorus Gadella, Including exporting metadata and automated conversion of all files in a directory.
//===========================================================================
// Version 0: 26-10-2018 Created by Dorus Gadella 26-10-2018
// Reads in metadata from FLIM-I FLIM stacks (.ics files) and converts them for input in ImageJ.
// Reference phase and modulation and frequency are imported, from which the phase and modulation images are calculated.
// Includes iris correction for 2x2 binned full frame FLIM stacks and 2x2 binned 256x256 pixel frame FLIM stacks.
// Optional can background subtraction per phase image. This can be a signed 16-bit correction or a rolling ball correction.
// Includes calculation of average image lifetime statistics based on a fixed background threshold.
//===========================================================================


Dialog.create("Convert FLIM-I ics stack to input for ImageJ Lifetimes15 macro");
Dialog.addChoice("What to convert ",newArray("Convert file","Convert current image","Convert entire directory"),"Convert file");
Dialog.addChoice("Background subtraction in phase images ", newArray("no background subtraction","subtract 32768","rolling ball with 500 pixels diameter","rolling ball with 250 pixels diameter","rolling ball with 125 pixels diameter", "sliding parabola with 500 pixels diameter","sliding parabola with 250 pixels diameter","sliding parabola with 125 pixels diameter"),"no background subtraction");
Dialog.addNumber("Gaussian blur data, 0=no, >0 is blur size", 1); //was1
Dialog.addCheckbox("Perform iris correction on phase lifetime:", true);
Dialog.addCheckbox("Calculate average image statistics:", true);
Dialog.addNumber("Background threshold for image stats", 100); //was 100
Dialog.addCheckbox("Write logfile with metadata and statistics:", true);
Dialog.addString("Output extension of filename", "_ij.tif");
Dialog.show();
choice_open= Dialog.getChoice();
bg_choice= Dialog.getChoice();
blur_choice = Dialog.getNumber();
iris_choice= Dialog.getCheckbox();
stat= Dialog.getCheckbox();
bg_stat= Dialog.getNumber();
log_choice= Dialog.getCheckbox();
extension = Dialog.getString();
//end input

//irislocation_PC="C:\\Program Files\\ImageJ\\FLIM macro\\";
irislocation_MAC="/Users/Franka/Desktop/FLIM_ImageJ_package/";
//irislocation_Parallels="C:\\Users\\gadella\\Desktop\\Image Proc\\";
irislocation=irislocation_MAC;
irisfile=irislocation+"new_iris.tif";
irisfile2=irislocation+"new_iris256x256.tif;
phasecorr=0;
modcorr=1;
logje=0;
v=newArray(160);
timepoints=newArray(10000);
timepoints[0]=0;
suffix=".ics";
getDateAndTime(year,month,dw,dm,hr,mi,sec,msec);
month=month+1;
//==============
setBatchMode(true);

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
	files=1;
	filedir = getDirectory("image"); filein=filedir+getTitle(); fileint=filein+".ics";fname=getTitle()+".ics";fileout=filein+extension;
	test2=endsWith(filein, suffix);
	if (test2==1) {
		fileint=filein;
		fname=getTitle;
		dotIndex = indexOf(filein, "."); 
           	 	filein = substring(filein, 0, dotIndex); 
		fileout=filein+extension;
	}


	filedate=File.dateLastModified(fileint);
//	print(filedate);
	rename("input");
}
if (choice_open=="Convert current image") {
	files=1;
	filedir = getDirectory("image"); filein=filedir+getTitle(); fileint=filein+".ics";fname=getTitle()+".ics";fileout=filein+extension;
	test2=endsWith(filein, suffix);
	if (test2==1) {
		fileint=filein;
		fname=getTitle;
		dotIndex = indexOf(filein, "."); 
           	 	filein = substring(filein, 0, dotIndex); 
		fileout=filein+extension;
	}

	filedate=File.dateLastModified(fileint);
//	print(filedate);
	rename("input");
}
open(irisfile);

rename("iris_source");
selectWindow("iris_source");
run("32-bit");
open(irisfile);

rename("iris_source2");
selectWindow("iris_source2");
run("32-bit");
for (ifile=0;ifile<files;ifile++) {
	if (files>1) {
//		open(filedir+list[ifile]);
		run("Bio-Formats (Windowless)", "open=["+filedir+list[ifile]+"]");
		filein=filedir+getTitle(); fileint=filein+".ics";fname=getTitle()+".ics";fileout=filein+"_ij.tiff";
		test2=endsWith(filein, suffix);
			if (test2==1) {
			fileint=filein;
			fname=getTitle;
			dotIndex = indexOf(filein, "."); 
           	 		filein = substring(filein, 0, dotIndex); 
			fileout=filein+extension;
		}

		filedate=File.dateLastModified(fileint);
//		print(filedate);
		rename("input");
	
	}
	
	selectWindow("input");
	metadata2=getMetadata("Info");
	metadata3="";

	z=nSlices();
	y=getHeight();
	x=getWidth();


	xx=0;yy=0;
	for (i=0; i<160; i++) {
		xx=floor(i/y);
		yy=i-xx*y;
		v[i]=getPixel(xx,yy);
	}
		test=1;
	phase=retrieve_phase(phase);
	if (phase<0) test=0;
	if (phase>360) test=0;
	//print("Calibration phase= "+phase+" degrees");
	mod=retrieve_mod(mod);
	if (mod<0) test=0;
	if (mod>360) test=0;
	//print("Calibration modulation="+mod+" %");
	bleachflag=retrieve_bleachflag(bleachflag);
	//print("Bleachflag="+bleachflag);
	orderflag=retrieve_orderflag(orderflag);
	//print("Orderflag="+orderflag);
	exptime=retrieve_exptime(exptime);
	//print("Exptime="+exptime+" ms");
	laser=retrieve_laser(laser);
	if (laser<0) test=0;
	if (laser>10) test=0;
	//print("Laser="+laser);
	lens=retrieve_lens(lens);
	if (lens<0) test=0;
	if (lens>10) test=0;
	//print("Lens="+lens);
	lfreq=retrieve_lfreq(lfreq);
	//print("Lfreq="+lfreq);
	hfreq=retrieve_hfreq(hfreq);
	//print("Hfreq="+hfreq);
	numstacks=retrieve_numstacks(numstacks);
	deltatime=retrieve_deltatime(deltatime);
	//print("Numstacks="+numstacks);
	//rename("input-bg");
	for (i=0;i<numstacks;i++){
		timepoints[i]=i*deltatime;
	}
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
		if (bg_choice=="subtract 32768") run("Subtract...", "value=32768 stack");
		rename("input-bg");
		nt=z/numstacks;
		ot=3*numstacks;
		//print(nt,z,numstacks,ot);
		w=hfreq*2*3.1415927/1000000000000;
		p=phase*2*3.1415927/360;
		m=200/mod;
		f=hfreq/7510000000*3.1415927/180;
		selectWindow("iris_source");
		if (x==256) selectWindow("iris_source2");
		run("Duplicate...", " ");
		rename("iris");
		run("Multiply...", "value="+f);
		newImage("Output", "32-bit", x, y, ot);
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
			run("Copy");
			selectWindow("Output");
			run("Paste");
			run("Divide...", "value="+nt+" slice");	
			oslice=oslice+1;
			setSlice(oslice);	
			imageCalculator("Divide create 32-bit", "cossum","sinsum");
			rename("tau_phi");
			selectWindow("tau_phi");
			run("Macro...", "code=v=atan(v)");
			//run("Add...", "value=3.1415927");	
			p=(phase-phasecorr)*2*3.1415927/360;
			run("Subtract...", "value="+p);	
			if (iris_choice==1){
				imageCalculator("Subtract 32-bit", "tau_phi","iris");
			}
			run("Macro...", "code=v=tan(v)");
			run("Divide...", "value="+w);	
			selectWindow("tau_phi");
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
			imageCalculator("Divide create 32-bit", "sum","sin2cos2");
			rename("tau_mod");
			selectWindow("tau_mod");
			m=200/(mod*modcorr);  
			run("Divide...", "value="+m);	
			run("Square");
			run("Subtract...", "value=1");	
			run("Square Root");
			run("Divide...", "value="+w);	
			selectWindow("tau_mod");
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
			selectWindow("tau_mod");
			close();
			selectWindow("tau_phi2");
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
					print("Directory: "+filedir+"\tConversion start at: "+dm+"-"+month+"-"+year+"  "+hr+":"+mi+":"+sec);
					print("filename \tfile creation date\tstack#\ttime (s)\texposure time (ms) \ttau(phi) (ns)\tsd_tau(phi) (ns)\ttau(mod) (ns)\tsd_tau(mod) (ns)\tIntensity threshold\tCalibration phase (deg)\tCalibration modulation (%)\tnumber of phase images\tbleachflag\torderflag\tlaser\tlens\tLfreq (Hz)\tH_freq (Hz)\tPhasecorr (deg)\tModcorr (fold)"); 
					logje=1;
				}
				kk=j-1;
				time=timepoints[kk];	
				print(fname+"\t"+filedate+"\t"+j+"\t"+time+"\t"+exptime+"\t"+taup_mean+"\t"+taup_sd+"\t"+taum_mean+"\t"+ taum_sd+"\t"+ bg_stat+"\t"+phase+"\t"+ mod+"\t"+ nt+"\t"+bleachflag+"\t"+ orderflag+"\t"+ laser+"\t"+ lens+"\t"+ lfreq+"\t"+ hfreq+"\t"+ phasecorr+"\t"+ modcorr);
			}

		}
		selectWindow("iris");
		close();
		selectWindow("Output");
		fre=hfreq/1000000;
		metadata2=metadata2+"------------------------Original FLIM-I metadata-------------------------\n";
		metadata2=metadata2+"Original FLIM stack name: "+fname+"\n";
		metadata2=metadata2+"Original FLIM stack creation date: "+filedate+"\n";
		metadata2=metadata2+"Calibration phase: "+phase+" degrees\n";
		metadata2=metadata2+"Calibration modulation: "+mod+" %\n";
		metadata2=metadata2+"Number of stacks or time points: "+numstacks+" \n";
		metadata2=metadata2+"Number of phase images per FLIM recording: "+nt+" \n";
		metadata2=metadata2+"Homodyne modulation frequency (MHz) ="+fre+" \n";
		metadata2=metadata2+"Bleachflag:  "+bleachflag+"\n";
		metadata2=metadata2+"Orderflag:  "+orderflag+"\n";
		metadata2=metadata2+"Exptime:  "+exptime+" ms\n";	
		if(numstacks>1) {
			metadata2=metadata2+"Deltatime:  "+deltatime+" s\n";	
		}	
		metadata2=metadata2+"Laser: "+laser+"\n";
		metadata2=metadata2+"Lens:  "+lens+"\n";
		metadata2=metadata2+"Lfreq:  "+lfreq+" Hz\n";
		metadata2=metadata2+"Hfreq:  "+hfreq+" Hz\n";
		metadata2=metadata2+"----------------Metadata used for conversion and statistics-----------------\n";
		metadata2=metadata2+"File converted using Convert_flim-I-to-ij_v2 ImageJ macro \n";
		metadata2=metadata2+"Iris_choice: "+iris_choice+" \n";
		metadata2=metadata2+"Additional phase correction (degrees): "+phasecorr+" \n";
		metadata2=metadata2+"Additional modulation correction (fold factor): "+modcorr+" \n";
		metadata2=metadata2+"Gaussian blur size of sample data: "+blur_choice+" \n";
		metadata2=metadata2+"Background_choice: "+bg_choice+" \n";
		metadata2=metadata2+"Statistics_choice: "+stat+" \n";
		metadata2=metadata2+"Background int for stats: "+bg_stat+" \n";
		if(numstacks>1) {
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
selectWindow("iris_source");	
close();
selectWindow("iris_source2");
close();
setBatchMode("exit and display");

function retrieve_phase(phase) {
k=0;
for (i=0; i<16; i++) {
	j=v[i];
	k=k+pow(2,i)*2*(j/2-floor(j/2));
	}
phase=(k-0.5)*360/65535;
return phase;
}

function retrieve_mod(mod) {
k=0;
for (i=0; i<16; i++) {
	j=v[i+16];
	k=k+pow(2,i)*2*(j/2-floor(j/2));
	}
mod=(k-0.5)*200/65535;
return mod;
}
function retrieve_bleachflag(bleachflag) {
k=0;
	j=v[32];
	k=k+2*(j/2-floor(j/2));
	bleachflag=k;
return bleachflag;
}
function retrieve_orderflag(orderflag) {
k=0;
	j=v[33];
	k=k+2*(j/2-floor(j/2));
	orderflag=k;
return orderflag;
}
function retrieve_exptime(exptime) {
k=0;
for (i=0; i<16; i++) {
	j=v[i+34];
	k=k+pow(2,i)*2*(j/2-floor(j/2));
	}
exptime=k;
return exptime;
}
function retrieve_deltatime(deltatime) {
k=0;
for (i=0; i<16; i++) {
	j=v[i+50];
	k=k+pow(2,i)*2*(j/2-floor(j/2));
	}
deltatime=k;
return deltatime;
}

function retrieve_laser(laser) {
k=0;
for (i=0; i<16; i++) {
	j=v[i+66];
	k=k+pow(2,i)*2*(j/2-floor(j/2));
	}
laser=k;
return laser;
}

function retrieve_lens(lens) {

k=0;
for (i=0; i<16; i++) {
	j=v[i+82];
	k=k+pow(2,i)*2*(j/2-floor(j/2));
	}
lens=k;
return lens;
}


function retrieve_lfreq(lfreq) {
k=0;
for (i=0; i<16; i++) {
	j=v[i+98];
	k=k+pow(2,i)*2*(j/2-floor(j/2));
	}
lfreq=k;
lfreq=lfreq*10000;
return lfreq;
}

function retrieve_hfreq(hfreq) {
k=0;
for (i=0; i<16; i++) {
	j=v[i+114];
	k=k+pow(2,i)*2*(j/2-floor(j/2));
	}
hfreq=k;
hfreq=hfreq*10000;
return hfreq;
}

function retrieve_numstacks(numstacks) {
k=0;
for (i=0; i<16; i++) {
	j=v[i+130];
	k=k+pow(2,i)*2*(j/2-floor(j/2));
	}
numstacks=k;
return numstacks;
}
 //    if (bleachflag==1) *(image+32)=*(image+32) | 1;
  //   if ((n_deltatime & (unsigned short)pow(2,i))>0) *(image+50+i)=*(image+50+i) | 1;
  // if ((n_exptime & (unsigned short)pow(2,i))>0) *(image+34+i)=*(image+34+i) | 1;
   //   if ((n_hfreq & (unsigned short)pow(2,i))>0) *(image+114+i)=*(image+114+i) | 1;hfreq=hfreq/10000;
   // if ((n_laser & (unsigned short)pow(2,i))>0) *(image+66+i)=*(image+66+i) | 1;
   //if ((n_lens & (unsigned short)pow(2,i))>0) *(image+82+i)=*(image+82+i) | 1;
    // if ((n_lfreq & (unsigned short)pow(2,i))>0) *(image+98+i)=*(image+98+i) | 1;lfreq=lfreq/10000;
//  if ((n_stack & (unsigned short)pow(2,i))>0) *(image+130+i)=*(image+130+i) | 1;
   //    if (orderflag==1) *(image+33)=*(image+33) | 1;
  
   
  


