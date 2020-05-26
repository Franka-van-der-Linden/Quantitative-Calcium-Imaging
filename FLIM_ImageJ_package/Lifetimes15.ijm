// macro display coloured composite lifetime image with histograms
// Dorus Gadella 23-02-19 Lifetimes15.txt updated with processing of entire directory, FLIM-II intensity rescaling, automatic metadata import/export
// Dorus Gadella 07-07-15 Lifetimes14.txt updated with automated modal background estimation
// Dorus Gadella 28-06-15 Lifetimes13.txt updated with metadata from FLIM screen
// Dorus Gadella 28-01-14 Lifetimes12.txt updated with optional output stack of lifetime images
// Dorus Gadella 18-12-13 Lifetimes11.txt updated with better time-lapse thresholding and timing display
// Dorus Gadella 21-06-13 Lifetimes10c.txt updated with optional average tau as line in histogram
// Dorus Gadella 21-06-13 Lifetimes10b.txt updated with optional tau(phi) or tau(mod) only display
// Dorus Gadella 20-06-13 Lifetimes10a.txt updated with colortable for intensity image
// Dorus Gadella 12-10-09 Lifetimes10.txt updated for white/black background & fonts
// Dorus Gadella 13-08-08 Lifetimes9.txt updated for use of many colortables
// Dorus Gadella 12-08-08 Lifetimes8.txt updated for polar plot
// Dorus Gadella 11-08-08 Lifetimes7.txt updated for fire-colortable intensity-lifetime maps
// Dorus Gadella 22-11-06 Lifetimes6.txt updated for FLIM time series
// Dorus Gadella 10-07-06 Lifetimes5.txt included 2D histograms
// Dorus Gadella 05-07-06 Lifetimes4.txt updated thresholding
// Dorus Gadella 01-07-06 Lifetimes3.txt added colored lifetime histograms
// Dorus Gadella 29-06-06 Lifetimes2.txt updated thresholding
// Dorus Gadella 29-06-06 Lifetimes.txt


//-------------------------------------------------------------------------------
//==========================================================
//open file & input parameters for display

Dialog.create("Input for Lifetime Display");
Dialog.addChoice("Open new file or work on current image",newArray("Open new file","Work on current image","Open entire directory"),"Open new file");
Dialog.addNumber("Minimal Tau for display (ns)  :", 0);
Dialog.addNumber("Maximal Tau for display (ns) :", 5);
Dialog.addNumber("Low threshold_for_Intensity  :", 1000);
Dialog.addNumber("High threshold_for_Intensity :", 250);
Dialog.addNumber("Gamma_for_Intensity            :", 0.7);
//Dialog.addNumber("Modulation frequency (MHz)  :", 40);
Dialog.addCheckbox("Threshold intensity image:", false);
Dialog.addCheckbox("Exclude analysis for I>Ihigh:", false);
Dialog.addCheckbox("Automatic determination of Ihigh:", true);
Dialog.addCheckbox("Include 2D-histograms:", true);
Dialog.addCheckbox("Display average lifetime as line in histogram:", false);
Dialog.addCheckbox("Make additional output image stack:", false);
//Dialog.addNumber("Time between FLIM stacks  :", 5);
//Dialog.addChoice("Time lapse units:",newArray("sec","min","hour"),"sec");
Dialog.addChoice("Display only Tau(phi), only Tau(mod), or both:",newArray("Both", "Tau(phi)","Tau(Mod)"),"Both");
Dialog.addChoice("Colortable for Lifetime Images:",newArray("Fire","Grays","Ice","Rainbow RGB","Red Hot", "Royal","16_colors","Green Fire Blue", "Phase"),"Fire");
Dialog.addChoice("Colortable for Intensity Images:",newArray("Grays","Blue", "Cyan","Green","Yellow","Red"),"Grays");
Dialog.addChoice("Background:",newArray("Black","White"),"Black");
Dialog.addChoice("Foreground :",newArray("Yellow","White","Black"),"Yellow");

Dialog.show();
choice_open= Dialog.getChoice();
taumin = Dialog.getNumber();
taumax = Dialog.getNumber();
imin2 = Dialog.getNumber();
imax2=Dialog.getNumber();
igamma= Dialog.getNumber();
//freq= Dialog.getNumber();
thres_yes= Dialog.getCheckbox();
up_yes=Dialog.getCheckbox();
ihigh_yes=Dialog.getCheckbox();
his2D_yes=Dialog.getCheckbox();
line_yes=Dialog.getCheckbox();
extrastack=Dialog.getCheckbox();
//dtime=Dialog.getNumber();
//dtimeunit=Dialog.getChoice();
mod_or_phi=Dialog.getChoice();
colortable=Dialog.getChoice();
colortableint=Dialog.getChoice();
background=Dialog.getChoice();
foreground=Dialog.getChoice();
// end input

dtimeunit="s";
//==========================================================
// calculation of scaling factors
if (taumin<0) taumin=0;
if (taumax<taumin) taumax=taumin+1;
if (imin2<0) imin2=0;
if (imax2<imin2) imax2=100;
if (ihigh_yes==true) up_yes=false;
suffix="_ij.tif";
setBatchMode(true);
timepoints=newArray(10000);
timepoints[0]=0;
if (choice_open=="Open entire directory") {
//	open();
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
if (choice_open=="Open new file") {
	open();
//	run("Bio-Formats (Windowless)");
	files=1;
	filedir = getDirectory("image"); filein=filedir+getTitle();fileout=filein+"_out.tif";
	test2=endsWith(filein, suffix);
	if (test2==1) {
		fileint=filein;
		fname=getTitle;
		dotIndex = indexOf(filein, "."); 
           	 	filein = substring(filein, 0, dotIndex); 
		fileout=filein+"_out.tif";
		if (extrastack==true) fileout2=filein+"_2.tif";
	}


	filedate=File.dateLastModified(fileint);
//	print(filedate);
	rename("input");
}
if (choice_open=="Work on current image") {
	files=1;
	filedir = getDirectory("image"); filein=filedir+getTitle(); fileout=filein+"_out.tif";
	test2=endsWith(filein, suffix);
	if (test2==1) {
		fileint=filein;
		fname=getTitle;
		dotIndex = indexOf(filein, "."); 
           	 	filein = substring(filein, 0, dotIndex); 
		fileout=filein+"_out.tif";
		if (extrastack==true) fileout2=filein+"_2.tif";
	}

	filedate=File.dateLastModified(fileint);
//	print(filedate);
	rename("input");
}

for (ifile=0;ifile<files;ifile++) {
	if (files>1) {
		open(filedir+list[ifile]);
//		run("Bio-Formats (Windowless)", "open=["+filedir+list[ifile]+"]");
		filein=filedir+getTitle(); fileout=filein+"_out.tif";
		test2=endsWith(filein, suffix);
			if (test2==1) {
			fileint=filein;
			fname=getTitle;
			dotIndex = indexOf(filein, "."); 
           	 		filein = substring(filein, 0, dotIndex); 
			fileout=filein+"_out.tif";
			if (extrastack==true) fileout2=filein+"_out2.tif";
		}

		filedate=File.dateLastModified(fileint);
//		print(filedate);
		rename("input");
	
	}
	selectWindow("input");
	metadata2=getMetadata("Info");
//	print(metadata2);
	metadata=split(metadata2,"\n");
	for (qq=0;qq<metadata.length;qq++) {
		if(startsWith(metadata[qq],"Homodyne modulation frequency (MHz) =")==true) {
			sinfo=split(metadata[qq],"=");
			freq=sinfo[1];
			freq=parseFloat(freq);
		}
	}
	wfreq=freq*2*3.14159265/1000000;

//filedir = "C:\\Users\\Richard\\"; filein=filedir+getTitle(); fileint=filein+".ics";fname=getTitle()+".ics";fileout=filein+".tiff";
rename("input");
run(colortableint);
y=getHeight();
x=getWidth();
numstacks=nSlices()/3;
if (x<150) {
	rename("inputje");
	sscale=150/x; 
	run("Scale...", "x="+sscale+" y="+sscale+" interpolate process create title=input");
	selectWindow("inputje");
	close();
	selectWindow("input");
	y=getHeight();
	x=getWidth();
}
if (ihigh_yes==true) {
	selectWindow("input");
	imax2=0;
	for (i=1;i<=numstacks;i++) {
		j=3*(i-1)+1;
		setSlice(j);
		getStatistics(area,mean,min,max);
		if (max>imax2) imax2=max;
	}
	if(imax2>3000) {
// This has been added 23-2-2019 to accommodate 16-bit data from FLIM-II
		for (i=1;i<=numstacks;i++) {
			j=3*(i-1)+1;
			setSlice(j);
			run("Divide...", "value=16 slice");
		}
		imax2=imax2/16;
		imin=imin2/16;
	}else{
		imin=imin2;
	}
}
fs=x*0.115; fs=round(fs);
run("Set Measurements...", "decimal=9");
x_1=x-1;y_1=y-1;
fss=fs/2; fss=round(fss); fss_1=fss-1;
hx=y;
if (x<y) hx=x;
hy=(12+2*y+2*fs)/3;
if(hy<hx) hx=hy;
hy=hx;ihx=hx-fs;ihy=hy-fs;hscale=ihy/255;
ihx_1=ihx-1;ihy_1=ihy-1;
ntau=taumax-taumin;
nint=imax2-imin;
ni=1;
fspace=2.2*fss/3;
if (imax2>=10) fspace=4.4*fss/3;
if (imax2>=100) fspace=2.2*fss;
if (imax2>=1000) fspace=2.2*fss*1.33;
if (nint>10) ni= 5;
if (nint>25) ni= 10;
if (nint>100) ni=50;
if (nint>250) ni=100;
if(nint>1000) ni=500;
if (nint>2500) ni=1000;
ji=imax2/ni; ji=round(ji);
nt=0.1;ntt=0.5;
if (ntau>1) nt=0.5;
if (ntau>2.5) {
	nt=1;
	ntt=1;
}
jt=taumax/nt; jt=round(jt);jtt=taumax/ntt;jtt=round(jtt);
ox=4*x+70+2*hx+2*fspace;
if (his2D_yes==false) ox=4*x+50;
if(mod_or_phi!="Both")  ox=3*x+40;
oy=2*y+34+2*fs;
ooy=2.5*fs+10+y+10*fss+0.48*x;
oooy=y+2.5*fs+5*fss+10;
if (ooy>oy) oy=ooy;
if(mod_or_phi=="Tau(phi)")  oy=oooy;
if(mod_or_phi=="Tau(mod)")  oy=oooy;
// end of scalefactors


//==========================================================
// define output image(s)
newImage("output", "RGB Black", ox, oy, numstacks);
if (background=="White") setForegroundColor(255,255,255);
if (background=="Black") setForegroundColor(0,0,0);
run("Select All");
run("Fill", "stack");
if (foreground=="White") setForegroundColor(255,255,255);
if (foreground=="Yellow") setForegroundColor(255,255,0);
if (foreground=="Black") setForegroundColor(0,0,0);
if (background=="White") setBackgroundColor(255,255,255);
if (background=="Black") setBackgroundColor(0,0,0);
//added extra output stack in version 12
if (extrastack==true) {
	noextra=7;
	if(mod_or_phi=="Tau(phi)")  noextra=3;
	if(mod_or_phi=="Tau(mod)")  noextra=3;
	totalextra=numstacks*noextra;
	newImage("output2", "RGB Black", x, y, totalextra);
}

//==========================================================
// start of iloop (iloop>1 for flim stacks)
for(iloop=0;iloop<numstacks;iloop++){
	for (qq=0;qq<metadata.length;qq++) {
		if(startsWith(metadata[qq],"Time_after_start "+iloop+" =")==true) {
			sinfo=split(metadata[qq],"=");
			time=sinfo[1];
			timepoints[iloop]=parseFloat(time);
		}
	}
}
for (iloop=1;iloop<=numstacks;iloop++) {
i_input=(iloop-1)*3+1;
if (numstacks>1){
	showProgress(iloop, numstacks+1);
}



//==========================================================
// process DC image
orx=10;
ory=10;
selectWindow("input");
setSlice(i_input);
// create threshold images
run("Duplicate...", "title=input-1");
getStatistics(area,mean,min,max);
upthres=max+100;
imm=imax2+1;
setMinAndMax(imin,imm);
if (up_yes==false) setMinAndMax(imin, upthres);
run("8-bit");
setThreshold(1, 254);
run("Threshold", "thresholded remaining black slice");
rename("threshold2");
run("Duplicate...", "title=threshold1");
run("Invert");
run("Divide...", "value=255");
run("16-bit");
getStatistics(npix,tau_div,p,q,taup_sd);
tau_area=npix*tau_div;
selectWindow("threshold2");
run("Divide...", "value=255");
run("Multiply...", "value=2");
run("16-bit");
upthres=imax2;
selectWindow("input");
setSlice(i_input);
run("Duplicate...", "title=input-1");
setMinAndMax(0, upthres);
run("8-bit");
run("Gamma...", "value="+igamma);
run("16-bit");
if (thres_yes==true) {
   	run("Image Calculator...", "image1=input-1 operation=Multiply image2=threshold1");
    	selectWindow("input-1");
}
setMinAndMax(0, 255);
run("8-bit");
selectWindow("input-1");
run("RGB Color");
selectWindow("input-1");
setColor(255,255,255);
if (foreground=="Black") setColor(0,0,0);
drawLine(0,0,x_1,0);drawLine(0,0,0,y_1);drawLine(x_1,0,x_1,y_1);drawLine(0,y_1,x_1,y_1);
run("Copy");
selectWindow("output");
setSlice(iloop);
makeRectangle(orx, ory, x, y);
run("Paste");
//addition in lifetime12
if (extrastack==true) {
	selectWindow("output2");
	ioutput2=(iloop-1)*noextra+1;
	setSlice(ioutput2);
	makeRectangle(0, 0, x, y);
	run("Paste");
}
selectWindow("input-1");
close();

//==========================================================
// process tau_phi image
orx=20+x;
ory=10;
selectWindow("input");
setSlice(i_input+1);
run("Duplicate...", "title=input-1");
run("Image Calculator...", "image1=input-1 operation=Multiply image2=threshold1");
selectWindow("input-1");
nn=1000*taumin;
nm=1000*taumax;
setThreshold(1,20000);
getHistogram(values,histaup,126,nn,nm);
getStatistics(npix,taup_mean,p,q,taup_sd);
x2=npix*(taup_sd*taup_sd+taup_mean*taup_mean);
taup_mean=taup_mean/tau_div;
taup_sd=sqrt((x2-tau_area*taup_mean*taup_mean)/tau_area);
taup_mean=round(taup_mean);
taup_sd=round(taup_sd);
taup_mean=taup_mean/1000; taup_sd=taup_sd/1000;
run("Subtract...","value="+nn);
nn=1000*(taumax-taumin); nn=nn/255;
run("Divide...", "value="+nn);
run("Image Calculator...", "image1=input-1 operation=Multiply image2=threshold1");
selectWindow("input-1");
run("Image Calculator...", "image1=input-1 operation=Add image2=threshold2");
selectWindow("input-1");
setMinAndMax(0, 255);
run("8-bit");
selectWindow("input-1");
run(colortable);
run("RGB Color");
selectWindow("input-1");
setColor(255,255,255);
if (foreground=="Black") setColor(0,0,0);
drawLine(0,0,x_1,0);drawLine(0,0,0,y_1);drawLine(x_1,0,x_1,y_1);drawLine(0,y_1,x_1,y_1);
run("Copy");
selectWindow("output");
makeRectangle(orx, ory, x, y);
run("Paste");
//addition in lifetimes12
if (extrastack==true) {
	selectWindow("output2");
	ioutput2=(iloop-1)*noextra+2;
	setSlice(ioutput2);
	makeRectangle(0, 0, x, y);
	run("Paste");
}
selectWindow("input-1");
close();

//==========================================================
//make tau_phi histogram
max_level=0;
orx=40+3*x;
n2=4;
if (mod_or_phi=="Tau(phi)") {
	orx=30+2*x;
	n2=3;
}
ory=10;
for (i=1; i<=125; i++) {
      if (histaup[i]>max_level) max_level=histaup[i];
}
max_level=max_level/(0.95*y);
newImage("taup_his", "8-bit black", 126, y, 1);
if (background=="White") setForegroundColor(255,255,255);
if (background=="Black") setForegroundColor(0,0,0);
run("Select All");
run("Fill", "stack");
for (i=1; i<=125; i++) {
	ic=i*255/125;
	iin=histaup[i]/max_level+0.05*y;
	setColor(ic,ic,ic);
	drawLine(i,0,i,iin);
}
setColor(255,255,255);
if (foreground=="Black") setColor(0,0,0);
iin=0.05*y;
drawLine(0,iin,125,iin);
run("Size...", "width="+x+" height="+y);
for (i=0; i<=jt;i++) {
	k=i*nt;
	if (k>=taumin) {
		if (k<=taumax) {
			zx=(k-taumin)/(taumax-taumin)*x;
			drawLine(zx,0,zx,iin);
		}
	}
}
setColor(255,255,0);
if (foreground=="Black") setColor(0,0,0);
if (foreground=="White") setColor(255,255,255);
selectWindow("taup_his");
run(colortable);
run("RGB Color");
if (line_yes==true) {
	iin=0.05*y;
	zx=(taup_mean-taumin)/(taumax-taumin)*x;
	drawLine(zx,iin,zx,y);
}
run("Flip Vertically");
selectWindow("taup_his");
run("Copy");
selectWindow("output");
makeRectangle(orx, ory, x, y);
if (mod_or_phi!="Tau(mod)")  {
	run("Paste");
	if (extrastack==true) {
		selectWindow("output2");
		ioutput2=(iloop-1)*noextra+n2;
		setSlice(ioutput2);
		makeRectangle(0, 0, x, y);
		run("Paste");
	}
}
selectWindow("taup_his");
close();
ory=ory+y+fss+1;
selectWindow("output");
setFont("SansSerif" , fss, "bold");
if (mod_or_phi!="Tau(mod)")  {
for (i=0; i<=jtt;i++) {
	k=i*ntt;
	if (k>=taumin) {
		if (k<=taumax) {
			zx=(k-taumin)/(taumax-taumin)*x;
			zx=orx+zx-0.3*fss;
			drawString(k, zx, ory);
			}
		}
	}
	zx=orx+0.5*x-fss*3;ory=ory+fss+1;
	drawString("Tau(phi) ns",zx,ory);
}	


//==========================================================
//process tau_mod image

orx=20+x;
ory=22+fs+y;
n2=5;
n3=7;
if (mod_or_phi=="Tau(mod)") {
	ory=10;
	n2=2;
	n3=3;
}
selectWindow("input");
setSlice(i_input+2);
nn=1000*taumin;
nm=1000*taumax;
run("Duplicate...", "title=input-1");
run("Image Calculator...", "image1=input-1 operation=Multiply image2=threshold1");
selectWindow("input-1");
setThreshold(1,20000);
getHistogram(values,histaum,126,nn,nm);
getStatistics(npix,taum_mean,p,q,taum_sd);
x2=npix*(taum_sd*taum_sd+taum_mean*taum_mean);
taum_mean=taum_mean/tau_div;
taum_sd=sqrt((x2-tau_area*taum_mean*taum_mean)/tau_area);
taum_mean=round(taum_mean);
taum_sd=round(taum_sd);
taum_mean=taum_mean/1000; taum_sd=taum_sd/1000;
nn=1000*taumin;
run("Subtract...","value="+nn);
nn=1000*(taumax-taumin); nn=nn/256;
run("Divide...", "value="+nn);
run("Image Calculator...", "image1=input-1 operation=Multiply image2=threshold1");
selectWindow("input-1");
run("Image Calculator...", "image1=input-1 operation=Add image2=threshold2");
selectWindow("input-1");
setMinAndMax(0, 255);
run("8-bit");
selectWindow("input-1");
run(colortable);
selectWindow("input-1");
setColor(255,255,255);
if (foreground=="Black") setColor(0,0,0);
drawLine(0,0,x_1,0);drawLine(0,0,0,y_1);drawLine(x_1,0,x_1,y_1);drawLine(0,y_1,x_1,y_1);
run("Copy");
selectWindow("output");
makeRectangle(orx, ory, x, y);
if (mod_or_phi!="Tau(phi)")  {
	run("Paste");
	if (extrastack==true) {
		selectWindow("output2");
		ioutput2=(iloop-1)*noextra+n2;
		setSlice(ioutput2);
		makeRectangle(0, 0, x, y);
		run("Paste");
	}
}
selectWindow("input-1");
close();


//==========================================================
//make tau_mod histogram
max_level=0;
orx=40+3*x;
ory=22+fs+y;
if (mod_or_phi=="Tau(mod)") {
	ory=10;
	orx=30+2*x;
	}
for (i=1; i<=125; i++) {
      if (histaum[i]>max_level) max_level=histaum[i];
}
max_level=max_level/(0.95*y);
newImage("taum_his", "8-bit black", 126, y, 1);
if (background=="White") setForegroundColor(255,255,255);
if (background=="Black") setForegroundColor(0,0,0);
run("Select All");
run("Fill", "stack");
for (i=1; i<=125; i++) {
	ic=i*255/128;
	iin=histaum[i]/max_level+0.05*y;
	setColor(ic,ic,ic);
	drawLine(i,0,i,iin);
}
setColor(255,255,255);
if (foreground=="Black") setColor(0,0,0);
iin=0.05*y;
drawLine(0,iin,125,iin);
run("Size...", "width="+x+" height="+y);
for (i=0; i<=jt;i++) {
	k=i*nt;
	if (k>=taumin) {
		if (k<=taumax) {
			zx=(k-taumin)/(taumax-taumin)*x;
			drawLine(zx,0,zx,iin);
		}
	}
}
setColor(255,255,0);
if (foreground=="Black") setColor(0,0,0);
if (foreground=="White") setColor(255,255,255);
if (line_yes==true) {
	iin=0.05*y;
	zx=(taum_mean-taumin)/(taumax-taumin)*x;
	drawLine(zx,iin,zx,y);
}
selectWindow("taum_his");
run(colortable);
run("RGB Color");
if (line_yes==true) {
	iin=0.05*y;
	zx=(taum_mean-taumin)/(taumax-taumin)*x;
	drawLine(zx,iin,zx,y);
}
run("Flip Vertically");
selectWindow("taum_his");
run("Copy");
selectWindow("output");
makeRectangle(orx, ory, x, y);
if (mod_or_phi!="Tau(phi)")  {
	run("Paste");
// added in lifetimes12
	if (extrastack==true) {
		selectWindow("output2");
		ioutput2=(iloop-1)*noextra+n3;
		setSlice(ioutput2);
		makeRectangle(0, 0, x, y);
		run("Paste");
	}
}
selectWindow("taum_his");
close();
ory=ory+y+fss+1;
setFont("SansSerif" , fss, "bold");
selectWindow("output");
if (mod_or_phi!="Tau(phi)")  {
	for (i=0; i<=jtt;i++) {
		k=i*ntt;
		if (k>=taumin) {
			if (k<=taumax) {
				zx=(k-taumin)/(taumax-taumin)*x;
				zx=orx+zx-0.3*fss;
				drawString(k, zx, ory);
			}
		}
	}
	zx=orx+0.5*x-fss*3;ory=ory+fss+1;
	drawString("Tau(mod) ns",zx,ory);
}


//==========================================================
//copy intensity image
selectWindow("input");
setSlice(i_input);
run("Duplicate...", "title=input-2");
setMinAndMax(0, upthres);
run("8-bit");
run("Gamma...", "value="+igamma);


//==========================================================
//add text to output image
selectWindow("output");
setFont("SansSerif" , fs, "bold");
setColor(255,255,0);
if (foreground=="Black") setColor(0,0,0);
if (foreground=="White") setColor(255,255,255);
orx=10;
ory=12+fs+y;
drawString("Intensity", orx, ory);
orx=20+x;
ory=12+fs+y;
if(mod_or_phi=="Tau(mod)") {
	drawString("Tau(mod)", orx, ory);
}else{
	drawString("Tau(phi)", orx, ory);
}
orx=30+2*x;
ory=12+fs+y;
if(mod_or_phi=="Both") drawString("Tau(phi)_Int", orx, ory);
orx=20+x;
ory=24+2*fs+2*y;
drawString("Tau(mod)", orx, ory);
orx=30+2*x;
ory=24+2*fs+2*y;
drawString("Tau(mod)_Int", orx, ory);
setFont("SansSerif" , fss, "bold");
orx=10;
ory=2.5*fs+y;
drawString(fname, orx, ory);
getDateAndTime(year,month,dw,dm,hr,mi,sec,msec);
month=month+1;
ory=ory+fss+2;
if (numstacks==1){
	drawString(dm+"-"+month+"-"+year+"  "+hr+":"+mi+":"+sec,orx,ory);
}else{
	kaka=iloop-1;
	dd=timepoints[kaka]; 
	drawString(dm+"-"+month+"-"+year+"  time="+dd+" "+dtimeunit,orx,ory);
}
ory=ory+fss+2;
drawString("Mean Tau(phi)  ="+taup_mean+"+-"+taup_sd+" ns",orx,ory);
ory=ory+fss+2;
drawString("Mean Tau(mod)="+taum_mean+"+-"+taum_sd+" ns",orx,ory);
ory=ory+fss+2;
drawString("Frequency="+freq+" MHz",orx,ory);
if(mod_or_phi=="Both") {
	//draw intensity scalebar
	rx=0.8*x; rx=round(rx);
	ry=rx/10; ry=round(ry);
	newImage("ramp","8bit ramp", rx, ry, 1);
	run("Gamma...", "value="+igamma);
	setMinAndMax(0, 255);
	run(colortableint);
	run("RGB Color");
	orx=10+0.1*x;orx=round(orx);
	ory=y+2.5*fs+5*fss+10;
	setLineWidth(1);
	setColor(255,255,255);
	if (foreground=="Black") setColor(0,0,0);
	drawLine(0,0,rx-1,0);drawLine(0,0,0,ry-1);drawLine(rx-1,0,rx-1,ry-1);drawLine(0,ry-1,rx-1,ry-1);
	setColor(255,255,0);
	if (foreground=="Black") setColor(0,0,0);
	if (foreground=="White") setColor(255,255,255);
	selectWindow("ramp");
	run("Select All");
	run("Copy");
	selectWindow("output");
	makeRectangle(orx, ory, rx, ry);
	run("Paste");
	run("Select None");
	selectWindow("ramp");
	close();
	selectWindow("output");
	ory=ory+fss+1+ry;
	drawString("0", orx, ory);
	orx=orx+rx-0.375*fs;
	drawString(upthres, orx, ory);
	orx=10+0.5*x-1.20*fs;
	drawString("Intensity", orx, ory);


//==========================================================
//put in 2D-histograms (optional)
	if (his2D_yes==true) {
		setColor(255,255,0);

//copy intensity image
		selectWindow("input");
		setSlice(i_input);
		run("Duplicate...", "title=DC");
		run("Image Calculator...", "image1=DC operation=Multiply image2=threshold1");
		imax3=imax2-1;
		setMinAndMax(imin, imax3);
		run("8-bit");
		nn=1000*taumin;
		nm=1000*taumax;
// process tau_phi image
		selectWindow("input");
		setSlice(i_input+1);
		run("Duplicate...", "title=taup");
		run("Image Calculator...", "image1=taup operation=Multiply image2=threshold1");
		setMinAndMax(nn, nm);
		run("8-bit");
// process tau_mod image
		selectWindow("input");
		setSlice(i_input+2);
		run("Duplicate...", "title=taum");
		run("Image Calculator...", "image1=taum operation=Multiply image2=threshold1");
		setMinAndMax(nn, nm);
		run("8-bit");
//==========================================================
//make tau(phi) intensity 2D histogram
		setColor(255,255,255);
//if (foreground=="Black") setColor(0,0,0);
//if (foreground=="White") setColor(255,255,255);
		run("Image Correlator", "image1=taup image2=DC");
		rename("I_taup");
		drawLine(0,0,0,255);
		drawLine(0,0,255,0);
		drawLine(0,255,255,255);
		drawLine(255,0,255,255);
		resetMinAndMax();
		setColor(0,0,0);
		run("8-bit");
		selectImage("I_taup");
		if (background=="White") {
			run("Select All");
			run("Invert");
		}
		run("Scale...", "x="+hscale+" y="+hscale+" create title=I_taup_1");
		ihy=getHeight();
		ihx=getWidth();
		ihx_1=ihx-1;
		ihy_1=ihy-1;
		setForegroundColor(0,0,0);
		if (foreground=="Black") setForegroundColor(255,255,255);
		selectWindow("I_taup_1");
		for (i=0; i<=ji;i++) {
			k=i*ni;
			if (k>imin) {
				zy=ihy-(k-imin)/(imax2-imin)*ihy;
				zx=ihx/25;
				drawLine(0,zy,zx,zy);
			}
		}
		for (i=0; i<=jt;i++) {
			k=i*nt;
			if (k>taumin) {
				if (k<taumax) {
					zx=(k-taumin)/(taumax-taumin)*ihx;
					zy=ihy-ihy/25;
					drawLine(zx,ihy_1,zx,zy);
				}
			}
		}
		drawLine(0,ihy_1,ihx_1,ihy_1);drawLine(0,0,0,ihy_1);drawLine(0,0,ihx_1,0);drawLine(ihx_1,0,ihx_1,ihy_1);
		run(colortable);
		run("RGB Color");
		makeRectangle(0,0,ihx,ihy);
		run("Copy");
		selectWindow("output");
		orx=50+4*x+fs+fspace;
		ory=10+y-ihy;
		makeRectangle(orx, ory, ihx, ihy);
		run("Paste");
		selectWindow("I_taup_1");
		close();
		selectImage("I_taup");
		close();
		setColor(255,255,0);
		if (foreground=="Black") setColor(0,0,0);
		if (foreground=="White") setColor(255,255,255);
		zy=ory+ihy+fss+1;
		setFont("SansSerif" , fss, "bold");
		selectWindow("output");
		for (i=0; i<=jtt;i++) {
			k=i*ntt;
			if (k>=taumin) {
				if (k<=taumax) {
					zx=(k-taumin)/(taumax-taumin)*ihx;
					zx=orx+zx-0.3*fss;
					drawString(k, zx, zy);
				}
			}
		}
		zx=orx+0.5*ihx-fss*3;zy=zy+fss-1;
		drawString("Tau(phi) ns",zx,zy);
		for (i=0; i<=ji;i++) {
			k=i*ni;
			if (k>=imin) {
				if (k<=imax2) {
					zy=ory+ihy-(k-imin)/(imax2-imin)*ihy+fss/2;
					if (zy<(fss+ory)) zy=fss+ory;
					zx=orx-fspace;
					drawString(k, zx, zy);
				}
			}
		}
		xs=105*fss/20;
		xs=round(xs);xs_1=xs-1;
		newImage("textwindow", "RGB Black", xs, fss_1, 1);
		if (background=="White") setForegroundColor(255,255,255);
		if (background=="Black") setForegroundColor(0,0,0);
		run("Select All");
		run("Fill", "stack");
		setFont("SansSerif" ,fss, "bold");
		setColor(255,255,0);
		if (foreground=="Black") setColor(0,0,0);
		if (foreground=="White") setColor(255,255,255);
		drawString("Intensity",0,fss);
		run("Rotate 90 Degrees Left");
		makeRectangle(0,0,fss_1,xs_1);
		run("Copy");
		close();
		zx=50+4*x+fss; zy=ory+0.5*ihy-0.5*xs;
		selectWindow("output");
		makeRectangle(zx,zy,fss_1,xs_1);
		run("Paste");

//==========================================================
//make tau(mod) intensity 2D histogram
		setColor(255,255,255);
//if (foreground=="Black") setColor(0,0,0);
		run("Image Correlator", "image1=taum image2=DC");
		rename("I_taum");
		drawLine(0,0,0,255);
		drawLine(0,0,255,0);
		drawLine(0,255,255,255);
		drawLine(255,0,255,255);
		resetMinAndMax();
		run("8-bit");
		selectImage("I_taum");
		if (background=="White") {
			run("Select All");
			run("Invert");
		}
		run("Scale...", "x="+hscale+" y="+hscale+" create title=I_taum_1");
		ihy=getHeight();
		ihx=getWidth();
		ihx_1=ihx-1;
		ihy_1=ihy-1;
		selectWindow("I_taum_1");
		setForegroundColor(0,0,0);
		if (foreground=="Black") setForegroundColor(255,255,255);
		for (i=0; i<=ji;i++) {
			k=i*ni;
			if (k>imin) {
				zy=ihy-(k-imin)/(imax2-imin)*ihy;
				zx=ihx/25;
				drawLine(0,zy,zx,zy);
			}
		}
		for (i=0; i<=jt;i++) {
			k=i*nt;
			if (k>taumin) {
				if (k<taumax) {
					zx=(k-taumin)/(taumax-taumin)*ihx;
					zy=ihy-ihy/25;
					drawLine(zx,ihy_1,zx,zy);
				}
			}
		}
		drawLine(0,ihy_1,ihx_1,ihy_1);drawLine(0,0,0,ihy_1);drawLine(0,0,ihx_1,0);drawLine(ihx_1,0,ihx_1,ihy_1);
		run(colortable);
		run("RGB Color");
		makeRectangle(0,0,ihx,ihy);
		run("Copy");
		selectWindow("output");
		orx=50+4*x+fs+fspace;
		ory=2*(12+2*fs+2*y)/3+10;
		makeRectangle(orx, ory, ihx, ihy);
		run("Paste");
		selectWindow("I_taum_1");
		close();
		selectImage("I_taum");
		close();
		zy=ory+ihy+fss+1;
		setColor(255,255,0);
		if (foreground=="Black") setColor(0,0,0);
		if (foreground=="White") setColor(255,255,255);
		setFont("SansSerif" , fss, "bold");
		selectWindow("output");
		for (i=0; i<=jtt;i++) {
			k=i*ntt;
			if (k>=taumin) {
				if (k<=taumax) {
					zx=(k-taumin)/(taumax-taumin)*ihx;
					zx=orx+zx-0.3*fss;
					drawString(k, zx, zy);
				}
			}
		}
		zx=orx+0.5*ihx-fss*3;zy=zy+fss-1;
		drawString("Tau(mod) ns",zx,zy);
		for (i=0; i<=ji;i++) {
			k=i*ni;
			if (k>=imin) {
				if (k<=imax2) {
					zy=ory+ihy-(k-imin)/(imax2-imin)*ihy+fss/2;
					if (zy<(fss+ory)) zy=fss+ory;
					zx=orx-fspace;
					drawString(k, zx, zy);
				}
			}
		}
		xs=105*fss/20;
		xs=round(xs);xs_1=xs-1;
		newImage("textwindow", "RGB Black", xs, fss_1, 1);
		if (background=="White") setForegroundColor(255,255,255);
		if (background=="Black") setForegroundColor(0,0,0);
		run("Select All");
		run("Fill", "stack");
		setFont("SansSerif" ,fss, "bold");
		setColor(255,255,0);
		if (foreground=="Black") setColor(0,0,0);
		if (foreground=="White") setColor(255,255,255);
		drawString("Intensity",0,fss);
		run("Rotate 90 Degrees Left");
		makeRectangle(0,0,fss_1,xs_1);
		run("Copy");
		close();
		zx=50+4*x+fss; zy=ory+0.5*ihy-0.5*xs;
		selectWindow("output");
		makeRectangle(zx,zy,fss_1,xs_1);
		run("Paste");
// make tau(mod) vs tau(phi) 2D histogram
		setColor(255,255,255);
//if (foreground=="Black") setColor(0,0,0);
		run("Image Correlator", "image1=taup image2=taum");
		rename("taup_taum");
		drawLine(0,0,0,255);
		drawLine(0,0,255,0);
		drawLine(0,255,255,255);
		drawLine(255,0,255,255);
		setColor(0,0,0);
		resetMinAndMax();
		drawLine(0,255,255,0);
		run("8-bit");
		selectImage("taup_taum");
		if (background=="White") {
			run("Select All");
			run("Invert");
		}
		run("Scale...", "x="+hscale+" y="+hscale+" create title=I_taump_1");
		ihy=getHeight();
		ihx=getWidth();
		ihx_1=ihx-1;
		ihy_1=ihy-1;
		selectWindow("I_taump_1");
		setForegroundColor(0,0,0);
		if (foreground=="Black") setForegroundColor(255,255,255);
		for (i=0; i<=jt;i++) {
			k=i*nt;
			if (k>taumin) {
				if (k<taumax) {
					zy=ihy-(k-taumin)/(taumax-taumin)*ihy;
					zx=ihx/25;
					drawLine(0,zy,zx,zy);
				}
			}
		}
		for (i=0; i<=jt;i++) {
			k=i*nt;
			if (k>taumin) {
				if (k<taumax) {
					zx=(k-taumin)/(taumax-taumin)*ihx;
					zy=ihy-ihy/25;
					drawLine(zx,ihy_1,zx,zy);
				}
			}
		}
		drawLine(0,ihy_1,ihx_1,ihy_1);drawLine(0,0,0,ihy_1);drawLine(0,0,ihx_1,0);drawLine(ihx_1,0,ihx_1,ihy_1);
		run(colortable);
		run("RGB Color");
		makeRectangle(0,0,ihx,ihy);
		run("Copy");
		selectWindow("output");
		orx=60+4*x+fs+2*fspace+hx;
		ory=2*(12+2*fs+2*y)/3+10;
		makeRectangle(orx, ory, ihx, ihy);
		run("Paste");
		selectWindow("I_taump_1");
		close();
		selectImage("taup_taum");
		close();
		selectWindow("DC");
		close();
		selectWindow("taup");
		close();
		selectWindow("taum");
		close();
		setColor(255,255,0);
		if (foreground=="Black") setColor(0,0,0);
		if (foreground=="White") setColor(255,255,255);
		zy=ory+ihy+fss+1;
		setFont("SansSerif" , fss, "bold");
		selectWindow("output");
		for (i=0; i<=jtt;i++) {
			k=i*ntt;
			if (k>=taumin) {
				if (k<=taumax) {
					zx=(k-taumin)/(taumax-taumin)*ihx;
					zx=orx+zx-0.3*fss;
					drawString(k, zx, zy);
				}
			}
		}
		fspace2=2.2*fss/3;
		if (taumax>=10) fspace2=4.4*fss/3;
		if (ntt==0.5) fspace2=5*fss/3;
		zx=orx+0.5*ihx-fss*3;zy=zy+fss-1;
		drawString("Tau(phi) ns",zx,zy);
		for (i=0; i<=jtt;i++) {
			k=i*ntt;
			if (k>=taumin) {
				if (k<=taumax) {
					zy=ory+ihy-(k-taumin)/(taumax-taumin)*ihy+fss/2;
					if (zy<(fss+ory)) zy=fss+ory;
					zx=orx-fspace2;
					drawString(k, zx, zy);
				}
			}
		}		
		xs=134*fss/20;
		xs=round(xs);xs_1=xs-1;
		newImage("textwindow", "RGB Black", xs, fss_1, 1);
		if (background=="White") setForegroundColor(255,255,255);
		if (background=="Black") setForegroundColor(0,0,0);
		run("Select All");
		run("Fill", "stack");
		setFont("SansSerif" ,fss, "bold");
		setColor(255,255,0);
		if (foreground=="Black") setColor(0,0,0);
		if (foreground=="White") setColor(255,255,255);
		drawString("Tau(mod) ns",0,fss);
		run("Rotate 90 Degrees Left");
		makeRectangle(0,0,fss_1,xs_1);
		run("Copy");
		close();
		zx=50+4*x+2*fss+hx+fspace; zy=ory+0.5*ihy-0.5*xs;
		selectWindow("output");
		makeRectangle(zx,zy,fss_1,xs_1);
		run("Paste");

//==========================================================
// make polar plot
		orx=60+4*x+fs+2*fspace+hx;
		selectWindow("input");
		setSlice(i_input+1);
		run("Duplicate...", "title=tauphi");
		run("Image Calculator...", "image1=tauphi operation=Multiply image2=threshold1");
		selectWindow("tauphi");
		run("32-bit");
		selectWindow("tauphi");
		run("Multiply...", "value="+wfreq);
		run("Duplicate...", "title=tauphi-1");
		run("Image Calculator...", "image1=tauphi operation=Multiply image2=tauphi");
		selectWindow("tauphi");
		run("Add...", "value=1.0");
		selectWindow("input");
		setSlice(i_input+2);
		run("Duplicate...", "title=taumod");
		run("Image Calculator...", "image1=taumod operation=Multiply image2=threshold1");
		selectWindow("taumod");
		run("32-bit");
		selectWindow("taumod");
		run("Multiply...", "value="+wfreq);
		run("Image Calculator...", "image1=taumod operation=Multiply image2=taumod");
		selectWindow("taumod");
		run("Add...", "value=1.0");
		run("Image Calculator...", "image1=tauphi operation=Multiply image2=taumod");
		selectWindow("tauphi");
		run("Square Root");
		selectWindow("tauphi");
		run("Reciprocal");
		run("Image Calculator...", "image1=tauphi-1 operation=Multiply image2=tauphi");
		selectWindow("tauphi-1");
		rename("S");
		run("Multiply...", "value=255");
		setMinAndMax(0,255);
		run("8-bit");
		selectWindow("tauphi");
		rename("G");
		run("Multiply...", "value=255");
		setMinAndMax(0,255);
		run("8-bit");
		run("Image Correlator", "image1=G image2=S");
		rename("taup_taum");
		setColor(255,255,255);
//if (foreground=="Black") setColor(0,0,0);
		drawLine(0,0,0,255);
		drawLine(0,0,255,0);
		drawLine(0,255,255,255);
		drawLine(255,0,255,255);
		getStatistics(area,mean,min,max);
		setMinAndMax(min,max);
		run("8-bit");
		selectImage("taup_taum");
		if (background=="White") {
			run("Select All");
			run("Invert");
		}
		run("Scale...", "x="+hscale+" y="+hscale+" create title=I_taump_1");
		ihy=getHeight();
		ihx=getWidth();
		ihx_1=ihx-1;
		ihy_1=ihy-1;
		selectImage("I_taump_1");
		setForegroundColor(0,0,0);
		if (foreground=="Black") setForegroundColor(255,255,255);
		ory=10+y-ihy;
		drawLine(0,ihy_1,ihx_1,ihy_1);drawLine(0,0,0,ihy_1);drawLine(0,0,ihx_1,0);drawLine(ihx_1,0,ihx_1,ihy_1);
		run(colortable);
		run("RGB Color");
		v=256*256*256-1;
		if (foreground=="Black") v=0;
		for (i=0; i<1000;i++) {
			k=ihx-ihx/2*sin(i*180/1000);
			kk=ihx/2*(cos(i*180/1000)+1);
			setPixel(kk,k,v);
		}	
		makeRectangle(0,0,ihx,ihy);
		run("Copy");
		selectWindow("output");
		makeRectangle(orx, ory, ihx, ihy);
		run("Paste");
		zy=ory+ihy+fss+1;
		zx=orx+0.5*ihx-fss*3;zy=zy+fss-1;
		setFont("SansSerif" , fss, "bold");
		setColor(255,255,0);
		if (foreground=="Black") setColor(0,0,0);
		if (foreground=="White") setColor(255,255,255);
		drawString("Polar plot",zx,zy);
		selectWindow("G");
		close();
		selectWindow("S");
		close();
		selectWindow("taup_taum");
		close();
		selectWindow("I_taump_1");
		close();
		selectWindow("taumod");
		close();
	}
// end of optional display of 2D histograms


//==========================================================
// process tau_phi_int image
	orx=30+2*x;
	ory=10;
	selectWindow("input-2");
	run("16-bit");
	if (thres_yes==true) {
 	 	run("Image Calculator...", "image1=input-2 operation=Multiply image2=threshold1");
 		selectWindow("input-2");
	}
	setMinAndMax(0, 255);
	run("8-bit");
	selectWindow("input");
	setSlice(i_input+1);
	run("Duplicate...", "title=input-1");
	nn=1000*taumin;
	run("Subtract...","value="+nn);
	nn=1000*(taumax-taumin); nn=nn/255;
	run("Divide...", "value="+nn);
	setMinAndMax(0, 255);
	run("8-bit");
	selectWindow("input-1");
	run(colortable);
	run("RGB Color");
	selectWindow("input-1");
	run("RGB Split");
	imageCalculator("Multiply create 32-bit", "input-1 (blue)","input-2");
	setMinAndMax(0, 65535);
	run("8-bit");
	rename("input-3 (blue)");
	selectWindow("input-1 (blue)");
	close();
	imageCalculator("Multiply create 32-bit", "input-1 (green)","input-2");
	setMinAndMax(0, 65535);
	run("8-bit");
	rename("input-3 (green)");
	selectWindow("input-1 (green)");
	close();
	imageCalculator("Multiply create 32-bit", "input-1 (red)","input-2");
	setMinAndMax(0, 65535);
	run("8-bit");
	rename("input-3 (red)");
	selectWindow("input-1 (red)");
	close();
	run("RGB Merge...", "red=[input-3 (red)] green=[input-3 (green)] blue=[input-3 (blue)] gray=*None*");
	rename("input-1");
	selectWindow("input-1");
	setColor(255,255,255);
	if (foreground=="Black") setColor(0,0,0);
	drawLine(0,0,x_1,0);drawLine(0,0,0,y_1);drawLine(x_1,0,x_1,y_1);drawLine(0,y_1,x_1,y_1);
	run("Copy");
	selectWindow("output");
	makeRectangle(orx, ory, x, y);
	run("Paste");
	if (extrastack==true) {
		selectWindow("output2");
		ioutput2=(iloop-1)*noextra+3;
		setSlice(ioutput2);
		makeRectangle(0, 0, x, y);
		run("Paste");
	}
	selectWindow("input-1");
	close();

//==========================================================
// process tau_mod_int image
	orx=30+2*x;
	ory=22+fs+y;
	selectWindow("input");
	setSlice(i_input+2);
	run("Duplicate...", "title=input-1");
	nn=1000*taumin;
	run("Subtract...","value="+nn);
	nn=1000*(taumax-taumin); nn=nn/255;
	run("Divide...", "value="+nn);
	setMinAndMax(0, 255);
	run("8-bit");
	selectWindow("input-1");
	run(colortable);
	run("RGB Color");
	selectWindow("input-1");
	run("RGB Split");
	imageCalculator("Multiply create 32-bit", "input-1 (blue)","input-2");
	setMinAndMax(0, 65535);
	run("8-bit");
	rename("input-3 (blue)");
	selectWindow("input-1 (blue)");
	close();
	imageCalculator("Multiply create 32-bit", "input-1 (green)","input-2");
	setMinAndMax(0, 65535);
	run("8-bit");
	rename("input-3 (green)");
	selectWindow("input-1 (green)");
	close();
	imageCalculator("Multiply create 32-bit", "input-1 (red)","input-2");
	setMinAndMax(0, 65535);
	run("8-bit");
	rename("input-3 (red)");
	selectWindow("input-1 (red)");
	close();
	run("RGB Merge...", "red=[input-3 (red)] green=[input-3 (green)] blue=[input-3 (blue)] gray=*None*");
	rename("input-1");
	selectWindow("input-1");
	setColor(255,255,255);
	if (foreground=="Black") setColor(0,0,0);
	drawLine(0,0,x_1,0);drawLine(0,0,0,y_1);drawLine(x_1,0,x_1,y_1);drawLine(0,y_1,x_1,y_1);
	run("Copy");
	selectWindow("output");
	makeRectangle(orx, ory, x, y);
	run("Paste");
	if (extrastack==true) {
		selectWindow("output2");
		ioutput2=(iloop-1)*noextra+6;
		setSlice(ioutput2);
		makeRectangle(0, 0, x, y);
		run("Paste");
	}
	selectWindow("input-1");
	close();

//==========================================================
//draw int_tau scalebar
	orx=10+0.1*x;orx=round(orx);
	ory=2.5*fs+10+y+7*fss+0.08*x;ory=round(ory);
	newImage("input-1", "8-bit Ramp", rx, rx, 1);
	run("Duplicate...", "title=input-4");
	run("Gamma...", "value="+igamma);
	selectWindow("input-1");
	run("Flip Horizontally");
	selectWindow("input-1");
	run(colortable);
	run("RGB Color");
	selectWindow("input-1");
	run("Rotate 90 Degrees Left");
	run("RGB Split");
	imageCalculator("Multiply create 32-bit", "input-1 (blue)","input-4");
	setMinAndMax(0, 65535);
	run("8-bit");
	rename("input-3 (blue)");
	selectWindow("input-1 (blue)");
	close();
	imageCalculator("Multiply create 32-bit", "input-1 (green)","input-4");
	setMinAndMax(0, 65535);
	run("8-bit");
	rename("input-3 (green)");
	selectWindow("input-1 (green)");
	close();
	imageCalculator("Multiply create 32-bit", "input-1 (red)","input-4");
	setMinAndMax(0, 65535);
	run("8-bit");
	rename("input-3 (red)");
	selectWindow("input-1 (red)");
	close();
	run("RGB Merge...", "red=[input-3 (red)] green=[input-3 (green)] blue=[input-3 (blue)] gray=*None*");
	rename("input-1");
	selectWindow("input-4");
	close();
	selectWindow("input-1");
	run("Rotate 90 Degrees Left");
	setColor(255,255,255);
	if (foreground=="Black") setColor(0,0,0);
	run("Scale...", "x=1 y=0.5 interpolate create title=I_T");
	ry=getHeight();
	rx=getWidth();
	selectWindow("input-1");
	close();
	selectWindow("I_T");
	setLineWidth(1);
	setColor(255,255,255);
	if (foreground=="Black") setColor(0,0,0);
	drawLine(0,0,rx-1,0);drawLine(0,0,0,ry-1);drawLine(rx-1,0,rx-1,ry-1);drawLine(0,ry-1,rx-1,ry-1);
	setColor(255,255,0);
	if (foreground=="Black") setColor(0,0,0);
	if (foreground=="White") setColor(255,255,255);
	run("Select All");
	run("Copy");
	selectWindow("output");
	makeRectangle(orx, ory, rx, ry);
	run("Paste");
	ory=ory+fss+1+ry;
	drawString(taumin, orx, ory);
	orx=orx+rx-0.375*fs;
	drawString(taumax, orx, ory);
	orx=10+0.5*x-1.05*fs;
	drawString("Tau (ns)", orx, ory);
	ory=ory-1-fss-ry;
	orx=10+0.1*x;
	orx=orx-0.75*fs;
	if (orx<0) orx=0;
	drawString(upthres, orx, ory);
	ory=ory+0.5*ry;
	drawString("Int",orx,ory);
	ory=ory+0.5*ry;
	drawString("  0",orx,ory);
	selectWindow("I_T");
	close();
}
//end of tau(phi) or tau(mod) option




//==========================================================
//close temporary images
selectWindow("input-2");
close();
selectWindow("threshold1");
close();
selectWindow("threshold2");
close();

}
// end iloop (not indented with one tab)

//==========================================================
//Display and exit Batch mode, save results
if (extrastack==true) {
	selectWindow("output2");
	setSlice(1);
}
selectWindow("output");
//close();
setSlice(1);
selectWindow("input");
close();
run("Select None");

selectWindow("output");
//Update metadata with display settings==================================

metadata2=metadata2+"--------------Lifetimes15 display settings for generating output-----------\n";
metadata2=metadata2+"Minimal Tau for display (ns)  :"+taumin+" \n";
metadata2=metadata2+"Maximal Tau for display (ns) :"+taumax+" \n";
metadata2=metadata2+"Low threshold_for_Intensity  :"+imin2+" \n";
metadata2=metadata2+"High threshold_for_Intensity :"+imax2+" \n";
metadata2=metadata2+"Gamma_for_Intensity            :"+igamma+" \n";
metadata2=metadata2+"Threshold intensity image:"+thres_yes+" \n";
metadata2=metadata2+"Exclude analysis for I>Ihigh:"+up_yes+" \n";
metadata2=metadata2+"Automatic determination of Ihigh:"+ihigh_yes+" \n";
metadata2=metadata2+"Include 2D-histograms:"+his2D_yes+" \n";
metadata2=metadata2+"Display average lifetime as line in histogram:"+line_yes+" \n";
metadata2=metadata2+"Make additional output image stack:"+extrastack+" \n";
metadata2=metadata2+"Display only Tau(phi), only Tau(mod), or both:"+mod_or_phi+" \n";
metadata2=metadata2+"Colortable for Lifetime Images:"+colortable+" \n";
metadata2=metadata2+"Colortable for Intensity Images:"+colortableint+" \n";
metadata2=metadata2+"Background:"+background+" \n";
metadata2=metadata2+"Foreground :"+foreground+" \n";
metadata2=metadata2+"--------------------------------------------------------------------\n";

setMetadata("Info",metadata2);

saveAs("Tiff", fileout);


if (files>1) close();
if (extrastack==true) {
	selectWindow("output2");
	setMetadata("Info",metadata2);
	saveAs("Tiff", fileout2);
	if (files>1) close();

}

}
setBatchMode("exit and display");



 

