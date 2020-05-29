# Quantitative-Calcium-Imaging
Custom code used for the analysis of quantitative calcium imaging is stored here.
These scripts were used for the publication on the development of Tq-Ca-FLITS.

### FLIM ImageJ Package
This folder contains all ImageJ scripts created by Dorus Gadella, to convert lifetime stacks (obtained with one of our FLIM setups) to lifetime images (Convert_flim-I-to-ij_v2.ijm and Convert_flim-II-to-ij.ijm) and to convert these lifetime images into a multipanel display (Lifetimes14.ijm and Lifetimes15.ijm). All other files in this folder are required for these scripts to function properly.

### R scripts for modelling
The scripts in this folder created by FH van der Linden, were used to determine the calcium sensitivity of Tq-Ca-FLITS and its pH sensitivity.
The CSV files in the folder contain the raw measurements needed to repeat the analysis.

### ImageJ masking and concentration
The scripts in this folder, created by FH van der Linden, were used to mask organoid and HUVEC lifetime stacks (a timelapse and single image variant) and to calculate the lifetime data into a calcium concentration based in the in vivo calibration. The required input lifetime data is in the format of G and S coordinates. See the scripts for the specific required input.
