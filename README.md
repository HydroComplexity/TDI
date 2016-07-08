# Topographic Depression Identification

### Install
Install Matlab and mex with the following toolboxes:
- Image Processing, 
- Image Accquisition,
- Mapping.

### Compile
Compile C functions to mex executable functions
- All functions are store in VER15/FUNCTIONS folder
- in Matlab, go the folder VER15, type mexfile

### Run
The main code controlling all the processes is Main_BPNM.m
- This file is setup for an example in Bird Point New Madrid area. The lidar data is stored in "LiDAR_Data/BPNM/BPNM_DEM2.tif"
- Modify this info in section "LOAD TOPOGRAPHIC DATA" to change the input lidar data for your study area.
- After that, click Run arrow or type Main_BPNM to run the code.
- The results are stored in RESULTS folder.
- You need to write your own script to extract information from result file, depending on the purposes you want.    

### Note
Depending on the size and the number of large depressions on the landscape, the complete time can be very short or long due the loop search algorithm.





