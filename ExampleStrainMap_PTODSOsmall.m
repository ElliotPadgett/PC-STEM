%This script demonstrates use of the PC-STEM package for strain mapping
%using scanning nanobeam electron diffraction.
%Written by Elliot Padgett in the Muller Group at Cornell University.  
%Last updated July 18, 2019.


%% Load Data 

fnm = '../PAD_Matlab_ExampleData/PTODSO_small_128x130x64x64_float_x64_y64.raw';
filetag = 'PTODSO_small'; %name stub for automatically generated files.

[data4d] = single(readSTEM4D(fnm)); %data ordered [k1,k2,x1,y2]

%k1, x1 always refer to the FIRST x or k dimension, k2, x2 the SECOND
[N_k1,N_k2,N_x1,N_x2] = size(data4d);

%make data non-negative, for math convenience
data4d = data4d-min(data4d(:)); 

%% Browse data
browseSTEM4D(data4d)

%% Make 2D images for various imaging modes

%Make ADF image
center = [63,63];
adfOutR = 62; adfInR = 50; 
adfmask = makeVirtualDetector('adf',124,center,[adfInR,adfOutR]);
adfim=mapVirtualDetector(data4d,adfmask,true);

%Make DF image
dfmask = makeVirtualDetector('bf',124,[74,80],3);
dfim=mapVirtualDetector(data4d,dfmask,true);


%% Set up map: select ROI to crop and fit window positions
interactive = true; %If true, this will wait for you to choose or edit the 
                    %roi and fit windows and then exit the GUI.  If false, 
                    %it will load saved values and proceed without waiting.
savename = [filetag, '.mat']; %your roi and window choices will be saved 
                              %in this file

%run the GUI to set up a strain map                              
[spotList,data4d_roi,dfim_roi,roi,wins] = ...
    setupEWPCMap(data4d,dfim,savename,interactive);
%spotList is a struct that identifies the spots you want to map
%data4d_roi and dfim_roi are cropped to the ROI specified

valid = ones(size(dfim_roi)); %Identify what area of ROI should be mapped

%% Make maps of peak positions

tol = 1e-2; %Sets tolerance for peak finding.  A smaller value may improve 
            %precision but will take longer.
[ spotMaps ] = ewpc_mapSpots( data4d_roi,spotList,tol );  
%spotMaps contains the results of the spot fitting.  For each spot, it
%includes: Q1map and Q2map, which are the indices of the spot in
%non-zero-centered EWPC space; VectorX1 and VectorX2, which are the vector
%components of the spot in zero-centered EWPC space; VectorLength and
%VectorAngle, which are the length and angle of the vector to the spot in
%zero-centered EWPC space.

% Plot Peak-track results 
plotSpotMaps(spotMaps)

%% Map strain tensor with absolute reference
%This section calculates a strain tensor using user-defined reference
%lattice vectors for a known lattice structure

%Set up known lattice vectors for reference
cal = 0.2420; %in Angstrom/pixel.  Calibration depends on camera length
cubic_a=4.028/cal; % Average lattice parameter in pixels
ref1 = cubic_a*[1,0]; %100 lattice vector
ref2 = cubic_a*[0,1]; %010 lattice vector

%Build a spotReferences struct containing reference vectors for all spots
%in spotMaps
spotReferences = struct('id',{'100','010'},'point',{ref1,ref2});

%Identify strain coordinate system.  If true, the strain directions match
%follow the given reference vectors.  If false, the strain directions match
%the directions of the diffraction pattern.
latticeCoords = true;

%Calculate strain tensor
[ StrainComponents,StrainTensors ] = ...
           calculateStrainMap( spotMaps, spotReferences, latticeCoords );

%Plot strain map results
plotStrainMap(StrainComponents);

%% Map strain tensor with relative reference
%This section calculates a strain tensor using an internal reference
%determined by averaging spot positions in some area

%Define image-space region to use for relative reference. This version
%averages the entire ROI.
refWin_x1 = 1:size(dfim_roi,1);
refWin_x2 = 1:size(dfim_roi,2);

%Generate the spotReferences struct from the average positions
[ spotReferences ] = ...
    makeRelativeSpotReferences( spotMaps, refWin_x1,refWin_x2 );

%Identify strain coordinate system.  
latticeCoords = true;

%Calculate strain tensor
[ StrainComponents,StrainTensors ] = ...
         calculateStrainMap( spotMaps, spotReferences, latticeCoords );

%Plot strain map results
plotStrainMap(StrainComponents);

