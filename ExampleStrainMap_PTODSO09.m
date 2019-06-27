%This script demonstrates use of the PC-STEM package for strain mapping
%using scanning nanobeam electron diffraction.
%Written by Elliot Padgett in the Muller Group at Cornell University.  
%Last updated June 26, 2019.


%% Load Data 

fnm = '../PAD_Matlab_ExampleData/09_PTO_640kX_940mm_10um_1p5mrad_50x_50y_50z_108step_x256_y256.raw';
filetag = '09_PTO'; %name stub for automatically generated files.
%filetag = mfilename;  %mfilename gets the name of the currently running
                       %script

[data4d] = single(readSTEM4D(fnm)); %data ordered [k1,k2,x1,y2]

%k1, x1 always refer to the FIRST x or k dimension, k2, x2 the SECOND
[N_k1,N_k2,N_x1,N_x2] = size(data4d);

%make data non-negative, for math convenience
data4d = data4d-min(data4d(:)); 

%% Make 2D images for various imaging modes

% Measure central beam center  for mask positioning
meanCBED = mean(mean(data4d,4),3);
[ center, ~ ] = measureCentralBeam( meanCBED , 1);

%Make ADF image
adfOutR = 62; adfInR = 50; 
adfmask = makeVirtualDetector('adf',124,center,[adfInR,adfOutR]);
adfim=mapVirtualDetector(data4d,adfmask,true);

%% Select ROI to crop and fit window positions for mapping
interactive = true;
savename = [filetag, '.mat'];

[spotList,data4d_roi,adfim_roi,roi,wins] = ...
    setupEWPCMap(data4d,adfim,savename,interactive);

valid = ones(size(adfim_roi)); %Identify what area of ROI should be mapped

%% Make maps of peak positions

tol = 1e-2;
[ spotMaps ] = ewpc_mapSpots( data4d_roi,spotList,tol );  

% Plot Peak-track results 
plotSpotMaps(spotMaps)

%% Map strain tensor with absolute reference

%Set up known lattice vectors for reference
cal = 0.2420; %in Angstrom/pixel.  Calibration depends on camera length
cubic_a=mean([4.152,3.904])/cal;
ref100 = cubic_a*[-1,0];
ref010 = cubic_a*[0,1];
spotReferences = struct('id',{'100','010'},'point',{ref100,ref010});

%Calculate strain
[ StrainComponents,StrainTensors ] = ...
                        calculateStrainMap( spotMaps, spotReferences, 1 );

%Plot strain map results
plotStrainMap(StrainComponents);

%% Map strain tensor with relative reference
%
%use overall mean for relative reference
refWin_x1 = 1:size(adfim_roi,1);
refWin_x2 = 1:size(adfim_roi,2);

[ spotReferences ] = makeRelativeSpotReferences( spotMaps, refWin_x1,refWin_x2 );

%Calculate strain
[ StrainComponents,StrainTensors ] = ...
                        calculateStrainMap( spotMaps, spotReferences, 0 );

%Plot strain map results
plotStrainMap(StrainComponents);

