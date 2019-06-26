function [spotList,data4d_roi,im2d_roi,roi,wins] = setupEWPCMap(data4d,im2d,savename,interactive)
%setupEWPCMap runs the window-picking GUI and saves/loads user-set windows
% Inputs:
%   data4d -- 4D STEM dataset to be inspected for EWPC lattice mapping.
%             Dimensions ordered [k1,k2,x1,x2].
%   im2d -- A 2D real-space image, e.g. ADF, to assist in navigating the
%           data.
%   savename -- .mat filename for saving window and ROI settings, or to
%               check for pre-existing settings
%   interactive -- (optional) choose whether to run program
%                  interactively (default). If false, the program will plot
%                  the windows loaded from the save file (if available) 
%                  and not wait for user interaction.
% Outputs:
%   spotList -- a struct array identifying the spots to be mapped.
%                   Fields are: 
%                   'id' -- a name identifying the spot. (defaults as
%                           numbers)
%                   'spotRangeX1' -- index range around the spot in X1.
%                   'spotRangeX2' -- index range around the spot in X2.
%   data4d_roi -- 4D STEM dataset cropped to the real-space ROI
%   im2d_roi -- 2D real-space image cropped to the real-space ROI
%   roi -- boundaries of real-space region of interest, in order
%          [x1min,x1max,x2min,x2max].
%   wins -- boundaries of EWPC-space windows selecting peaks to map, in
%           NX4 array for N windows, with each row as
%           [x1min,x1max,x2min,x2max].
%
%This function is part of the PC-STEM Package by Elliot Padgett in the 
%Muller Group at Cornell University.  Last updated June 26, 2019.

if interactive
 %get rectangular crop area from user with helper GUI
    if isfile(savename) %check if save file already exists
        fprintf('Loading saved ROI and windows from: %s\n',savename)
        load(savename);
        [roi,wins] = windowPickingHelper(data4d,im2d,roi,wins);
    else
        [roi,wins] = windowPickingHelper(data4d,im2d);
    end
    ROI_x1 = roi(1):roi(2); ROI_x2 = roi(3):roi(4);
    
    % Make list of spots formatted for mapping function
    spotList = makeSpotList(wins);
    fprintf('Saving ROI and windows to: %s\n',savename)
    save(savename,'roi','wins','spotList','ROI_x1','ROI_x2');
else
    %load saved fit windows
    fprintf('Loading saved ROI and windows from: %s\n',savename)
    load(savename);
    windowPickingHelper(data4d,im2d,roi,wins,false);
end

data4d_roi=data4d(:,:,ROI_x1,ROI_x2);
im2d_roi = im2d(ROI_x1,ROI_x2);

end

