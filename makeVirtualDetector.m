function [ detectormask ] = makeVirtualDetector(detectortype, pixels, center, detectorprops)
%makeVirtualDetector creates a reciprocal-space virtual detector mask
%   Use with mapVirtualDetector to create a 2D image from a 4D dataset
%   input:
%       detectortype -- type of virtual detector desired.  Options include 
%           'adf', 'bf', 'dpc', 'com'
%       Ksize -- number of pixels in reciprocal space detector
%       center -- center position [x1,x2] of virtual detector
%       detectorprops -- vector describing parameters of mask.  Form
%           depends on detector type:
%               'adf': [innerRadius, outerRadius].  An annular detector.
%               'bf': [radius]. A circular detector.
%               'dpc': []. Differential phase contrast detectors filling 
%                          the full plane. Combine with adf or bf for 
%                          other shapes.
%               'com': []. 'Center of Mass' detectors.
%   output: 
%       detectormask --  diffraction space detector mask (matches [k1,k2]),
%                 may be binary or positive/negative real values. For DPC
%                 and COM detectors, the x1 and x2 masks are returned as a
%                 stack.
%
%This function is part of the PC-STEM Package by Elliot Padgett in the 
%Muller Group at Cornell University.  Last updated June 25, 2019.

switch detectortype
    case 'adf' % Annular detector
        
        adfInR = detectorprops(1);
        adfOutR = detectorprops(2); 
        detectormask = makeDisk(center,adfOutR,1,[pixels,pixels])...
                       - makeDisk(center,adfInR,1,[pixels,pixels]);
                   
    case 'bf' % Circular detector
     
        R = detectorprops;
        detectormask = makeDisk(center,R,1,[pixels,pixels]);
        
    case 'dpc' % DPC detector: positive on one half, negative on the other
        %This is a 'crude' version with a hard edge, even if a non-integer
        %center is given.
        [X2,X1]=meshgrid(1:pixels,1:pixels);
        dpc1 = double(X1 > center(1)) - double(X1 < center(1));
        dpc2 = double(X2 > center(2)) - double(X2 < center(2));
        detectormask = cat(3,dpc1,dpc2);
        
    case 'com' % COM detector: ramp weighting pixels by position
        %This version ignores the given center value. to "center" the COM
        %images, simply subtract the known center position.
        [X2,X1]=meshgrid(1:pixels,1:pixels);
        detectormask = cat(3,X1,X2);
        
    otherwise
        error('Unrecognized detector type.')
end
    
end


function [ diskim ] = makeDisk( center, radius, interpfactor, imsize )
%makeDisk Makes an image with a disk in it for making virtual detectors.
% inputs:
%   center -- [row,col] of disk center in final image.
%   radius -- radius of disk in final image.
%   interpfactor -- factor to scale image by to give disk smooth edges.
%   imsize -- [row,col] size of final image.
% outputs:
%   diskim -- final image containing disk.
%
% Example usage to make a virtual ADF detector:
%   adfmask = MakeDisk(center,outRad,1,[N_k1,N_k2]) - MakeDisk(center,inRad,1,[N_k1,N_k2]);
%
% Elliot Padgett -- Muller Group, Cornell University -- December 2018

if nargin<4
    imsize=[128,128];
end
if nargin<3
    interpfactor=4;
end

%diskim=zeros(interpfactor*imsize);
Nr=interpfactor*imsize(1); Nc=interpfactor*imsize(2);  
w=1/interpfactor;
[r,c]=meshgrid(linspace(0.5*(1+w),imsize(1)+0.5*(1-w),Nr), ...
               linspace(0.5*(1+w),imsize(2)+0.5*(1-w),Nc));

diskim=double((r-center(1)).^2+(c-center(2)).^2<=radius^2);

%Bin back down to imsize
diskim=blockproc(diskim,[interpfactor,interpfactor],@(A) mean2(A.data));

end