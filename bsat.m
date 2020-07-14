function [ im ] = bsat( im, bsize, center)
%bblank saturates image im above the maximum value outside specified region
%   input:
%       im -- 1d-44 image to have intensity corrected
%       bsize -- edge half-length in pixels (default 2) of region
%       center -- center of region (defaults to image center)
%   output:
%       im -- intensity-corrected image
%
%This function is part of the PC-STEM Package by Elliot Padgett in the 
%Muller Group at Cornell University.  Last updated June 26, 2019.

%default values
switch nargin
    case 2
        center = 0.5*(size(im)+1);
    case 1
        center = 0.5*(size(im)+1); bsize = 3;
end

Nd=ndims(im);
mask = true(size(im));

if Nd==2 & (size(im,1)==1)
    mask(round(center(2)-bsize):round(center(2)+bsize))=0;
elseif Nd==2 & (size(im,2)==1)
    mask(round(center(1)-bsize):round(center(1)+bsize))=0;
elseif Nd==2
    mask(round(center(1)-bsize):round(center(1)+bsize),...
        round(center(2)-bsize):round(center(2)+bsize)) = 0;    
elseif Nd==3 
    mask(round(center(1)-bsize):round(center(1)+bsize),...
        round(center(2)-bsize):round(center(2)+bsize),:) = 0;
elseif Nd==4
    mask(round(center(1)-bsize):round(center(1)+bsize),...
        round(center(2)-bsize):round(center(2)+bsize),:,:) = 0;
end

satval = max(abs(im(mask)));
im(im>satval) = satval;
im(im<-1*satval) = -1*satval;

end

