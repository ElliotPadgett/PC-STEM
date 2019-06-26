function [ im ] = mapVirtualDetector( data4d, detectormask, doPlot )
%mapVirtualDetector makes an image from a 4D cbed with an arbitrary mask
%   input:
%       cbed4d -- 4D STEM dataset with dimensions ordered [k1,k2,x1,x2].
%       detectormask --  diffraction space detector mask (matches [k1,k2]),
%                 can be binary or positive/negative real values.
%       doPlot -- True or false value indicating if function should plot 
%       mask and resulting image (optional, default false)
%   output: 
%       im -- resultant 2D image.
%
%This function is part of the PC-STEM Package by Elliot Padgett in the 
%Muller Group at Cornell University.  Last updated June 25, 2019.

if nargin<3
    doPlot = false;
end

[N_k1,N_k2,N_x1,N_x2]=size(data4d);

%arrange in 1D 
cbed_lin=reshape(data4d,N_k1,N_k2,N_x1*N_x2); 
im_lin=zeros(N_x1*N_x2,1);

%if a stack of masks was given, iterate through all
numMasks = size(detectormask,3);
im = zeros(N_x1,N_x2,numMasks);
for i=1:numMasks
    
    %calculate for all diffraction patterns
    for j=1:(N_x1*N_x2)
        im_lin(j)=sum(sum(detectormask(:,:,i).*cbed_lin(:,:,j)));
    end
    
    %reshape to 2D image
    im(:,:,i) = reshape(im_lin,N_x1,N_x2);
end

if doPlot
    %plot mask and image
    
    meanCBED = mean(mean(data4d,4),3);
    figure,
    
    for i = 1:numMasks
        subplot(numMasks,3,3*(i-1) + 1),
            plotIM(detectormask(:,:,i)), title('detector mask')
        subplot(numMasks,3,3*(i-1) + 2),
            plotIM(log(meanCBED)), title('mean pattern')
            hold on; visboundaries(detectormask(:,:,i)>0,'EnhanceVisibility',false)
        subplot(numMasks,3,3*(i-1) + 3),
            plotIM(im(:,:,i)), title('image')
    end
    drawnow
end

end

