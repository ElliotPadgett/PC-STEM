function [ center, alpha ] = measureCentralBeam( cbed , doPlot)
%measureCentralBeam Measures center and convergence angle from central beam
% This is not a highly accurate method, but is typically sufficient for
% placement of virtual detector masks. 
%input:
%   cbed -- 2D diffraction pattern used for measurement. PACBED or average 
%           of entire dataset, is recommended.
%   doPlot -- logical indicating whether fit should be plotted when
%             complete, for verification.  Optional -- default is 0.
%output:
%   center -- beam center position [x, y] in units of pixels
%   alpha -- convergence semiangle in pixels
%
%This function is part of the PC-STEM Package by Elliot Padgett in the 
%Muller Group at Cornell University.  Last updated June 26, 2019.

if nargin < 2
    doPlot=0;
end

%guess threshold to identify central beam
cbed_binary= imbinarize(mat2gray(cbed));

%get flat disk intensity (guess from 90th percentile above threshold)
flat_int=prctile(cbed(cbed_binary),90);
half_int=flat_int/2;

%Estimate convergence angle from 90% radius
beam_bw=cbed>=half_int;
beamArea=sum(sum(beam_bw));
alpha=sqrt(beamArea/pi);

%calculate centroid of binarized beam
stats=regionprops(beam_bw,'centroid');
center=stats.Centroid;

if doPlot
    figure, hold on
    plotIM(cbed)
    PlotCircle(center(1),center(2),alpha,100,'r')
    title('Central Beam Fit')
end
end

function PlotCircle(Center_col,Center_row,R,N,Color)

t=(0:N)*2*pi/N;
x=R*cos(t)+Center_col;
y=R*sin(t)+Center_row;
plot(x,y,Color);

end