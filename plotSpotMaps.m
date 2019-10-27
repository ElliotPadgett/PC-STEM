function [] = plotSpotMaps(spotMaps)
%ewpc_mapSpots Produces a figure showing the fitted EWPC peaks
%   input:
%       spotMaps -- struct array containing maps of maximum spot index Q1,Q2
%                 for each spot in spotList. Fields are:
%                 'id' -- a name identifying the spot.
%                 'Q1map' -- map of Q1 index value at peak maximum.
%                 'Q2map' -- map of Q2 index value at peak maximum.
%                 'x1map' -- map of x1 vector component for the spot.
%                 'x2map' -- map of x2 vector component for the spot.
%                 'spotlength' -- map of the spot vector length.
%                 'spotangle' -- map of the spot vector angle in degrees.
%
%This function is part of the PC-STEM Package by Elliot Padgett in the 
%Muller Group at Cornell University.  Last updated June 26, 2019.

figure
numSpots = length(spotMaps);
for i=1:numSpots
    subplot(numSpots,4,i*4-3)
    plotIM(spotMaps(i).VectorX1,'AlphaData',~isnan(spotMaps(i).VectorX1));
    title(sprintf('spot %d: X1,',i)),colorbar
    
    subplot(numSpots,4,i*4-2)
    plotIM(spotMaps(i).VectorX2,'AlphaData',~isnan(spotMaps(i).VectorX2));
    title(sprintf('X2,')),colorbar
    
    subplot(numSpots,4,i*4-1)
    plotIM(spotMaps(i).VectorLength,'AlphaData',~isnan(spotMaps(i).VectorLength));
    title(sprintf('length,')),colorbar
    
    subplot(numSpots,4,i*4)
    plotIM(spotMaps(i).VectorAngle,'AlphaData',~isnan(spotMaps(i).VectorAngle));
    title(sprintf('angle (^o)')),colorbar
end
colormap parula
end

