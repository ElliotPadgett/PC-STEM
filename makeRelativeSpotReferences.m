function [ spotReferences ] = makeRelativeSpotReferences( spotMaps, refWin_x1,refWin_x2 )
%MakeRelativeSpotReferences Makes reference list for StrainTensorMap
%   inputs:
%       spotMaps -- struct array containing maps of maximum spot index Q1,Q2
%                 for each spot in spotList. Fields are:
%                 'id' -- a name identifying the spot.
%                 'Q1map' -- map of Q1 index value at peak maximum.
%                 'Q2map' -- map of Q2 index value at peak maximum.
%                 'x1map' -- map of x1 vector component for the spot.
%                 'x2map' -- map of x2 vector component for the spot.
%                 'spotlength' -- map of the spot vector length.
%                 'spotangle' -- map of the spot vector angle in degrees.
%       refWin_x1,refWin_x2 -- index ranges identifying area to be averaged 
%                              to make reference spots. If only one
%                              argument is supplied, it is assumed to be a
%                              2D mask identifying the averaging area.
%   outputs:
%       spotReferences -- struct array containing reference points
%                         corresponding to points in Qmaps. Fields are:
%                         'id' -- a name identifying the spot.
%                         'point' -- reference spot location [q1,q2]
%
%This function is part of the PC-STEM Package by Elliot Padgett in the 
%Muller Group at Cornell University.  Last updated June 26, 2019.


numspots=length(spotMaps);
ids = {};
refs = {};

for i=1:numspots
    ids{i} = spotMaps(i).id;
    
    refs{i} = [mean(mean(spotMaps(i).VectorX1(refWin_x1))),...
        mean(mean(spotMaps(i).VectorX2(refWin_x1)))];
end

spotReferences = struct('id',ids,'point',refs);

end

