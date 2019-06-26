function spotList = makeSpotList(wins)
%makeSpotList formats an array of window boundaries for EWPC_mapPQ
% inputs:
%   wins -- boundaries of EWPC-space windows selecting peaks to map, in
%           NX4 array for N windows, with each row as
%           [x1min,x1max,x2min,x2max].
% outputs:
%       spotList -- a struct array identifying the spots to be mapped.
%                   Fields are: 
%                   'id' -- a name identifying the spot. (defaults as
%                           numbers)
%                   'spotRangeX1' -- index range around the spot in X1.
%                   'spotRangeX2' -- index range around the spot in X2.
%
%This function is part of the PC-STEM Package by Elliot Padgett in the 
%Muller Group at Cornell University.  Last updated June 26, 2019.

numSpots = size(wins,1);

%make lists of ranges for each spot
id = {};
spotRangeQ1 = {};
spotRangeQ2 = {};

for i=1:numSpots
    id{i} = num2str(i);
    spotRangeQ1{i} = wins(i,1):wins(i,2);
    spotRangeQ2{i} = wins(i,3):wins(i,4);
end

% format struct array
spotList = struct('id',id,'spotRangeQ1',spotRangeQ1,'spotRangeQ2',spotRangeQ2);


end
