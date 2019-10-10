function [ spotMaps ] = ewpc_mapSpots( data4d,spotList,tol,valid)
%ewpc_mapSpots Produces maps of peak fits in CBEDFFT
%   input:
%       data4d -- the  4D dataset, ordered k1,k2,x1,x2
%       spotList -- a struct array identifying the spots to be mapped.
%                   Fields must be:
%                   'id' -- a name identifying the spot.
%                   'spotRangeX1' -- index range around the spot in X1.
%                   'spotRangeX2' -- index range around the spot in X2.
%                   Additional Field may be supplied:
%                   'valid' -- a mask for each spot that identifies the
%                               subregion in which to perform the fitting
%                               (for that spot)
%       tol -- (optional) convergence tolerance. Default is 1e-3.
%       valid -- (optional) a mask to identify a subregion in which to
%                perform the fitting procedure (applies to all spots)
%   output:
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

%initialize data info
[N_k1,N_k2,N_x1,N_x2]=size(data4d);
win=window2(N_k1,N_k2,@hann);
minval=min(data4d(:));

if nargin<3
    tol=1e-3;
end
if nargin<4
    valid = ones(N_x1,N_x2);
end
flag_spotDepValid = isfield(spotList,'valid');


%Set up Q indexing
q1range = 1:N_k1; q2range = 1:N_k2;
[Q2,Q1]=meshgrid(q2range,q1range);

%Prealocate Q maps
spotMaps = struct('id',{spotList.id},'Q1map',zeros(N_x1,N_x2),'Q2map',zeros(N_x1,N_x2));

fprintf('Starting CBEDFFT spot-tracking map.\n'); tic
totalspots = length(spotMaps)*sum(valid(:)); spotsfit = 0; % For tracking progress
if flag_spotDepValid
    totalspots = 0;
    for s=1:length(spotList)
        totalspots = totalspots+sum(spotList(s).valid(:).*valid(:));
    end
end
%iterate through spatial positions
for j=1:N_x1
    for k=1:N_x2
        if valid(j,k)
            
            %select local CBED
            CBED = data4d(:,:,j,k);
            EWPC = ewpc(CBED);
            %define continuous Fourier transform
            PeakFun = @(x) -abs(cft2(win.*log(CBED-minval+0.1),x(1),x(2),1)); %x is [q1,q2]
            
            %iterate through spots of interest
            for s = 1:length(spotList)
                
                if flag_spotDepValid & ~spotList(s).valid(j,k)
                    spotMaps(s).Q1map(j,k) = nan;
                    spotMaps(s).Q2map(j,k) = nan;
                    continue
                else
                                    
                    %Get spot locations from input struct
                    spot_ROI_q1 = spotList(s).spotRangeQ1;
                    spot_ROI_q2 = spotList(s).spotRangeQ2;
                    spotNeighborhood = EWPC(spot_ROI_q1,spot_ROI_q2);
                    
                    %Find rough maximum of peak
                    [~,maxidx] = max(spotNeighborhood(:));
                    Q1_roi = Q1(spot_ROI_q1,spot_ROI_q2);
                    Q2_roi = Q2(spot_ROI_q1,spot_ROI_q2);
                    Q1max = Q1_roi(maxidx); Q2max = Q2_roi(maxidx);
                    
                    %Search for spot peak in continuous Fourier transform
                    constrainedPeakFun = @(x) ConstrainedFun(x,PeakFun,...
                        [spot_ROI_q1(1),spot_ROI_q1(end)],[spot_ROI_q2(1),spot_ROI_q2(end)]);
                    [peakQ] = fminsearch(constrainedPeakFun,[Q1max,Q2max],optimset('TolX',tol,'TolFun',tol));
                    
                    %Assign in maps
                    spotMaps(s).Q1map(j,k) = peakQ(1);
                    spotMaps(s).Q2map(j,k) = peakQ(2);
                    
                    spotsfit = spotsfit+1;
                    
                end
            end
        else
            for s = 1:length(spotList)
                spotMaps(s).Q1map(j,k) = nan;
                spotMaps(s).Q2map(j,k) = nan;
            end
        end
    end
    fprintf('Completed %.1f percent of map. About %d seconds remain.\n',...
        100*spotsfit/totalspots,round(toc/spotsfit*(totalspots-spotsfit)));
end

t=toc;
fprintf('\nDone. Total time: %.1f s. Time per peak fit: %.3f s.\n',t,t/totalspots)

% Add vector-form maps to spotMaps struct
spotMaps = calculateSpotMapVectors( spotMaps, N_k1, N_k2);

end

%%
function y = ConstrainedFun(x,fun,win1,win2)
%adds constraint to objective function fun, which is assumed to be always
%negative, by adding a positive "cone of shame" outside the specified
%window
if x(1) < win1(1) || x(1) > win1(2) || x(2) < win2(1) || x(2) > win2(2)
    cent = [mean(win1),mean(win2)];
    y = norm(x-cent);
else
    y = fun(x);
end
end

%%
function w = window2(N,M,w_func)
%Makes a 2D window function for FFT

wc=window(w_func,N);
wr=window(w_func,M);
[maskr,maskc]=meshgrid(wr,wc);

w=maskr.*maskc;

end

%%
function [ spotMaps_updated ] = calculateSpotMapVectors( spotMaps, N_k1, N_k2)
%calculateSpotMapVectors Calculates vector components, length, and angle
%   with detector with N_k1,  N_k2 pixels, i.e. zero=(N_k1,N_k2)/2+1. 
%   These are added to the spotMaps struct
center_1 = N_k1/2+1;
center_2 = N_k2/2+1;
numSpots = length(spotMaps);
spotMaps_updated=spotMaps;
for i=1:numSpots
    x1map=spotMaps(i).Q1map;
    x1map=x1map-center_1; 
    x1map(x1map>(center_1-1)) = x1map(x1map>(center_1-1))- N_k1;
    
    x2map=spotMaps(i).Q2map;
    x2map=x2map-center_2; 
    x2map(x2map>(center_2-1)) = x2map(x2map>(center_2-1))- N_k2;
    
    spotlength=sqrt(x1map.^2+x2map.^2);
    spotangle=atan2d(x1map,x2map);
    
    spotMaps_updated(i).VectorX1 = x1map;
    spotMaps_updated(i).VectorX2 = x2map;
    spotMaps_updated(i).VectorLength = spotlength;
    spotMaps_updated(i).VectorAngle = spotangle;
    
end

end

