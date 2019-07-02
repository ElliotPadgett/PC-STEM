function [ roi, wins ] = windowPickingHelper( data4d , im2d, p_roi, p_wins, interactive)
%windowPickingHelper is a GUI to pick a real-space ROI and EWPC-space spot
%   windows for 4D STEM EWPC lattice mapping.
% Inputs:
%   data4d -- 4D STEM dataset to be inspected for EWPC lattice mapping.
%             Dimensions ordered [k1,k2,x1,x2].
%   im2d -- A 2D real-space image, e.g. ADF, to assist in navigating the
%           data.
%   p_roi -- (optional) preset roi value, enter to update or view
%             selections rather than create from scratch.
%   p_wins -- (optional) preset wins list, enter to update or view
%             selections rather than create from scratch.
%   interactive -- (optional) choose whether to run program
%                  interactively (default). If false, the program will plot
%                   the given presets and not wait for user interaction.
% Outputs:
%   roi -- boundaries of real-space region of interest, in order
%          [x1min,x1max,x2min,x2max].
%   wins -- boundaries of EWPC-space windows selecting peaks to map, in
%           NX4 array for N windows, with each row as
%           [x1min,x1max,x2min,x2max].
%
%This function is part of the PC-STEM Package by Elliot Padgett in the 
%Muller Group at Cornell University.  Last updated June 26, 2019.

%Set default inputs
if nargin>2
    preset = 1;
    p_roi = [p_roi(3:4), p_roi(1:2)] ; %follow index convention
else
    preset = 0;
end
if nargin < 5
    interactive = 1;
end

%Set up figure
drawnow
f=figure();
resizeFig(f,2,2)

%Set up guidata structure for shared data
gdata = struct('d4d',data4d);
gdata.EwpcWindows = [];
gdata.Im2D = im2d;
guidata(f,gdata)

%% Set up full image and ROI selection box
%show adf image
gdata.ImageAx = subplot(2,3,1);
plotIM(im2d)
[N_k1,N_k2,N_x1,N_x2] = size(data4d);
axis([1,N_x2,1,N_x1])
title({'2D image:','Select ROI'})

guidata(gcf,gdata)

%make ADF ROI selection box
pause(0.01) % clears some display timing bug
if preset %use provided preset value for roi box position
    pos = [p_roi(1),p_roi(3),p_roi(2)-p_roi(1),p_roi(4)-p_roi(3)];
    gdata.ImageRoiBox = imrect(gdata.ImageAx, pos, ...
        'PositionConstraintFcn',...
        makeConstrainToRectFcn('imrect', get(gca,'XLim'), get(gca,'YLim')));
else %Ask user to place box position
    gdata.ImageRoiBox = imrect(gdata.ImageAx,...
        'PositionConstraintFcn',...
        makeConstrainToRectFcn('imrect', get(gca,'XLim'), get(gca,'YLim')));
end

addNewPositionCallback(gdata.ImageRoiBox,@AdfRoiCallback);
gdata.ImageRoiBox.Deletable = false;


%Display cropped area in separate subfigure
gdata.CroppedImageAx = subplot(2,3,4);
gdata.CroppedImage = plotIM(im2d);
title({'ROI: Select','View Point'})
guidata(gcf,gdata)

pos = round(gdata.ImageRoiBox.getPosition);
pos(pos<1)=1;
setPosition(gdata.ImageRoiBox,pos)
gdata.adfRoi = [pos(1),pos(1)+pos(3),pos(2),pos(2)+pos(4)];
axis(gdata.adfRoi)
cmin = min(min(im2d(gdata.adfRoi(3):gdata.adfRoi(4),gdata.adfRoi(1):gdata.adfRoi(2))));
cmax = max(max(im2d(gdata.adfRoi(3):gdata.adfRoi(4),gdata.adfRoi(1):gdata.adfRoi(2))));
set(gca,'Clim',[cmin,cmax])

guidata(gcf,gdata)

%% Set up ROI-zoomed image and diffraction selection point
%wait until user clicks on ROI to show
gdata.CroppedImage.ButtonDownFcn = @RoiWindowButtonDownFcn;
if preset
    CreateDiffPoint( [pos(1)+pos(3)/2,pos(2)+pos(4)/2] )
else
    uiwait
end
gdata = guidata(gcf);
pause(0.01)


%% Set up EWPC selection boxes

%Declare function to create new boxes with click
gdata.EwpcIm.ButtonDownFcn = @EWPCWindowButtonDownFcn;

if preset
    subplot(2,3,5)
    for i=1:size(p_wins,1)
        w = p_wins(i,:);
        p = [w(3),w(1),w(4)-w(3),w(2)-w(1)];
        
        
        h = imrect(gca,p);
        setPositionConstraintFcn(h,makeConstrainToRectFcn('imrect', get(gca,'XLim'), get(gca,'YLim')))
        addNewPositionCallback(h,@UpdateEWPCWins);
        addlistener( h , 'ObjectBeingDestroyed', @EWPCWinDeleted );
        
        gdata.EwpcWindows = [gdata.EwpcWindows, h];
        guidata(gcf,gdata)
        
        
    end
    UpdateEWPCWins( h.getPosition )
else
    uiwait
end
gdata = guidata(gcf);

%% finish program when user hits space
if interactive % non-interactive skips to end
    gdata = guidata(gcf);
    roi = gdata.adfRoi;

    
    fprintf('When finished, hit any key to return.\n')
    
    waitforbuttonpress;
    while isempty(f.CurrentCharacter)
        waitforbuttonpress;
    end
end

%User has ended function, format return values
gdata = guidata(gcf);
roi = gdata.adfRoi;
roi = [roi(3:4), roi(1:2)] ; %follow index convention
wins = zeros(length(gdata.EwpcWindows),4);
for i=1:length(gdata.EwpcWindows)
    p = round(getPosition(gdata.EwpcWindows(i)));
    p(p<1)=1;
    wins(i,:) = [p(2),p(2)+p(4),p(1),p(1)+p(3)]; %follow index convention
end

pause(0.01)
drawnow
pause(0.01)
fprintf('Done.\n')
end



%% Helper functions and callback functions
function AdfRoiCallback( pos )
gdata = guidata(gcf);

pos = round(pos);
pos(pos<1)=1;
setPosition(gdata.ImageRoiBox,pos)
gdata.adfRoi = [pos(1),pos(1)+pos(3),pos(2),pos(2)+pos(4)];

subplot(gdata.CroppedImageAx)
axis(gdata.adfRoi)
im2d = gdata.Im2D;
cmin = min(min(im2d(gdata.adfRoi(3):gdata.adfRoi(4),gdata.adfRoi(1):gdata.adfRoi(2))));
cmax = max(max(im2d(gdata.adfRoi(3):gdata.adfRoi(4),gdata.adfRoi(1):gdata.adfRoi(2))));
set(gca,'Clim',[cmin,cmax])

guidata(gcf,gdata)

end

function RoiWindowButtonDownFcn(hObject, eventdata, handles)

gdata = guidata(gcf);
pt = get(gdata.CroppedImageAx,'CurrentPoint');

if isfield(gdata,'DiffPoint')
    %move point to mouse click location
    setPosition(gdata.DiffPoint,pt(1,1:2));
    guidata(gcf,gdata)
else
    %Create CBED/EWPC selection point
    CreateDiffPoint( pt(1,1:2) )
    
    uiresume
end


end

function DiffPointCallback( pos )
gdata = guidata(gcf);

diffPointPos = round(pos);
setPosition(gdata.DiffPoint,diffPointPos);

gdata.CbedIm.CData=log(gdata.d4d(:,:,diffPointPos(2),diffPointPos(1)));
gdata.EwpcIm.CData = bsat(ewpc(gdata.d4d(:,:,diffPointPos(2),diffPointPos(1))));

guidata(gcf,gdata)
UpdateEWPCWins( pos )
end

function EWPCWindowButtonDownFcn(hObject, eventdata, handles)
gdata = guidata(gcf);

uiresume

subplot(2,3,5)
h = imrect(gca);
setPositionConstraintFcn(h,makeConstrainToRectFcn('imrect', get(gca,'XLim'), get(gca,'YLim')))
addNewPositionCallback(h,@UpdateEWPCWins);
addlistener( h , 'ObjectBeingDestroyed', @EWPCWinDeleted );

gdata.EwpcWindows = [gdata.EwpcWindows, h];
guidata(gcf,gdata)

UpdateEWPCWins( h.getPosition )

end

function UpdateEWPCWins( pos )

gdata = guidata(gcf);
diffPointPos = round(getPosition(gdata.DiffPoint));

ewpcPointBoxes = gdata.EwpcWindows;

colors = get(gca,'colororder');

for i = 1:length(ewpcPointBoxes)
    h = ewpcPointBoxes(i);
    
    p = round(getPosition(h));
    p(p<1)=1;
    
    a=subplot(length(ewpcPointBoxes),3,3*i);
    PCim = bsat(ewpc(gdata.d4d(:,:,diffPointPos(2),diffPointPos(1))));
    plotIM(PCim);
    ewpcRoi = [p(1),p(1)+p(3),p(2),p(2)+p(4)];
    axis(ewpcRoi)
    imax = max(max(PCim(ewpcRoi(1):ewpcRoi(2),ewpcRoi(3):ewpcRoi(4))));
    imin = min(min(PCim(ewpcRoi(1):ewpcRoi(2),ewpcRoi(3):ewpcRoi(4))));
    caxis([imin,imax]);
    title(sprintf('EWPC spot %d',i))
    coloridx = mod(i-1,length(colors))+1;
    a.XColor = colors(coloridx,:); a.YColor = colors(coloridx,:); a.LineWidth = 4;
    setColor(h,colors(coloridx,:));
end

end

function EWPCWinDeleted( src,evt )
f=gcf;
if strcmp(f.BeingDeleted,'off')
    gdata = guidata(gcf);
    
    gdata.EwpcWindows = gdata.EwpcWindows(gdata.EwpcWindows~=src);
    
    guidata(gcf,gdata)
    UpdateEWPCWins( 1 )
end
end

function CreateDiffPoint( pt )

gdata = guidata(gcf);

%Create CBED/EWPC selection point
gdata.DiffPoint = impoint(gca,pt);
guidata(gcf,gdata)
addNewPositionCallback(gdata.DiffPoint,@DiffPointCallback);
setPositionConstraintFcn(gdata.DiffPoint,makeConstrainToRectFcn('impoint', get(gdata.ImageAx,'XLim'), get(gdata.ImageAx,'YLim')))
diffPointPos = round(getPosition(gdata.DiffPoint));
guidata(gcf,gdata)
setPosition(gdata.DiffPoint,diffPointPos);

%Show CBED at selected point
gdata.CbedAx = subplot(2,3,2);
gdata.CbedIm = plotIM(log(gdata.d4d(:,:,diffPointPos(2),diffPointPos(1))));
title('log(CBED)')

%Show EWPC at selected point
gdata.EwpcAx=subplot(2,3,5);
gdata.EwpcIm=plotIM(bsat(ewpc(gdata.d4d(:,:,diffPointPos(2),diffPointPos(1)))));
title({'EWPC: Pick', 'fit windows'})

guidata(gcf,gdata)

end
