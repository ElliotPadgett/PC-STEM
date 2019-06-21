function [] = browseSTEM4D( data4d )
%browseSTEM4D A browser for viewing 4D STEM real and diffraction space
% Inputs:
%   data4d -- 4D STEM dataset to be inspected for EWPC lattice mapping.
%             Dimensions ordered [k1,k2,x1,x2].
%
%This function is part of the PC-STEM Package by Elliot Padgett in the 
%Muller Group at Cornell University.  Last updated June 21, 2019.


[N_k1,N_k2,N_x1,N_x2] = size(data4d);

%Set up figure
f=figure('color','w');
f.KeyPressFcn = @WindowKeyPressFcn;
resizeFig(f,1.5,1)

%Set up guidata structure for shared data
gdata = struct('d4d',data4d);
gdata.Transf = @(x) x.^.001;
guidata(f,gdata)

%% Set up image and CBED display
%show adf image
StartPointDiff = ceil([N_x1,N_x2]/2);
StartPointReal = ceil([N_k1,N_k2]/2);

gdata.ImageAx = subplot(1,2,1);
gdata.Image = plotIM(squeeze(data4d(StartPointReal(1),StartPointReal(2),:,:)));
title(sprintf('Real Space (%d, %d)',StartPointDiff))
gdata.Image.ButtonDownFcn = @WindowButtonDownFcn;

gdata.CbedAx = subplot(1,2,2);
t = gdata.Transf;
gdata.Cbed = plotIM(t(squeeze(data4d(:,:,StartPointDiff(1),StartPointDiff(2)))));
title(sprintf('Diffraction Space (%d, %d)',StartPointReal))
gdata.Cbed.ButtonDownFcn = @WindowButtonDownFcn;

guidata(gcf,gdata)


gdata.DiffPoint = impoint(gdata.ImageAx,StartPointDiff);
setPositionConstraintFcn(gdata.DiffPoint,makeConstrainToRectFcn('impoint', get(gdata.ImageAx,'XLim')-[0,1], get(gdata.ImageAx,'YLim')-[0,1]))
addNewPositionCallback(gdata.DiffPoint,@DiffPointCallback);
setColor(gdata.DiffPoint,'r')

gdata.RealPoint = impoint(gdata.CbedAx,StartPointReal);
addNewPositionCallback(gdata.RealPoint,@RealPointCallback);
setPositionConstraintFcn(gdata.RealPoint,makeConstrainToRectFcn('impoint', get(gdata.CbedAx,'XLim')-[0,1], get(gdata.CbedAx,'YLim')-[0,1]))
setColor(gdata.RealPoint,'r')

guidata(gcf,gdata)
pause(0.01)

%% Make contrast controls
MakeUiControls()

end

function bselection(source,event)
gdata = guidata(gcf);


Gamma = gdata.s1.Value;
gdata.Transf = @(x) x.^Gamma;

t = gdata.Transf;
PointPos = round(getPosition(gdata.DiffPoint));
gdata.Cbed.CData = t(squeeze(gdata.d4d(:,:,PointPos(2),PointPos(1))));
guidata(gcf,gdata)

end

%% Helper functions and callback functions

function WindowKeyPressFcn(hObject, eventdata, handles)
gdata = guidata(gcf);

switch eventdata.Key
    case 'uparrow'
        delta = [0,-1];
    case 'downarrow'
        delta = [0,1];
    case 'leftarrow'
        delta = [-1,0];
    case 'rightarrow'
        delta = [1,0];
end

if gca == gdata.ImageAx
    setConstrainedPosition(gdata.DiffPoint,getPosition(gdata.DiffPoint)+delta)
elseif gca == gdata.CbedAx
    setConstrainedPosition(gdata.RealPoint,getPosition(gdata.RealPoint)+delta)
end
guidata(gcf,gdata)
end

function WindowButtonDownFcn(hObject, eventdata, handles)
gdata = guidata(gcf);

pos = eventdata.IntersectionPoint(1:2);

if gca == gdata.ImageAx
    setPosition(gdata.DiffPoint,pos)
elseif gca == gdata.CbedAx
    setPosition(gdata.RealPoint,pos)
end

guidata(gcf,gdata)
end

function DiffPointCallback( pos )
gdata = guidata(gcf);
t = gdata.Transf;

PointPos = round(pos);
gdata.Cbed.CData = t(squeeze(gdata.d4d(:,:,PointPos(2),PointPos(1))));
title(gdata.ImageAx,sprintf('Real Space (%d, %d)',PointPos(2),PointPos(1)))

guidata(gcf,gdata)
end

function RealPointCallback( pos )
gdata = guidata(gcf);

PointPos = round(pos);
gdata.Image.CData = squeeze(gdata.d4d(PointPos(2),PointPos(1),:,:));
title(gdata.CbedAx,sprintf('Diffraction Space (%d, %d)',PointPos(2),PointPos(1)))

guidata(gcf,gdata)
end

function [] = MakeUiControls()
gdata = guidata(gcf);
bg = uibuttongroup('Title','Diffraction Space Contrast',...
                  'Visible','off',...
                  'Position',[.6 0.01 .28 .15],...
                  'FontWeight','bold',...
                  'FontSize',16,...
                  'BackgroundColor','w',...
                  'SelectionChangedFcn',@bselection);
              

              
gdata.s1 = uicontrol(bg,'Style','slider',...
                  'Position',[60 10 100 20],...
                  'BackgroundColor','w',...
                  'Min',0.001,'Max',1,...
                  'Value',0.001,...
                  'Callback',@bselection,...
                  'HandleVisibility','off');     

uicontrol(bg,'Style','text',...
                  'String','log',...
                  'FontSize',16,...
                  'HorizontalAlignment','right',...
                  'BackgroundColor','w',...
                  'Position',[0 10 50 25],...
                  'HandleVisibility','off');
uicontrol(bg,'Style','text',...
                  'String','linear',...
                  'FontSize',16,...
                  'HorizontalAlignment','left',...
                  'BackgroundColor','w',...
                  'Position',[170 10 150 25],...
                  'HandleVisibility','off');

bg.Visible = 'on';
guidata(gcf,gdata)
end


function [  ] = resizeFig( F,scale1,scale2 )

p=get(F,'position');
set(F,'position',[p(1:2),p(3)*scale1,p(4)*scale2])

end

function varargout=plotIM(varargin)
% Specialized version of imagesc using nice-looking defaults (grayscale, no
% ticks if scale undefined, proper aspect ratio, big font) for  electron 
% microscopy images. Takes an argument string in the same form as imagesc.

im = imagesc(varargin{:});

if ~isvector(varargin{1})
    %if x and y scales are not specified, don't show axis ticks
    set(gca,'XTick',[],'YTick',[]);
end
set(gca,'FontSize',18,'FontWeight','Bold');
axis image;
colormap(gray);

if nargout
    varargout{1} = im;
end

end
