function varargout=plotIM(varargin)
% Specialized version of imagesc using nice-looking defaults (grayscale, no
% ticks if scale undefined, proper aspect ratio, big font) for  electron 
% microscopy images. Takes an argument string in the same form as imagesc.

im = imagesc(varargin{:});

if ~isvector(varargin{1})
    %if x and y scales are not specified, don't show axis ticks
    set(gca,'XTick',[],'YTick',[]);
end
set(gca,'FontSize',16,'FontWeight','Bold');
axis image;
colormap(gray);

if nargout
    varargout{1} = im;
end

end
