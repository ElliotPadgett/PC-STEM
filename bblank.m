function [ im ] = bblank( im, bsize, shape, center, val )
%bblank adds "beam blank" to image im
%   input:
%       im -- 2d image to add blanked region to
%       bsize -- blank spot radius or edge half-length in pixels (default 2)
%       shape -- square 's' or circle 'c' (optional, default 's')
%       val -- value to fill in blanked region (optional, default 0)
%   output: 
%       blankedim -- image with beam blank added
%
%This function is part of the PC-STEM Package by Elliot Padgett in the 
%Muller Group at Cornell University.  Last updated June 26, 2019.

%default values
switch nargin
    case 4
        val = 0;
    case 3
        val = 0; center = 0.5*(size(im)+1);
    case 2
        val = 0; center = 0.5*(size(im)+1); shape = 's';
    case 1
        val = 0; center = 0.5*(size(im)+1); shape = 's'; bsize = 2;
end

Nd=ndims(im);

switch shape
    case 's'
        if Nd==2
        im(round(center(1)-bsize):round(center(1)+bsize),...
                  round(center(2)-bsize):round(center(2)+bsize)) = val;
        elseif Nd==4
        im(round(center(1)-bsize):round(center(1)+bsize),...
                  round(center(2)-bsize):round(center(2)+bsize),:,:) = val;
        end
    case 'c'
        [N1,N2,N3,N4] = size(im);
        [X,Y] = meshgrid(1:N1,1:N2);
        blank = repmat((X-center(1)).^2 + (Y-center(2)).^2 <= bsize^2,1,1,N3,N4);
        im = im .* single(~blank);
        im = im + val*blank;
    otherwise
        error('bblank: Unrecognized shape')
end


end

