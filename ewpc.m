function [ FFTdata ] = ewpc( data , useWindow)
%ewpc calculates the EWPC transform fft(log(data)) for 2d or 4d data
%   input:
%       data -- 2d or 4d scanning electron diffraction data, ordered (kx,
%               ky, x, y)
%       useWindow -- logical indicating whether or not to apply hann window
%               before FFT.  Default is 1.  The window is useful to prevent
%               FFT artifacts from non-periodicity.  This is especially
%               important it the diffraction patterns have significant
%               intensity at their edges.

if nargin==1
    useWindow = 1;
end

[N_kx,N_ky,N_x,N_y]=size(data);
minval=min(data(:));

if useWindow
    win=window2(N_kx,N_ky,@hann);
else
    win=ones(N_kx,N_ky);
end

fftMag = @(x) abs(fftshift(fft2(win.*log(x-minval+0.1))));

% Convert to FFT of CBED map
FFTdata = data;
for x=1:N_x
    for y=1:N_y
        FFTdata(:,:,x,y) = fftMag(FFTdata(:,:,x,y));
    end
end

end

function w=window2(N,M,w_func)

wc=window(w_func,N);
wr=window(w_func,M);
[maskr,maskc]=meshgrid(wr,wc);

w=maskr.*maskc;

end
