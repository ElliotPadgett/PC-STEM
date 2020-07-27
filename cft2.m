function [ F ] = cft2( f, q1, q2 ,zeroCentered)
%cft2 2D continuous Fourier tranform of f evaluated at point q1, q2
%   inputs:
%       f -- the 2D array the fourier transform is calculated from
%       q1,q2 -- indices where the transform is to be evaluated, following
%              the same convention as fft2.  q1,q2 can be non-integers.
%       zeroCentered -- boolean indicating the q index corresponding to
%                       zero: 0 - default, zero is at index 1,1 (same as
%                                 fft2(f))
%                             1 - zero is at the image center,
%                                 corresponding to fftshift(fft2(f))
%   outputs:
%       F -- value of the fourier transform of f at q1,q2.  This is a complex
%            number, rather than an array.
%
%This function is part of the PC-STEM Package by Elliot Padgett in the 
%Muller Group at Cornell University.  Last updated by Megan Holtz for speed July 27, 2020.

if nargin<4
    zeroCentered = 0;
end

[m,n]=size(f);
jgr = 0:m-1;
kgr = 0:n-1;

if zeroCentered
    q1=q1+m/2;
    q2=q2+n/2;
end

F = sum(sum( f.* (exp(-2*pi*1i*jgr'*(q1-1)/m)*exp(-2*pi*1i*kgr*(q2-1)/n)) ));

end

