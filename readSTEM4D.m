function [ data ] = readSTEM4D( fnm )
%readSTEM4D Loads 4D EMPAD data in order k1,k2,x1,x2
% input:
%   fnm -- name of raw data file to read.  x and y sizes should be listed
%          at end, e.g. 'scan_1_x256_y256.raw' is a valid file name. This
%          function assumes a data is square in real space.
% output:
%   data -- 4D data in order k1, k2, x1, x2 in single data type.
%
%This function is part of the PC-STEM Package by Elliot Padgett in the 
%Muller Group at Cornell University.  Last updated June 21, 2019

fprintf('Loading file: %s \n',fnm), tic
% Parse filename
fnm_parts=strsplit(fnm,'_');
x_size=str2double(fnm_parts{end}(2:end-4)); %real space size
k1_size = 128; %Assumed diffraction space sizes for EMPAD data
k2_size = 130;

%read data
data = single(fread(fopen(fnm,'r'),k1_size*k2_size*x_size*x_size,'float32'));

%make 4D array
data = reshape(data,k1_size,k2_size,x_size,x_size);

%Crop off junk rows
data = data(3:126,3:126,:,:);

%Set data orientation
data = permute(data,[1,2,4,3]);
data = flip(data,3);
data = flip(data,2);

fprintf('Data loaded successfully.\n'), toc
end

