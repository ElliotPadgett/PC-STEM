function [  ] = writeSTEM4D( data, fnmstub )
%readSTEM4D Writes 4D EMPAD data in orinal-equivalent raw format
% input:
%   data -- 4D data in order k1, k2, x1, x2 in single data type.
%   fnmstub -- stub for filename.  full filename will include parameters
%       for reading raw file.
%
%This function is part of the PC-STEM Package by Elliot Padgett in the 
%Muller Group at Cornell University.  Last updated July 18, 2019

data = flip(data,2);
data = flip(data,3);
data = permute(data,[1,2,4,3]);

%Add on junk rows
data = padarray(data,[2,2,0,0],0,'pre');
data = padarray(data,[2,4,0,0],0,'post');

writeraw(data,fnmstub);

end

function fullfilename = writeraw(II,varargin)
% writeraw, writeraw(II) or writeraw(II,filename)
% this function write a input matrix into a binary format .dat file
% II = matrix that u want to write out
% filename = the name u want to name it
% by Huolin Xin

if nargin==0
   varlist = evalin('base','who');
   varselnum = listdlg('ListString',varlist);
   varname = varlist{varselnum};
   II = evalin('base',varname);
   filename = inputdlg('Input the name for this file (w/o suffix):');
   filename = filename{1};
elseif nargin==1
   filename = inputdlg('Input the name for this file (w/o suffix):');
   filename = filename{1};
else
   filename = varargin{1};
end


[pathstr, name, ext] = fileparts(filename);
if isempty(pathstr)
   dirname = pwd;
else
   dirname = '';
end
suffix = [];
matsize = size(II);
for i=1:length(matsize)
   if i==1
       suffix = num2str(matsize(i));
   else
       suffix = [suffix,'x',num2str(matsize(i))];
   end
end

% Added ending for compatibility with readSTEM4D -- Elliot
fullfilename=fullfile(dirname,sprintf('%s_%s_float_x%d_y%d.raw',filename,suffix,matsize(3),matsize(4)));

fp=fopen(fullfilename,'w');
fwrite(fp,II,'float');
fclose(fp);
end