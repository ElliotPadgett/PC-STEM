function [ data4d_binSpace ] = binSpace( data4d, scale)
%binSpace scales 4d cbed data to bin in spatial dimensions only
%   input:
%       data4d -- 4D STEM dataset with dimensions ordered [k1,k2,x1,x2].
%       scale -- real-space scaling factor, equal to 1/binning.
%   output: 
%       data4d_binSpace -- re-scaled 4D dataset.
%
% Elliot Padgett -- Muller Group, Cornell University -- October 2019

fprintf('Binning by %d...\n',1/scale)

[N_kx,N_ky,N_x,N_y]=size(data4d);

data4d_binSpace=zeros([N_kx,N_ky,N_x*scale,N_y*scale]);
for i=1:N_kx
    for j=1:N_ky
        data4d_binSpace(i,j,:,:)=imresize(squeeze(data4d(i,j,:,:)),scale,'bilinear');
    end
end

fprintf('Done binning.\n')
end

