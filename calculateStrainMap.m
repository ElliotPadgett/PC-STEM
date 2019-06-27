function [ StrainComponents,StrainTensors] = calculateStrainMap( spotMaps, spotReferences, latticeCoords )
%StrainTensorMap Makes maps of strain tensor components from peak maps
%   inputs:
%       spotMaps -- struct array containing maps of maximum spot index Q1,Q2
%                 for each spot in spotList. Fields are:
%                 'id' -- a name identifying the spot.
%                 'Q1map' -- map of Q1 index value at peak maximum.
%                 'Q2map' -- map of Q2 index value at peak maximum.
%                 'x1map' -- map of x1 vector component for the spot.
%                 'x2map' -- map of x2 vector component for the spot.
%                 'spotlength' -- map of the spot vector length.
%                 'spotangle' -- map of the spot vector angle in degrees.
%       spotReferences -- struct array containing reference points
%                         corresponding to points in Qmaps. Fields are:
%                         'id' -- a name identifying the spot.
%                         'point' -- reference spot location [q1,q2]
%       latticeCoords -- choose coordinate system to reference
%                        strain to.  Value 0 selects "image 
%                        coordinates", with the 1,2 directions 
%                        corresponding to image x1,x2.  Value 1 selects
%                        "lattice coordinates", where the 1,2 directions
%                        refer to the reference directions.
%   outputs:
%       StrainComponents -- struct containing 2D arrays (images) of the
%               strain tensor components, rotation, and strain ellipse 
%               parameters, including:
%                   Eps11 -- The first diagonal element of the strain tensor
%                   Eps22 -- The second diagonal element of the strain tensor
%                   Eps12 -- The 2D shear component of the strain tensor (1,2) off
%                            diagonal
%                   Theta -- The rotation relative to the reference, decomposed from
%                            the strain tensor, in degrees
%                   minAx -- semi-minor axis length of the strain ellipse
%                   majAx -- semi-major axis length of the strain ellipse
%                   strainAngle -- angle of the strain ellipse semi-major
%                                  axis
%       StrainTensors -- struct array containing the strain tensor 'E' and
%                        rotation matrix 'R' describing the distortion at each
%                        real space position, along with the eigenvalues
%                        'evals' and eigenvectors 'evecs' of the strain tensor
%
%This function is part of the PC-STEM Package by Elliot Padgett in the 
%Muller Group at Cornell University.  Last updated June 26, 2019.

%initialize data info
[Nx1,Nx2] = size(spotMaps(1).Q1map);

%prealocate results structs
StrainComponents = struct('Eps11', nan(Nx1,Nx2), 'Eps22', nan(Nx1,Nx2),...
                          'Eps12', nan(Nx1,Nx2), 'Theta', nan(Nx1,Nx2),...
                          'majAx', nan(Nx1,Nx2), 'minAx', nan(Nx1,Nx2),...
                          'strainAngle', nan(Nx1,Nx2));
StrainTensors = struct('E',cell(Nx1,Nx2),'R',cell(Nx1,Nx2),...
    'evals',cell(Nx1,Nx2),'evecs',cell(Nx1,Nx2));

%prepare reference point list
referencePoints = [0,0];
for s = 1:length(spotReferences)
    referencePoints = [referencePoints; spotReferences(s).point];
end

%calculate strain across map
for j=1:Nx1
    for k=1:Nx2
        
        dataPoints = [0,0];
        for s = 1:length(spotReferences)
            %center spot
            q1c = spotMaps(s).VectorX1(j,k);
            q2c = spotMaps(s).VectorX2(j,k);
            
            %include in list for tranformation calculation
            dataPoints = [dataPoints; [q1c,q2c]];
        end
        if sum(isnan(dataPoints(:)))
            %nan values mean values are unknown and this point should
            %be skipped

            StrainTensors(j,k).R = nan(2);
            StrainTensors(j,k).E = nan(2);
            StrainTensors(j,k).evecs = nan(2,1);
            StrainTensors(j,k).evals = nan(2,1);

        else
            %Calculate transform
            %tform = fitgeotrans(movingPoints,fixedPoints,'affine');
            % we want to know what it would be like if we mapped the
            % reference lattice onto our actual data
            tform = fitgeotrans(referencePoints,dataPoints,'affine');
            D=tform.T; D = D(1:2,1:2);
            
            [R,U,V] = poldecomp(D); % D = R*U = V*R.  R is a rotation.
            %V is the strain in "world coordinates" and is independent of
            %the reference orientation. U is the strain in "lattice
            %coordinates", corresponding to the reference directions
            
            if latticeCoords
                E = U-eye(2); % strain tensor
            else
                E = V-eye(2); % strain tensor
            end
            
            StrainTensors(j,k).R = R;
            StrainTensors(j,k).E = E;
            
            %get strain tensor elements
            StrainComponents.Eps11(j,k) = E(1,1); %fractional stretch
            StrainComponents.Eps22(j,k) = E(2,2); %fractional stretch
            StrainComponents.Eps12(j,k) = E(1,2); %fractional shear
            StrainComponents.Theta(j,k) = atan2d(R(2,1),R(1,1)); %rotation in degrees
            
            %get strain elipse parameters
            [v,d]=eig(E);
            [evals,sorting] = sort([d(1,1),d(2,2)]); %sorts in ascending order
            evecs = v(:,sorting);
            StrainComponents.minAx(j,k) = evals(1); 
            StrainComponents.majAx(j,k) = evals(2); 
            StrainComponents.strainAngle(j,k) = atand(evecs(2,2)/evecs(1,2));
            
            StrainTensors(j,k).evecs = evecs(:,2);
            StrainTensors(j,k).evals = evals';
            
        end
    end
end

end


function [R, U, V] = poldecomp(F)
%POLDECOMP  Performs the polar decomposition of a regular square matrix.
%   [R U V] = POLDECOMP(F) factorizes a non-singular square matrix F such
%   that F=R*U and F=V*R, where
%   U and V are symmetric, positive definite matrices and
%   R is a rotational matrix
%
%   See also EIG, DIAG, REPMAT


% This kind of decomposition is often used in continuum mechanics so it is
% convenient to comment the code that way. From now, we use the matrix
% formalism of tensors. C is the right Cauchy-Green deformation tensor,
% F is the deformation tensor, lambda is the stretch.
%
%Copyright (c) 2014, Zoltán Csáti

% Check input
[m n] = size(F);
if m ~= n
    error('Matrix must be square.');
end

C = F'*F;
[Q0, lambdasquare] = eig(C);
lambda = sqrt(diag((lambdasquare))); % extract the components
% Uinv is the inverse of U and is constructed with the help of Q0. Uinv is
% produced in the same base as F not in the base of its eigenvectors.
Uinv = repmat(1./lambda',size(F,1),1).*Q0*Q0';
% Using the definition, R, U and V can now be calculated
R = F*Uinv;
U = R'*F;
V = F*R';
end