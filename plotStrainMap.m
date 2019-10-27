function [] = plotStrainMap(StrainComponents)
%ewpc_mapSpots Produces a figure showing the fitted EWPC peaks
%   input:
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
%
%This function is part of the PC-STEM Package by Elliot Padgett in the 
%Muller Group at Cornell University.  Last updated June 26, 2019.

figure,
subplot(2,2,1),plotIM(StrainComponents.Eps11,'AlphaData',~isnan(StrainComponents.Eps11))
    colorbar,title('\epsilon_1_1')
subplot(2,2,2),plotIM(StrainComponents.Eps22,'AlphaData',~isnan(StrainComponents.Eps22))
    colorbar,title('\epsilon_2_2')
subplot(2,2,3),plotIM(StrainComponents.Eps12,'AlphaData',~isnan(StrainComponents.Eps12))
    colorbar,title('\epsilon_1_2')
subplot(2,2,4),plotIM(StrainComponents.Theta,'AlphaData',~isnan(StrainComponents.Theta))
    colorbar,title('\theta')
colormap parula

figure,
subplot(1,3,1),plotIM(StrainComponents.majAx,'AlphaData',~isnan(StrainComponents.majAx))
    colorbar,title('Major axis')
subplot(1,3,2),plotIM(StrainComponents.minAx,'AlphaData',~isnan(StrainComponents.minAx))
    colorbar,title('Minor axis')
subplot(1,3,3),plotIM(StrainComponents.strainAngle,'AlphaData',~isnan(StrainComponents.strainAngle))
    colorbar,title('Axis angle')

colormap parula
end