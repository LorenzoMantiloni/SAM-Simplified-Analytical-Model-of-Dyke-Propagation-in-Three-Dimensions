function [S3,S3dirx,S3diry,S3dirz] = HalfSpaceNormalForce(X,Y,Z,P,Pr,...
    nu,Xloc,Yloc)

%  "HalfSpaceNormalForce" calculates the magnitude and direction of sigma_3
%  at arbitrary points in a half-space, in the presence of a distribution
%  of normal forces on the free surface. The distribution simulates  
%  surface loading/unloading (due to e.g. a volcanic edifice or a caldera).
%  Solutions for the stress components due to a single normal force of
%  intensity N are taken by Jaeger et al. (2007), cap. 13.5 (see reference
%  below). 
%
%  Inputs: 
%          X,Y,Z:      coordinates of observation points
%          P:          loading/unloading pressure (Pa) (e.g. for a caldera
%                      of depth h, P = -Pr*g*h)
%          Pr:         rock density (kg/m^3)
%          nu:         Poisson's ratio
%          Xloc,Yloc:  x,y coordinates of individual normal forces on the 
%                      free surface.
%
%  Reference: 
% 
%  Jaeger, J. C., Cook, N. G. W., & Zimmermann, R. W. (2007). Fundamentals 
%  of rock mechanics, edited by Chapman and Hall. New York, United States 
%  of America, 1-593.
% 
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  Author: Lorenzo Mantiloni
%  Copyright (c) Lorenzo Mantiloni, 2023 ,Deutsche GeoForschungsZentrum 
%  GFZ\University of Potsdam
%
%  Permission is hereby granted, free of charge, to any person obtaining a 
%  copy of this software and associated documentation files 
%  (the "Software"), to deal in the Software without restriction, including 
%  without limitation the rights to use, copy, modify, merge, publish, 
%  distribute, sublicense, and/or sell copies of the Software, and to 
%  permit persons to whom the Software is furnished to do so, subject to 
%  the following conditions:
% 
%  1) The above copyright notice and this permission notice shall be  
%  included in all copies or substantial portions of the Software.
% 
%  2) The Software is provided "as is", without warranty of any kind, 
%  express or implied. Under no circumstances shall the authors or
%  copyright holders be liable for any claim, damages or other liability
%  arising from, out of or in connection with the Software or the use or
%  other dealings in the Software.
%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g = 9.81; %Acceleration of gravity (m/s^2)

if numel(Xloc) == 1
    %  The loading/unloading pressure coincides with the intensity of the 
    %  single normal force
    N = P; 
else
    %  The loading/unloading pressure is equally distributed among all
    %  normal forces
    N = P/(pi*(0.5*abs(max(Xloc)-min(Xloc)))^2);
end

C = N/(2*pi);

X = X(:);
Y = Y(:);
Z = -Z(:);

Sxx = -Pr*g*Z; Syy = -Pr*g*Z; Szz = -Pr*g*Z; Sxy = 0; Sxz = 0; Syz = 0;

for i=1:numel(Xloc)
    Xshift = X - Xloc(i);
    Yshift = Y - Yloc(i);
    R = (Xshift.^2 + Yshift.^2 + Z.^2).^(1/2);
    X2 = Xshift.^2;
    Y2 = Yshift.^2;
    Z2 = Z.^2;
    Sxx = C*(3*X2.*Z./R.^5 + (1-2*nu)*(Y2 + Z2)./((Z + R).*R.^3) - ...
        (1-2*nu)*Z./R.^3 - (1-2*nu)*X2./(((Z + R).^2).*R.^2)) + Sxx;
    Syy = C*(3*Y2.*Z./R.^5 + (1-2*nu)*(X2 + Z2)./((Z + R).*R.^3) - ...
        (1-2*nu)*Z./R.^3 - (1-2*nu)*Y2./(((Z + R).^2).*R.^2)) + Syy;
    Szz = C*(3*Z.^3)./R.^5 + Szz;
    Sxy = C*(3*Xshift.*Yshift.*Z./R.^5 - (1-2*nu)*Xshift.*Yshift.*...
        (Z + 2*R)./(((Z + R).^2).*R.^3)) + Sxy;
    Sxz = C*(3*Xshift.*Z2./R.^5) + Sxz;
    Syz = C*(3*Yshift.*Z2./R.^5) + Syz;
end

[S3,~,~,S3dir,~,~] = EigCalc3d(Sxx,Syy,Szz,Sxy,Sxz,Syz);

S3dirx = S3dir(:,1);
S3diry = S3dir(:,2);
S3dirz = S3dir(:,3);




