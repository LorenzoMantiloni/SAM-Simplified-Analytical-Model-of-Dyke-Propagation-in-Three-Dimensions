function[S3,S3dir] = Retrieve_sigma3(X,Y,Z,TectS,Dds,Dss,Dn,P1,P2,P3,Pr,...
    nu,mu,g)

%  "Retrieve_sigma3" calculates the magnitude and direction of sigma_3 at
%  points of coordinates X,Y,Z. First, it evaluates Sij components of the
%  cartesian stress tensor through the Boundary Element software
%  "Cut&Displace" (Davis et al., 2017; 2019, see references below), then it
%  diagonalizes the stress tensor and retrieves the sigma_3 eigenvector. 
%  For further details on inputs and outputs, see the Instruction Manual,
%  as well as the inline comments.
%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  References: 
%
%  Davis, T., Healy, D., Bubeck, A., & Walker, R. (2017). 
%  Stress concentrations around voids in three dimensions: 
%  The roots of failure. Journal of Structural Geology, 102, 193-207.
% 
%  Davis, T., Healy, D., & Rivalta, E. (2019). Slip on wavy frictional 
%  faults: Is the 3rd dimension a sticking point?. 
%  Journal of Structural Geology, 119, 33-49.
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

TectSxx = TectS(1);
TectSyy = TectS(2);
TectSxy = TectS(3);

%  If halfspace = 1, the calculations are performed in a half-space. This  
%  is not the case when applying Martel & Muller (2000).
halfspace = 0;  

lambda = 2*mu*nu/(1-2*nu);

%  Calculating cartesian stress components at observation points X,Y,Z
[StressTTotal,~,~,~,~,~]=CalculateStressOnSurroundingPoints3d(Dss,Dds,Dn,mu,lambda,X,Y,Z,TectSxx,...
TectSyy,0,TectSxy,0,0,P1,P2,P3,halfspace,nu);

[Sxx,Syy,Szz,Sxy,Sxz,Syz] = ExtractCols(StressTTotal);

%  Superimposing the background (lithostatic) stress
condfreesurface = Z > 0;
SxxG = Pr*g*Z;	
SyyG = Pr*g*Z;		
SzzG = Pr*g*Z;
SxxG(condfreesurface) = 0;
SyyG(condfreesurface) = 0;
SzzG(condfreesurface) = 0;

Fctr=1;
Sxx=Sxx+(SxxG(:).*Fctr);
Syy=Syy+(SyyG(:).*Fctr);
Szz=Szz+(SzzG(:).*Fctr);

[S3,~,~,S3dir,~,~] = EigCalc3d(Sxx(:),Syy(:),Szz(:),Sxy(:),Sxz(:),Syz(:));

