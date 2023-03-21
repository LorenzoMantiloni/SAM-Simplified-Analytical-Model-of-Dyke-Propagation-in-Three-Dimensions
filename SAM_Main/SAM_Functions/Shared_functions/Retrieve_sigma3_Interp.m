function[S3,S3dir] = Retrieve_sigma3_Interp(X,Y,Z,SxxInterp,SyyInterp,...
    SzzInterp,SxyInterp,SxzInterp,SyzInterp)

%  "Retrieve_sigma3_Interp" calculates the magnitude and direction of 
%  sigma_3 at points of coordinates X,Y,Z. First, it evaluates Sij 
%  components of the cartesian stress tensor by interpolating functions 
%  (SijInterp), then it diagonalizes the stress tensor and retrieves the 
%  sigma_3 eigenvector. For further details on inputs and outputs, see the 
%  Instruction Manual, as well as the inline comments.
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

Sxx = SxxInterp(X(:),Y(:),Z(:));
Syy = SyyInterp(X(:),Y(:),Z(:));
Szz = SzzInterp(X(:),Y(:),Z(:));
Sxy = SxyInterp(X(:),Y(:),Z(:));
Sxz = SxzInterp(X(:),Y(:),Z(:));
Syz = SyzInterp(X(:),Y(:),Z(:));   

[S3,~,~,S3dir,~,~] = EigCalc3d(Sxx(:),Syy(:),Szz(:),Sxy(:),Sxz(:),Syz(:));

