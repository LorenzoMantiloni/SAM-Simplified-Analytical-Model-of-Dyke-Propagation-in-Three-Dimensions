function[al] = alpha_finder(S3dir,theta)

%  "alpha_finder" retrieves the strike angle (al) of a crack oriented
%  perpendicularly to the direction of sigma_3 evaluated at its centre, 
%  knowing the dip angle (theta) of the crack's plane. "al" is calculated
%  counter-clockwise away from the positive x-axis.
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

S3dirmod = (S3dir(1)^2 + S3dir(2)^2 + S3dir(3)^2)^(1/2);

            if S3dir(1) == 0
                if theta == 0
                    al = 0;
                else
                    if S3dir(2) == 0
                        al = acos(S3dir(1)/(S3dirmod*sin(theta)));
                    else
                        al = sign(S3dir(2))*acos(S3dir(1)/(S3dirmod*...
                            sin(theta)));
                    end
                end
            else
                if S3dir(2) == 0
                    al = acos(S3dir(1)/(S3dirmod*sin(theta)));
                else
                    al = sign(S3dir(2))*acos(S3dir(1)/(S3dirmod*...
                        sin(theta)));
                end
            end
            if imag(al)~=0
                al=real(al);
            end

end