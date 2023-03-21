function[dike_path_tip,dike_path_dip,dike_path_strike,dike_path_K] = ...
    SAM_Forward_Interp(StartPoint,c,O,Pm,TopoInterp,MDT,options,...
    InterpOption,varargin)

%  "SAM_Forward_Analytical" simulates dike trajectories from starting
%  points in the subsurface up to either the free surface or a user-defined
%  minimum distance threshold (MDT) from it. Stresses are evaluated through 
%  interpolating functions for either the individual components of the 
%  cartesian stress tensor, or the magnitude and direction components of 
%  sigma_3. For further details on inputs and outputs, see the Instruction 
%  Manual, as well as the inline comments.
%  
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Author: Lorenzo Mantiloni
%  Copyright (c) Lorenzo Mantiloni, 2023, Deutsches GeoForschungsZentrum 
%  GFZ - University of Potsdam
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

%Store the "original" observation points:
O_orig = O;

%  Getting some additional parameters controlling whether to plot the dike 
%  trajectory or not, as well as the lateral distance and depth at which 
%  trajectories are cut off 
if ~isempty(options)
    display = options.behavior(1);
    RigidSAM = options.behavior(2);
    displast = options.behavior(3);
    if isfield(options,'cutoff')
        Traj_RadialCutoff = options.cutoff(1);
        Traj_VertCutoff = options.cutoff(2);
    else
        Traj_RadialCutoff = [];
        Traj_VertCutoff = [];
    end
    if isfield(options,'backtracking')
        backtracking = options.backtracking;
    else
        backtracking = 0;
    end
    if isfield(options,'ProjectToFreeSurface')
        ToFreeSurface = options.ProjectToFreeSurface(1);
        if numel(options.ProjectToFreeSurface) > 1
            PathStep = options.ProjectToFreeSurface(2);
        else
            PathStep = 50;
        end
    else
        ToFreeSurface = 0;
    end
else
    display = 1;
    Traj_RadialCutoff = [];
    Traj_VertCutoff = [];
    RigidSAM = 0;
    displast = 0;
    ToFreeSurface = 0;
end

%  Importing the interpolating functions for the external stress
if strcmp(InterpOption,'Cartesian')
    Cartesian = 1;
    SxxInterp = varargin{1}; SyyInterp = varargin{2}; SzzInterp = varargin{3};
    SxyInterp = varargin{4}; SxzInterp = varargin{5}; SyzInterp = varargin{6};
    %  Checking if dike volume V and rock fracture toughness Kc are provided
    if numel(varargin) > 6
        DikeVol = 1;
        V = varargin{7};
        Kc = varargin{8};
        mu = varargin{9};
        nu = varargin{10};
    else
        DikeVol = 0;
    end
else
    Cartesian = 0;
    S3Interp = varargin{1}; S3xInterp = varargin{2}; 
    S3yInterp = varargin{3}; S3zInterp = varargin{4};
    %  Checking if dike volume V and rock fracture toughness Kc are provided
    if numel(varargin) > 5
        DikeVol = 1;
        V = varargin{5};
        Kc = varargin{6};
        mu = varargin{7};
        nu = varargin{8};
    else
        DikeVol = 0;
    end
end

%  Initializing some arrays
Norm_Stress = zeros(1,size(O,2));
Buoyancy = Norm_Stress;
Norm_Stress_Grad = Norm_Stress;
Buoyancy_Grad = Norm_Stress;
dike_path_tip = [];
dike_path_dip = [];
dike_path_strike = [];
dike_path_K = [];

%  Setting some parameters
condstop=0; %If == 1, the dike stops

%  Setting StartPoint as the first F. Dike centers, i.e. the points 
%  defining the trajectory, are called F in accordance with 
%  Mantiloni et al., 2023
F = StartPoint; 

%  Main loop of the function: it goes on until either the dike reaches the
%  free surface, turns back, stops horizontally, or crosses a radial or 
%  vertical threshold (if they exist). If RigidSAM == 1, the dike will also 
%  stop if it steers, dips or ascends too abruptly.

while ~condstop
    
    clear Ry Rz K %Make sure these variables are cleared
    
    %  Defining the "dike": a penny-shaped crack with radius c, centered at
    %  F, opening against the local direction of sigma_3, with n
    %  observation points along its tip-line.
    
    %  Calculating sigma_3 direction at F
    if Cartesian
        [~,S3dir] = Retrieve_sigma3_Interp(F(1),F(2),F(3),SxxInterp,...
            SyyInterp,SzzInterp,SxyInterp,SxzInterp,SyzInterp);
    else
        S3dir = [S3xInterp(F(1),F(2),F(3)) S3yInterp(F(1),F(2),F(3)) ...
            S3zInterp(F(1),F(2),F(3))];
    end
    S3dirmod = (S3dir(1)^2 + S3dir(2)^2 + S3dir(3)^2)^(1/2); %Modulus of the sigma3 eigenvector
    theta = acos(S3dir(3)/S3dirmod); %Dip angle of sigma3
    al = alpha_finder(S3dir,theta); %Strike angle of sigma3
    
    %  Rotation matrices to make the observation points ring aligned 
    %  perpendicularly to sigma3
    Ry = [cos(theta) 0 sin(theta);0 1 0;-sin(theta) 0 cos(theta)]; 
    Rz = [cos(al) -sin(al) 0;sin(al) cos(al) 0;0 0 1];
    
    %  Shifting and aligning the observation points ring to the current 
    %  dike location and orientation
    O = Ry*O; 
    O = Rz*O;
    O(1,:) = O(1,:) + F(1);O(2,:) = O(2,:) + F(2);O(3,:) = O(3,:) + F(3);

    %  Calculating sigma3 on each of the observation points O (and 
    %  subsequently the normal stress and the fluid's hydrostatic pressure 
    %  at each O). Normal stresses are simply equal to sigma3, since the
    %  crack's surface is perpendicular to it).
    for i=1:size(O,2)
        Ox = O(1,i);
        Oy = O(2,i);
        Oz = O(3,i);
        if Cartesian
            [S3,~] = Retrieve_sigma3_Interp(Ox,Oy,Oz,SxxInterp,...
                SyyInterp,SzzInterp,SxyInterp,SxzInterp,SyzInterp);
        else
            S3 = S3Interp(Ox,Oy,Oz);
        end
        Norm_Stress(i) = S3;
        Buoyancy(i) =-Pm*g*Oz;
    end

    %  Calculating the normal stress and buoyancy gradient (considering
    %  each O with its respective antipodal point and taking 2*c as their 
    %  distance)
    for i=1:numel(Norm_Stress)
        if i<=numel(Norm_Stress)/2
            Norm_Stress_Grad(i) = (Norm_Stress(i) - ...
                Norm_Stress(i+numel(Norm_Stress)/2))/(2*c);
            Buoyancy_Grad(i) = (Buoyancy(i) - ...
                Buoyancy(i+numel(Norm_Stress)/2))/(2*c);
        else
            Norm_Stress_Grad(i) = (Norm_Stress(i) - ...
                Norm_Stress(i-numel(Norm_Stress)/2))/(2*c);
            Buoyancy_Grad(i) = (Buoyancy(i) - ...
                Buoyancy(i-numel(Norm_Stress)/2))/(2*c);
        end
    end

    %  Calculating the stress intensity factor at each observation point
    if DikeVol
        K = 3*mu*V/(4*(1-nu)*(c^2)*sqrt(pi*c)) + (4/(3*pi))*c*((pi*c)^(1/2)).*(Norm_Stress_Grad + Buoyancy_Grad);
    else
        K = (4/(3*pi))*c*((pi*c)^(1/2)).*(Norm_Stress_Grad + Buoyancy_Grad);
    end
    dike_path_K = vertcat(dike_path_K,K);
    %  Determining O where K is maximum and, subsequently, the next F 
    %  ("Fn")
    condK = find(K==max(K));
    Fn = [O(1,condK(1)),O(2,condK(1)),O(3,condK(1))]; 
    %  If the dike volume is also provided, checking if Kmax/Kc > 1. If 
    %  not, the dike will stop.
    if DikeVol && max(K)/Kc < 1
       condstop = 1;
    end
    
    %  Determining whether any O has hit the free surface and stop the
    %  dike if it does
    zmax = max(O(3,:));
    condzmax = find(O(3,:)==zmax);
    xzmax = O(1,condzmax(1));
    yzmax = O(2,condzmax(1));
    if zmax > TopoInterp(xzmax,yzmax) - MDT
       condstop = 1;
    end
    
    %  Checking if the current F falls outside the radial threshold, if it 
    %  exists
    if ~isempty(Traj_RadialCutoff)
        if (Fn(1)^2 + Fn(2)^2)^(1/2) > Traj_RadialCutoff
            condstop = 1;
        end
    end
    
    %  Checking if the current F falls below the vertical threshold, if it 
    %  exists
    if ~isempty(Traj_VertCutoff)
        if Fn(3) < Traj_VertCutoff
            condstop = 1;
        end
    end
    
    %  Check if the dike is turning back and stop if it does
    if numel(dike_path_dip) >= 2
        direction_old = F - dike_path_tip(end,:);
        direction_current = Fn - F;
        if dot(direction_current,direction_old) < 0
            condstop = 1;
        end
        %  If RigidSAM == 1, this checks if the dike steers, dips or 
        %  ascends abruptly (fixed threshold on the change in the strike 
        %  and dip angles between the current and the last trajectory step)
        if RigidSAM
            if abs(dot(direction_current,direction_old)) < 0.75
                condstop = 1;
            end
        end
    end

    %  Updating the trajectory and setting Fn as the current F
    dike_path_tip = vertcat(dike_path_tip,F);
    dike_path_dip = [dike_path_dip theta];
    dike_path_strike = [dike_path_strike al];
    F = Fn;
    if display && ~condstop
        patch = fill3(O(1,:),O(2,:),O(3,:),'r');
        patch.FaceAlpha = 0.5;
        pause(0.1)
        hold on
    end
    if displast && condstop
        patch = fill3(O(1,:),O(2,:),O(3,:),'y');
        patch.FaceAlpha = 0.5;
        pause(0.1)
        hold on
    end
    
    O = O_orig; %Set the observation points back to the original ones
    
    if backtracking && numel(dike_path_dip) == 1
        %  If we are using the function in SAM_Backtrack_Analytical and the 
        %  algorithm stops at the first step, this makes sure dike_path_tip 
        %  includes a point different from the starting one:
        dike_path_tip = vertcat(dike_path_tip,Fn); 
    end
    if backtracking && numel(dike_path_dip) == 2 
        %  If we are using the function in SAM_Backtrack_Analytical, this 
        %  stops the algorithm at the second step.
        condstop = 1;
    end
    
end

if ToFreeSurface && numel(dike_path_strike) > 1
    F = dike_path_tip(end,:);
    theta = dike_path_dip(end-1);
    al = dike_path_strike(end-1);
    while F(3) < TopoInterp(F(1),F(2)) && F(3) > -2*MDT 
        F(1) = F(1) - PathStep*cos(theta)*cos(al);
        F(2) = F(2) - PathStep*cos(theta)*sin(al);
        F(3) = F(3) + PathStep*sin(theta);
    if display 
        scatter3(F(1),F(2),F(3),'o','b');
    end
       dike_path_tip = vertcat(dike_path_tip,F);
    end
end





