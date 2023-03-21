function[dike_path_tip,dike_path_dip,dike_path_strike,BSP] = ...
    SAM_Backtrack_Analytical(StartPoint,cB,O,Pr,Pm,TopoInterp,...
    Stressfun,Stressfun_Par,options,varargin)

%  "SAM_Backtrack_Analytical" backtracks dike trajectories from starting
%  points, either on the surface or in the subsurface, down to a potential
%  magma storage volume. Stresses are evaluated through an analytical
%  function provided by the user via a function handle. For further details
%  on inputs and outputs, see the Instruction Manual, as well as the inline
%  comments.
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

%  Importing the analytical function for stresses and the known parameters
Narg = nargin(Stressfun);

if numel(Stressfun_Par) + 3 < Narg
    error(['The number of inputs of Stressfun and the number of '...
        'provided inputs do not match'])
end

%  Set the option to "double-check" the backtracked trajectory step by step
if isfield(options,'adjustBT')
    adjust = options.adjustBT(1);
    maxtrials = options.adjustBT(2);
else
    adjust = 1;
    maxtrials = 10;
end

%  Getting some additional parameters controlling whether to plot the dike 
%  trajectory or not, as well as the lateral distance and depth at which
%  trajectories are cut off 
if isfield(options,'behavior')
    display = options.behavior(1);
    RigidSAM = options.behavior(2);
    displast = options.behavior(3);
else
    display = 1;
    RigidSAM = 0;
    displast = 0;
end
Traj_RadialCutoff = options.cutoff(1);
Traj_VertCutoff = options.cutoff(2);

%  Setting some parameters
condstop=0; %If == 1, the dike stops

%  Setting StartPoint as the first B. The points defining the backtracked
%  trajectory are called B in accordance with Mantiloni et al., 2023
B = StartPoint; 

%  We first perform one step of forward SAM from B, seting the magma 
%  density to Pm+Pr: this way, the dike will be anti-buoyant. Once we find
%  the direction of propagation, we will use it to initialize "ScalProdFun" 
%  later

Pm = Pr + Pm;
options.behavior = [0,RigidSAM,0]; %We do not want to display SAM pathway
options.backtracking = 1; %This way forward SAM will stop after the first step

[~,theta,al] = SAM_Forward_Analytical(B,cB,O,Pm,TopoInterp,0,Stressfun,...
    Stressfun_Par,options);

dike_path_tip = [];
dike_path_dip = [];
dike_path_strike = [];
x0 = [theta(end),al(end)]; %Initializing fminsearch in "ScalProdFun"

Pm = Pm - Pr; %Setting Pm back to its original value, so that the dike is buoyant again

%  Main loop of the function: it goes on until either the backtracked dike
%  turns back, stops horizontally, or crosses the radial or vertical 
%  threshold. If RigidSAM == 1, the trajectory will also stop if it steers,
%  dips or ascends too abruptly.

while ~condstop
    
    fun = @(x)ScalProdFun(x,B,Stressfun,Stressfun_Par,cB);
    [dike_path_angles,~,~] = fminsearch(fun,x0);
    Bc_x = B(1) + cB*cos(dike_path_angles(2))*cos(dike_path_angles(1));
    Bc_y = B(2) + cB*sin(dike_path_angles(2))*cos(dike_path_angles(1));
    Bc_z = B(3) - cB*sin(dike_path_angles(1));
    Bc = [Bc_x Bc_y Bc_z];
    
    %  Here we have found a candidate (Bc) for the current step of the 
    %  backtracked dike pathway. Now we perform forward SAM on it and see 
    %  how well the predicted forward point compares to the real, previous 
    %  one (B).
    
    %  Update the vertical threshold, otherwise forward SAM will stop 
    %  immediately
    options.cutoff(2) = Bc(3); 

    if adjust
        nn = 1;
        while nn < maxtrials
            [Bprime,~,~] = SAM_Forward_Analytical(Bc,cB,O,Pm,...
                TopoInterp,0,Stressfun,Stressfun_Par,options);
            centrediff = Bprime(end,:) - B;
            %  We shift Bc by the difference vector between the real and 
            %  predicted B
            Bc = Bc - centrediff;
            nn = nn+1;
        end
    end
    
    zcurrent = Bc(3);

    %  Determining whether Bc has hit the free surface and stop the BT 
    %  if it does
    if zcurrent > TopoInterp(Bc(1),Bc(2))
       condstop = 1;
    end
    
    %  Checking if Bc falls outside the radial threshold
    if (Bc(1)^2 + Bc(2)^2)^(1/2) > Traj_RadialCutoff
        condstop = 1;
    end
    
    %  Checking if Bc falls below the depth threshold
    if zcurrent < Traj_VertCutoff
        condstop = 1;
    end
    
    %  Check if the BT is turning back and stop if it does
    if numel(dike_path_dip) >= 2
        direction_old = B - dike_path_tip(end-1,:);
        direction_current = Bc - B;
        if dot(direction_current,direction_old) < 0
            condstop = 1;
        end
        %  If RigidSAM == 1, this checks if the BT steers, dips or ascends
        %  abruptly (fixed threshold on the change in the strike and dip
        %  angles between the current and the last trajectory step)
        if RigidSAM
            if abs(dot(direction_current,direction_old)) < 0.5
                condstop = 1;
            end
        end
        %  Check if the BT is starting to ascend (i.e., the dike has 
        %  become horizontal) and stop if it does
        if zcurrent > B(3)
            condstop = 1;
        end
    end

    if ~condstop
        x0 = dike_path_angles;
        B = Bc;
        dike_path_tip = vertcat(dike_path_tip,Bc);
        dike_path_dip = [dike_path_dip dike_path_angles(1)];
        dike_path_strike = [dike_path_strike dike_path_angles(2)];
    end
    if display && ~condstop
       scatter3(B(1),B(2),B(3),'o','g','filled')
       pause(0.1)
       hold on
    end
    if displast && condstop
       scatter3(Bc(1),Bc(2),Bc(3),'o','b','filled')
       pause(0.1)
       hold on
    end
end    

BSP = dike_path_tip(end,:);

end


function ScalProd = ScalProdFun(dike_angles,B,Stressfun,Stressfun_Par,cB)
        
dikedip = dike_angles(1);
dikestrike = dike_angles(2);
Bc_x = B(1) + cB*cos(dikestrike)*cos(dikedip);
Bc_y = B(2) + cB*sin(dikestrike)*cos(dikedip);
Bc_z = B(3) - cB*sin(dikedip);
Bc = [Bc_x Bc_y Bc_z];
    
Stressfun_Input = Input_Collecting(Bc_x,Bc_y,Bc_z,Stressfun_Par);
[~,S3dirx,S3diry,S3dirz] = Stressfun(Stressfun_Input{:});
S3dir = [S3dirx S3diry S3dirz];

RadVecMod = sqrt((B(1)-Bc(1))^2 + (B(2)-Bc(2))^2 + (B(3)-Bc(3))^2);
RadVecCompx = (B(1)-Bc(1))/RadVecMod;
RadVecCompy = (B(2)-Bc(2))/RadVecMod;
RadVecCompz = (B(3)-Bc(3))/RadVecMod;
RadVec = [RadVecCompx RadVecCompy RadVecCompz];
ScalProd = abs(dot(S3dir,RadVec));

end

function Stressfun_Input = Input_Collecting(X,Y,Z,Stressfun_Par)

    Stressfun_Input = num2cell([X Y Z]);
    Stressfun_Input = horzcat(Stressfun_Input,Stressfun_Par);
    
end

