function[dike_path_tip,dike_path_dip,dike_path_strike,BSP] = ...
    SAM_Backtrack_Interp(StartPoint,cB,O,Pr,Pm,TopoInterp,options,...
    InterpOption,varargin)

%  "SAM_Backtrack_Interp" backtracks dike trajectories from starting 
%  points, either on the surface or in the subsurface, down to a potential
%  magma storage volume. Stresses are evaluated through interpolating 
%  functions for either the individual components of the cartesian stress
%  tensor, or the magnitude and direction components of sigma_3. For 
%  further details on inputs and outputs, see the Instruction Manual, as 
%  well as the inline comments.
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

%  Set the option to "double-check" the backtracked trajectory step by step
if isfield(options,'adjustBT')
    adjust = options.adjustBT(1);
    maxtrials = options.adjustBT(2);
else
    adjust = 1;
    maxtrials = 10;
end

%  Importing the interpolating functions for the external stress
if strcmp(InterpOption,'Cartesian')
    Cartesian = 1;
    SxxInterp = varargin{1}; SyyInterp = varargin{2}; SzzInterp = varargin{3};
    SxyInterp = varargin{4}; SxzInterp = varargin{5}; SyzInterp = varargin{6};
else
    Cartesian = 0;
    S3Interp = varargin{1}; S3xInterp = varargin{2}; 
    S3yInterp = varargin{3}; S3zInterp = varargin{4};
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
%  Are we going to backtrack dikes from projected or actual arrival points?
if isfield(options,'FromFreeSurface') 
    FromFreeSurface = options.FromFreeSurface(1);
    Method = options.FromFreeSurface(2);
    if Method ~= 1 && Method ~= 2
        error(['options.FromFreeSurface(2) must be set to either 1 or '...
            '2 (see the Instruction Manual).'])
    end
    if Cartesian
        MDT = varargin{7};
    else
        MDT = varargin{5};
    end
    MDT = cell2mat(MDT);
    if numel(options.ProjectToFreeSurface) > 1
        PathStep = options.ProjectToFreeSurface(2);
    else
        PathStep = 50;
    end
else
    FromFreeSurface = 0;
end
Traj_RadialCutoff = options.cutoff(1);
Traj_VertCutoff = options.cutoff(2);

%  Setting some parameters
condstop=0; %If == 1, the dike stops

if FromFreeSurface
    color = 'c'; %Setting the color of backtracked trajectory points
    if Cartesian
        B = FirstBacktrackPointFinder(StartPoint,TopoInterp,MDT,...
            maxtrials,PathStep,Cartesian,display,Method,SxxInterp,...
            SyyInterp,SzzInterp,SxyInterp,SxzInterp,SyzInterp);
    else
        B = FirstBacktrackPointFinder(StartPoint,TopoInterp,MDT,...
            maxtrials,PathStep,Cartesian,display,Method,S3xInterp,...
            S3yInterp,S3zInterp);
    end
else
    color = 'g';
    %  Setting StartPoint as the first B. The points defining the
    %  backtracked trajectory are called B in accordance with 
    %  Mantiloni et al., 2023
    B = StartPoint; 
end

%We first perform one step of forward SAM from B, seting the magma density
%to Pm+Pr: this way, the dike will be anti-buoyant. Once we find the
%direction of propagation, we will use it to initialize "ScalProdFun" later

Pm = Pr + Pm;
options.behavior = [0,RigidSAM,0]; %We do not want to display SAM pathway
options.backtracking = 1; %This way forward SAM will stop after the first step
if isfield(options,'ProjectToFreeSurface')
    options.ProjectToFreeSurface(1) = 0; %We don't want forward SAM to project the arrival point to the free surface
end

if Cartesian
    [~,theta,al] = SAM_Forward_Interp(B,cB,O,Pm,TopoInterp,0,options,...
        InterpOption,SxxInterp,SyyInterp,SzzInterp,SxyInterp,...
        SxzInterp,SyzInterp);
else
    [~,theta,al] = SAM_Forward_Interp(B,cB,O,Pm,TopoInterp,0,options,...
        InterpOption,S3Interp,S3xInterp,S3yInterp,S3zInterp);
end

dike_path_tip = [];
dike_path_dip = [];
dike_path_strike = [];
x0 = [theta,al]; %Initializing fminsearch in "ScalProdFun"

Pm = Pm - Pr; %Setting Pm back to its original value, so that the dike is buoyant again

%  Main loop of the function: it goes on until either the backtracked dike
%  turns back, stops horizontally, or crosses the radial or vertical 
%  threshold. If RigidSAM == 1, the trajectory will also stop if it steers,
%  dips or ascends too abruptly.

while ~condstop
    
    if Cartesian
        fun = @(x)ScalProdFun(x,B,cB,Cartesian,SxxInterp,SyyInterp,...
            SzzInterp,SxyInterp,SxzInterp,SyzInterp);
    else
        fun = @(x)ScalProdFun(x,B,cB,Cartesian,S3Interp,S3xInterp,...
            S3yInterp,S3zInterp);
    end
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
            if Cartesian
                [Bprime,~,~] = SAM_Forward_Interp(Bc,cB,O,Pm,TopoInterp,...
                    0,options,InterpOption,SxxInterp,SyyInterp,...
                    SzzInterp,SxyInterp,SxzInterp,SyzInterp);
            else
                [Bprime,~,~] = SAM_Forward_Interp(Bc,cB,O,Pm,TopoInterp,...
                    0,options,InterpOption,S3Interp,S3xInterp,...
                    S3yInterp,S3zInterp);
            end
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
            if abs(dot(direction_current,direction_old)) < 0.75
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
       scatter3(B(1),B(2),B(3),'o',color,'filled')
       pause(0.1)
       hold on
    end
    if displast && condstop
       scatter3(Bc(1),Bc(2),Bc(3),'o','b','filled')
       pause(0.1)
       hold on
    end
end    

if isempty(dike_path_tip)
    disp('Warning: the backtracked trajectory was stopped immediately')
    BSP = [];
else
    BSP = dike_path_tip(end,:);
end

end


function ScalProd = ScalProdFun(dike_angles,B,cB,Cartesian,varargin)
        
dikedip = dike_angles(1);
dikestrike = dike_angles(2);
Bc_x = B(1) + cB*cos(dikestrike)*cos(dikedip);
Bc_y = B(2) + cB*sin(dikestrike)*cos(dikedip);
Bc_z = B(3) - cB*sin(dikedip);
Bc = [Bc_x Bc_y Bc_z];
    
%  Importing the interpolating functions for the external stress
if Cartesian
    SxxInterp = varargin{1}; 
    SyyInterp = varargin{2}; 
    SzzInterp = varargin{3};
    SxyInterp = varargin{4}; 
    SxzInterp = varargin{5}; 
    SyzInterp = varargin{6};
else
    S3xInterp = varargin{1}; 
    S3yInterp = varargin{2}; 
    S3zInterp = varargin{3};
end

if Cartesian
     [~,S3dir] = Retrieve_sigma3_Interp(Bc_x,Bc_y,Bc_z,SxxInterp,...
         SyyInterp,SzzInterp,SxyInterp,SxzInterp,SyzInterp);
else
    S3dir = [S3xInterp(Bc_x,Bc_y,Bc_z) S3yInterp(Bc_x,Bc_y,Bc_z) ...
        S3zInterp(Bc_x,Bc_y,Bc_z)];
end

RadVecMod = sqrt((B(1)-Bc(1))^2 + (B(2)-Bc(2))^2 + (B(3)-Bc(3))^2);
RadVecCompx = (B(1)-Bc(1))/RadVecMod;
RadVecCompy = (B(2)-Bc(2))/RadVecMod;
RadVecCompz = (B(3)-Bc(3))/RadVecMod;
RadVec = [RadVecCompx RadVecCompy RadVecCompz];
ScalProd = abs(dot(S3dir,RadVec));

end

function ScalProd = ScalProdFun3FreePar(dikePrevStep,B,Cartesian,...
    MDT,varargin)
        
cB = dikePrevStep(1);
dikedip = dikePrevStep(2);
dikestrike = dikePrevStep(3);
Bc_x = B(1) + cB*cos(dikestrike)*cos(dikedip);
Bc_y = B(2) + cB*sin(dikestrike)*cos(dikedip);
Bc_z = B(3) - cB*sin(dikedip);
Bc = [Bc_x Bc_y Bc_z];
    
%  Importing the interpolating functions for the external stress
if Cartesian
    SxxInterp = varargin{1}; 
    SyyInterp = varargin{2}; 
    SzzInterp = varargin{3};
    SxyInterp = varargin{4}; 
    SxzInterp = varargin{5}; 
    SyzInterp = varargin{6};
else
    S3xInterp = varargin{1}; 
    S3yInterp = varargin{2}; 
    S3zInterp = varargin{3};
end

if Cartesian
     [~,S3dir] = Retrieve_sigma3_Interp(Bc_x,Bc_y,Bc_z,SxxInterp,...
         SyyInterp,SzzInterp,SxyInterp,SxzInterp,SyzInterp);
else
    S3dir = [S3xInterp(Bc_x,Bc_y,Bc_z) S3yInterp(Bc_x,Bc_y,Bc_z) ...
        S3zInterp(Bc_x,Bc_y,Bc_z)];
end

RadVecMod = sqrt((B(1)-Bc(1))^2 + (B(2)-Bc(2))^2 + (B(3)-Bc(3))^2);
RadVecCompx = (B(1)-Bc(1))/RadVecMod;
RadVecCompy = (B(2)-Bc(2))/RadVecMod;
RadVecCompz = (B(3)-Bc(3))/RadVecMod;
RadVec = [RadVecCompx RadVecCompy RadVecCompz];
if Bc_z > -0.5*MDT
    ScalProd = 1;
else
    ScalProd = abs(dot(S3dir,RadVec));
end

end

function Bc = FirstBacktrackPointFinder(StartPoint,TopoInterp,MDT,...
    maxtrials,PathStep,Cartesian,display,Method,varargin)

%  Importing the interpolating functions for the external stress
if Cartesian
    SxxInterp = varargin{1}; 
    SyyInterp = varargin{2}; 
    SzzInterp = varargin{3};
    SxyInterp = varargin{4}; 
    SxzInterp = varargin{5}; 
    SyzInterp = varargin{6};
else
    S3xInterp = varargin{1}; 
    S3yInterp = varargin{2}; 
    S3zInterp = varargin{3};
end

diffthr = 65;

Candidates = zeros(3,4);
    
for k=1:3

    switch Method

        case 1

            Bc = [StartPoint(1),StartPoint(2),StartPoint(3)-k*MDT];

            nn = 1;
           
            while nn < maxtrials
        
                if Cartesian
                    [~,S3dir] = Retrieve_sigma3_Interp(Bc(1),Bc(2),...
                        Bc(3),SxxInterp,SyyInterp,SzzInterp,SxyInterp,...
                        SxzInterp,SyzInterp);
                else
                    S3dir = [S3xInterp(Bc(1),Bc(2),Bc(3)) ...
                        S3yInterp(Bc(1),Bc(2),Bc(3)) ...
                        S3zInterp(Bc(1),Bc(2),Bc(3))];
                end
                S3dirmod = (S3dir(1)^2 + S3dir(2)^2 + S3dir(3)^2)^(1/2);
                theta = acos(S3dir(3)/S3dirmod);
                al = alpha_finder(S3dir,theta);
                Bc_start = Bc;
                Bc_orig = Bc;

                while Bc(3) < TopoInterp(Bc(1),Bc(2)) 
                    Bc(1) = Bc(1) - PathStep*cos(theta)*cos(al);
                    Bc(2) = Bc(2) - PathStep*cos(theta)*sin(al);
                    Bc(3) = Bc(3) + PathStep*sin(theta);
                    Bc_orig = vertcat(Bc_orig,Bc);
                end

                if size(Bc_orig,1) > 1
                    Bc = Bc_orig(end-1,:);
                end

                SPdiff = Bc - StartPoint;
                modSPdiff = sqrt(SPdiff(1)^2 + SPdiff(2)^2 + SPdiff(3)^2);
                if modSPdiff < diffthr
                    nn = maxtrials;
                else
                    nn = nn+1;
                end
                Bc = Bc_start - SPdiff;

            end

        case 2

            al_0 = atan2(StartPoint(2),StartPoint(1));
            x0 = [k*MDT,pi/2,al_0];

            nn = 1;

            while nn < maxtrials

                if Cartesian
                    fun = @(x)ScalProdFun3FreePar(x,StartPoint,...
                        Cartesian,MDT,SxxInterp,SyyInterp,SzzInterp,...
                        SxyInterp,SxzInterp,SyzInterp);
                else
                    fun = @(x)ScalProdFun3FreePar(x,StartPoint,...
                        Cartesian,MDT,S3Interp,S3xInterp,S3yInterp,...
                        S3zInterp);
                end
                [BcPar] = fminsearch(fun,x0);
                Bc = [StartPoint(1) + BcPar(1)*cos(BcPar(3))*...
                    cos(BcPar(2)),StartPoint(2) + BcPar(1)*...
                    sin(BcPar(3))*cos(BcPar(2)),StartPoint(3)...
                    - BcPar(1)*sin(BcPar(2))];

                nn = nn+1;
                x0 = BcPar;

                Bc_start = Bc;
                Bc_orig = Bc;

                while Bc(3) < TopoInterp(Bc(1),Bc(2)) 
                    Bc(1) = Bc(1) - PathStep*cos(BcPar(2))*cos(BcPar(3));
                    Bc(2) = Bc(2) - PathStep*cos(BcPar(2))*sin(BcPar(3));
                    Bc(3) = Bc(3) + PathStep*sin(BcPar(2));
                    Bc_orig = vertcat(Bc_orig,Bc);
                    if Bc(3) < Bc_orig(end-1,3)
                        Bc(3) = TopoInterp(Bc(1),Bc(2)) + 1;
                    end
                end

                if size(Bc_orig,1) > 1
                    Bc = Bc_orig(end-1,:);
                end

                SPdiff = Bc - StartPoint;
                modSPdiff = sqrt(SPdiff(1)^2 + SPdiff(2)^2 + SPdiff(3)^2);
                Bc = Bc_start - SPdiff;

            end

    end

    Candidates(k,:) = [Bc modSPdiff];

end

condAboveMDT = Candidates(:,3) > -MDT;
Candidates = Candidates(~condAboveMDT,:);
Bc_best = find(Candidates(:,4) == min(Candidates(:,4)));
Bc = Candidates(Bc_best,1:3);

if display
    scatter3(Bc(1),Bc(2),Bc(3),'o','c','filled')
    pause(0.1)
    hold on
end

end


