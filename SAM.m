function[dyke_path_tip,dyke_path_dip,dyke_path_strike,varargout] = ...
    SAM(Direction,StartPoints,c,n,Pr,Pm,TopoInterp,StressModel,...
    options,varargin)

%  SAM: Simplified Analytical Model of dyke Pathways in Three Dimensions. 
%  Developed by Mantiloni et al., 2023 (full reference below). The model 
%  can either simulate forward dyke trajectories from a starting point to a
%  surface vent, or backtrack a dyke starting from a vent down to a
%  probable nucleation point.
% 
%  A dyke (both forward and backtracked) will stop by default if it hits 
%  the free surface, turns back or stops horizontally. Additional options 
%  may force dykes to stop under further conditions (see "options" in 
%  "Inputs")
%  
%  Examples are provided in the "Examples" folder.
% 
%  General conventions: everything is described in a Cartesian reference 
%  frame; the vertical coordinate is positive upward. sigma_3 or S3 denotes 
%  the least-compressive principal stress magnitude, S3x,y,z the components 
%  of the respective eigenvector.
%
%  For a detailed description of inputs, outputs and sub-routines, see the
%  Instruction Manual.
%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  Reference: 
%
%  Mantiloni, L., Rivalta, E., & Davis, T (2023). Mechanical 
%  modeling of pre-eruptive magma propagation scenarios at calderas. 
%  Journal of Geophysical Research: Solid Earth.
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
% 
%  I appreciate any feedback or bug report. Reach out to me at 
%  lorenzo@gfz-potsdam.de
%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Check if a figure is already open, and create it otherwise
gcf

%  Forward or backtrack?
if Direction ~= 'F' && Direction ~= 'B'
    error(['Please set Direction to either "F" (forward) '...
    'or "B" (backtrack).']);
end

%  Create "TopoInterp"
if ~isa(TopoInterp,'scatteredInterpolant')
    if isscalar(TopoInterp)
        TopoInterp = CreateTopoInterp(1,StartPoints);
    else
        if ~isfield(TopoInterp,'Xsurf') || ~isfield(TopoInterp,'Ysurf') ...
                || ~isfield(TopoInterp,'Zsurf')
            error(['Please provide x,y,z grids of points describing '...
                'the topography (see "Example2.m"). For a simple '...
                'flat surface, set TopoInterp = 0 (or any scalar).'])
        else
            Xsurf = TopInterp.Xsurf;
            TopoInterp = CreateTopoInterp(0,Xsurf,Ysurf,Zsurf);
        end
    end
end

%  Checking if everything is fine with "StressModel"
StressOption = StressModel.StressOption;
if ~strcmp(StressOption,'Interpolate') && ...
        ~strcmp(StressOption,'Numerical') && ...
        ~strcmp(StressOption,'Analytical')
    error(['StressOption must be set either to "Interpolate", '...
        '"Numerical" or "Analytical".'])
end
%  Checking if everything is fine with "options"
if ~isfield(options,'behavior')
    error('options must include the field "behavior" (see the Manual)')
else
    behavior = options.behavior;
    if ~isequal(size(behavior),[1,3])
        error('options.behavior must be a 1x3 array (see the Manual).')
    end
end

%  Creating the n observation points "O" (equally spaced, centered around
%  0, on a horizontal surface)
theta = pi/2;
phi = (0:2*pi/n:2*pi*(n-1)/n);
O = c*[sin(phi);cos(phi);cos(theta)*cos(phi)];

%  Preallocating the dyke pathways data as 1xN cells, N being the number of
%  dykes
N = size(StartPoints,1);
dyke_path_tip = cell(1,N); 
dyke_path_dip = cell(1,N); 
dyke_path_strike = cell(1,N);
if Direction == 'F'
    dyke_path_K = cell(1,N);
else
    BSP = cell(1,N);
end

if Direction == 'F' %Forward SAM
    
    %  Check if dykes are going to be anti-buoyant and, if they are, check
    %  if "options" includes the ".cutoff" field.
    if Pm > Pr
        if ~isfield(options,'cutoff')
            error(['Since the magma is denser than the host rocks, '...
                '"options" must include .cutoff, with maximum '...
                'distance from the origin and depth a backtracked '...
                'trajectory can reach. Otherwise, SAM may never stop'])
        end
    end

    %  Check if dyke volume is provided and, if it is, check if the
    %  fracture toughness Kc of the host rock is also provided.
    if ~isempty(varargin) 
        V = varargin{1};
        if isa(V,'double')
            V = ones(size(StartPoints,1),size(StartPoints,2))*V;
        else
            if numel(V) ~= size(StartPoints,1)
                error(['if the dykes have diffrent volumes, the size '...
                    'of V must match that of StartPoints'])
            end
        end
        if ~strcmp(StressOption,'Numerical')
            if numel(varargin) < 4
                error(['The host rock fracture toughness Kc, rigidity '...
                    'modulus mu and Poissons ratio nu are also needed '...
                    'if the dyke volume is assigned. Please provide '...
                    'them in this oder as the last arguments of SAM '...
                    'function'])
            else
                Kc = varargin{2};
                mu = varargin{3};
                nu = varargin{4};
            end
        else
            if numel(varargin) < 2
                error(['The host rock fracture toughness Kc is '...
                    'also needed if the dyke volume is assigned. '...
                    'Please provide it as the last argument of SAM '...
                    'function'])
            end
        end
    else
        V = [];
    end
    
    %  Importing the stress model.

     if strcmp(StressOption,'Analytical')
        %  Importing the analytical function handle
        Stressfun = StressModel.Stressfun;
        %  Importing the function input variables 
        Stressfun_Par = StressModel.Stressfun_Par;
        if isfield(StressModel,'MDT')
            MDT = StressModel.MDT;
        else
            MDT = 0;
        end
        if ~isempty(V)
            for i=1:N
                [dyke_path_tip_current,dyke_path_dip_current,...
                dyke_path_strike_current,dyke_path_K_current] = ...
                SAM_Forward_Analytical(StartPoints(i,:),c,O,Pm,...
                TopoInterp,MDT,Stressfun,Stressfun_Par,options,V(i),...
                Kc,mu,nu);
                dyke_path_tip{i} = dyke_path_tip_current; 
                dyke_path_dip{i} = dyke_path_dip_current;
                dyke_path_strike{i} = dyke_path_strike_current;
                dyke_path_K{i} = dyke_path_K_current;
            end
        else
            for i=1:N
                [dyke_path_tip_current,dyke_path_dip_current,...
                dyke_path_strike_current,dyke_path_K_current] = ...
                SAM_Forward_Analytical(StartPoints(i,:),c,O,Pm,...
                TopoInterp,MDT,Stressfun,Stressfun_Par,options);
                dyke_path_tip{i} = dyke_path_tip_current; 
                dyke_path_dip{i} = dyke_path_dip_current;
                dyke_path_strike{i} = dyke_path_strike_current;
                dyke_path_K{i} = dyke_path_K_current;
            end
        end
        varargout{i} = dyke_path_K;
     end

    if strcmp(StressOption,'Interpolate')
        ObsPoints = StressModel.ObsPoints;
        if ~isfield(StressModel,'InterpOption') 
            error(['In "Interp" mode, please provide the '...
                'StressModel.InterpOption field.'])
        end
        InterpOption = StressModel.InterpOption;
        if ~strcmp(InterpOption,'Cartesian') && ~strcmp(InterpOption,'S3')
            error(['InterpOption must be either "Cartesian" '...
                '(interpolating cartesian stress tensor components '...
                'Sij) or "S3" (interpolating S3 and S3dir directly).'])
        end
        if strcmp(InterpOption,'Cartesian')
            Stensor = StressModel.Stensor;
            %  Interpolating Sij components over the observation points 
            %  grid
            [SxxInterp,SyyInterp,SzzInterp,SxyInterp,SxzInterp,...
                SyzInterp] = StressInterp(InterpOption,ObsPoints,Stensor);
        else
            S3 = StressModel.S3;
            S3dir = StressModel.S3dir;
            if ~isfield(options,'adjustdir')
                AdjustDir = logical(options.adjustdir);
            else
                AdjustDir = 1;
            end
             %  Interpolating S3 and its components over the observation
             %  points grid
            [S3Interp,S3xInterp,S3yInterp,S3zInterp] = ...
                StressInterp(InterpOption,ObsPoints,S3,S3dir,AdjustDir);
        end
        if isfield(StressModel,'MDT')
            MDT = StressModel.MDT;
        else
            MDT = 0;
        end
        for i=1:N
            if ~isempty(V)
                if strcmp(InterpOption,'Cartesian')
                    [dyke_path_tip_current,dyke_path_dip_current,...
                    dyke_path_strike_current,dyke_path_K_current] = ...
                    SAM_Forward_Interp(StartPoints(i,:),c,O,Pm,...
                    TopoInterp,MDT,options,InterpOption,SxxInterp,...
                    SyyInterp,SzzInterp,SxyInterp,SxzInterp,...
                    SyzInterp,V(i),Kc,mu,nu);
                else
                    [dyke_path_tip_current,dyke_path_dip_current,...
                    dyke_path_strike_current,dyke_path_K_current] = ...
                    SAM_Forward_Interp(StartPoints(i,:),c,O,Pm,...
                    TopoInterp,MDT,options,InterpOption,S3Interp,...
                    S3xInterp,S3yInterp,S3zInterp,V(i),Kc,mu,nu);
                end
                dyke_path_tip{i} = dyke_path_tip_current; 
                dyke_path_dip{i} = dyke_path_dip_current;
                dyke_path_strike{i} = dyke_path_strike_current;
                dyke_path_K{i} = dyke_path_K_current;
            else
                if strcmp(InterpOption,'Cartesian')
                    [dyke_path_tip_current,dyke_path_dip_current,...
                    dyke_path_strike_current,dyke_path_K_current] = ...
                    SAM_Forward_Interp(StartPoints(i,:),c,O,Pm,...
                    TopoInterp,MDT,options,InterpOption,SxxInterp,...
                    SyyInterp,SzzInterp,SxyInterp,SxzInterp,SyzInterp);
                else
                    [dyke_path_tip_current,dyke_path_dip_current,...
                    dyke_path_strike_current,dyke_path_K_current] = ...
                    SAM_Forward_Interp(StartPoints(i,:),c,O,Pm,...
                    TopoInterp,MDT,options,InterpOption,S3Interp,...
                    S3xInterp,S3yInterp,S3zInterp);
                end
                dyke_path_tip{i} = dyke_path_tip_current; 
                dyke_path_dip{i} = dyke_path_dip_current;
                dyke_path_strike{i} = dyke_path_strike_current;
                dyke_path_K{i} = dyke_path_K_current;
            end
        end 
        varargout{1} = dyke_path_K;
    end

    if strcmp(StressOption,'Numerical')
        %  Importing elastic parameters and density of host rock
        RockPar = StressModel.RockPar;
        %  Importing the vertices of triangular dislocations of mesh
        P1 = StressModel.P1;
        P2 = StressModel.P2;
        P3 = StressModel.P3; 
        DisplMat = StressModel.DisplMat;
        %  Getting the displacement components of triangular dislocations 
        %  from the Boundary Elements solution
        Dds = DisplMat(:,1);
        Dss = DisplMat(:,2);
        Dn = DisplMat(:,3);
        %  Importing tectonic stress components
        TectS = StressModel.TectS;
        if isfield(StressModel,'MDT')
            MDT = StressModel.MDT;
        else
            MDT = findMDT(P1,P2,P3);
        end
        if ~isempty(V)
            for i=1:N
                [dyke_path_tip_current,dyke_path_dip_current,...
                dyke_path_strike_current,dyke_path_K_current] = ...
                SAM_Forward_Numerical(StartPoints(i,:),c,O,RockPar,Pr,...
                Pm,TopoInterp,MDT,Dss,Dds,Dn,P1,P2,P3,TectS,options,...
                V(i),Kc);
                dyke_path_tip{i} = dyke_path_tip_current; 
                dyke_path_dip{i} = dyke_path_dip_current;
                dyke_path_strike{i} = dyke_path_strike_current;
                dyke_path_K{i} = dyke_path_K_current;
            end
        else
            for i=1:N
                [dyke_path_tip_current,dyke_path_dip_current,...
                dyke_path_strike_current,dyke_path_K_current] = ...
                SAM_Forward_Numerical(StartPoints(i,:),c,O,RockPar,Pr,...
                Pm,TopoInterp,MDT,Dss,Dds,Dn,P1,P2,P3,TectS,options);
                dyke_path_tip{i} = dyke_path_tip_current; 
                dyke_path_dip{i} = dyke_path_dip_current;
                dyke_path_strike{i} = dyke_path_strike_current;
                dyke_path_K{i} = dyke_path_K_current;
            end
        end
        varargout{1} = dyke_path_K;
    end

end

if Direction == 'B' %Backtrack SAM
    
    %  Importing the stress model.

    if strcmp(StressOption,'Analytical')
        %  Importing the analytical function handle
        Stressfun = StressModel.Stressfun;
        %  Importing the function input variables 
        Stressfun_Par = StressModel.Stressfun_Par;
        for i=1:N
            [dyke_path_tip_current,dyke_path_dip_current,...
                dyke_path_strike_current,BSP_current] = ...
                SAM_Backtrack_Analytical(StartPoints(i,:),c,O,Pr,Pm,...
                TopoInterp,Stressfun,Stressfun_Par,options,varargin);
            dyke_path_tip{i} = dyke_path_tip_current; 
            dyke_path_dip{i} = dyke_path_dip_current;
            dyke_path_strike{i} = dyke_path_strike_current;
            BSP{i} = BSP_current;
        end 
        varargout{1} = BSP;
    end

    if strcmp(StressOption,'Interpolate')
        ObsPoints = StressModel.ObsPoints;
        InterpOption = StressModel.InterpOption;
        if ~strcmp(InterpOption,'Cartesian') && ~strcmp(InterpOption,'S3')
            error(['InterpOption must be either "Cartesian" '...
                '(interpolating cartesian stress tensor components '...
                'Sij) or "S3" (interpolating S3 and S3dir directly).'])
        end
        if strcmp(InterpOption,'Cartesian')
            Stensor = StressModel.Stensor;
            %  Interpolating Sij components over the observation points
            %  grid
            [SxxInterp,SyyInterp,SzzInterp,SxyInterp,SxzInterp,...
                SyzInterp] = StressInterp(InterpOption,ObsPoints,Stensor);
        else
            S3 = StressModel.S3;
            S3dir = StressModel.S3dir;
            if ~isfield(options,'adjustdir')
                AdjustDir = logical(options.adjustdir);
            else
                AdjustDir = 1;
            end
             %  Interpolating S3 and its components over the observation 
             %  points grid
            [S3Interp,S3xInterp,S3yInterp,S3zInterp] = ...
                StressInterp(InterpOption,ObsPoints,S3,S3dir,AdjustDir);
        end
        if isfield(options,'FromFreeSurface')   
            if ~isfield(options,'ProjectToFreeSurface')
                error(['When backtracking dykes from the free surface, '...
                    '"option.ToFreeSurface" is also required as it '...
                    'would be in forward mode.'])
            end
            if isfield(StressModel,'MDT')
                MDT = StressModel.MDT;
                varargin{1} = MDT;
            else
                error(['When backtracking dykes from the free surface '...
                    'with a numerical stress model, StressModel.MDT '...
                    'is required'])
            end
        end
        if ~isfield(options,'cutoff')
            error(['"options" must include .cutoff, with maximum '...
                'distance from the origin and depth a backtracked '...
                'trajectory can reach. Otherwise, SAM will never stop'])
        end
        for i=1:N
            if strcmp(InterpOption,'Cartesian')
                [dyke_path_tip_current,dyke_path_dip_current,...
                    dyke_path_strike_current,BSP_current] = ...
                    SAM_Backtrack_Interp(StartPoints(i,:),c,O,Pr,Pm,...
                    TopoInterp,options,InterpOption,SxxInterp,...
                    SyyInterp,SzzInterp,SxyInterp,SxzInterp,...
                    SyzInterp,varargin);
            else
                [dyke_path_tip_current,dyke_path_dip_current,...
                    dyke_path_strike_current,BSP_current] = ...
                    SAM_Backtrack_Interp(StartPoints(i,:),c,O,Pr,Pm,...
                    TopoInterp,options,InterpOption,S3Interp,...
                    S3xInterp,S3yInterp,S3zInterp,varargin);
            end
            dyke_path_tip{i} = dyke_path_tip_current; 
            dyke_path_dip{i} = dyke_path_dip_current;
            dyke_path_strike{i} = dyke_path_strike_current;
            BSP{i} = BSP_current;
        end 
        varargout{1} = BSP;
    end

    if strcmp(StressOption,'Numerical')
        %  Importing elastic parameters and density of host rock
        RockPar = StressModel.RockPar;
        %  Importing the vertices of triangular dislocations of mesh
        P1 = StressModel.P1;
        P2 = StressModel.P2;
        P3 = StressModel.P3; 
        DisplMat = StressModel.DisplMat;
        %  Getting the displacement components of triangular dislocations 
        %  from the Boundary Elements solution
        Dds = DisplMat(:,1);
        Dss = DisplMat(:,2);
        Dn = DisplMat(:,3);
        %  Importing tectonic stress components
        TectS = StressModel.TectS;
        if isfield(options,'FromFreeSurface')    
            if ~isfield(options,'ToFreeSurface')
                error(['When backtracking dykes from the free surface, '...
                    '"option.ToFreeSurface" is also required as it '...
                    'would be in forward mode.'])
            end
            if isfield(StressModel,'MDT')
                MDT = StressModel.MDT;
                varargin{1} = MDT;
            else
                error(['When backtracking dykes from the free surface '...
                    'with a numerical stress model, StressModel.MDT '...
                    'is required'])
            end
        end
        for i=1:N
            [dyke_path_tip_current,dyke_path_dip_current,...
                dyke_path_strike_current,BSP_current] = ...
                SAM_Backtrack_Numerical(StartPoints(i,:),c,O,RockPar,...
                Pm,TopoInterp,Dss,Dds,Dn,P1,P2,P3,TectS,options,varargin);
            dyke_path_tip{i} = dyke_path_tip_current; 
            dyke_path_dip{i} = dyke_path_dip_current;
            dyke_path_strike{i} = dyke_path_strike_current;
            BSP{i} = BSP_current;
        end 
        varargout{1} = BSP;
    end
    
end

end

function [varargout] = StressInterp(InterpOption,ObsPoints,varargin)

% Extracting the observation points 
Xstressgrid = ObsPoints(:,1); 
Ystressgrid = ObsPoints(:,2); 
Zstressgrid = ObsPoints(:,3);

if strcmp(InterpOption,'Cartesian')

    Stensor = varargin{1};
    %  Extracting stress components
    Sxx = Stensor(:,1); Syy = Stensor(:,2); Szz = Stensor(:,3);
    Sxy = Stensor(:,4); Sxz = Stensor(:,5); Syz = Stensor(:,6);

    SxxInterp = scatteredInterpolant(Xstressgrid,Ystressgrid,...
        Zstressgrid,Sxx);
    SyyInterp = scatteredInterpolant(Xstressgrid,Ystressgrid,...
        Zstressgrid,Syy);
    SzzInterp = scatteredInterpolant(Xstressgrid,Ystressgrid,...
        Zstressgrid,Szz);
    SxyInterp = scatteredInterpolant(Xstressgrid,Ystressgrid,...
        Zstressgrid,Sxy);
    SxzInterp = scatteredInterpolant(Xstressgrid,Ystressgrid,...
        Zstressgrid,Sxz);
    SyzInterp = scatteredInterpolant(Xstressgrid,Ystressgrid,...
        Zstressgrid,Syz);

    varargout{1} = SxxInterp; 
    varargout{2} = SyyInterp; 
    varargout{3} = SzzInterp;
    varargout{4} = SxyInterp; 
    varargout{5} = SxzInterp; 
    varargout{6} = SyzInterp;

end

if strcmp(InterpOption,'S3')

    S3 = varargin{1};
    S3dir = varargin{2};
    AdjustDir = varargin{3};

    %  Adjusting the provided S3 components so that it points upward 
    %  everywhere (not doing so would result in artifacts in the 
    %  interpolation)
    if AdjustDir
        benchmarkvectorZ = S3dir*0 + [0 0 1]; %Generic vector parallel to z-axis and pointing upward
        benchmarkvectorZ = benchmarkvectorZ'; 
        condflipZ = dot(S3dir',benchmarkvectorZ) < 0;
        condflipZ = condflipZ';
        S3dir(condflipZ,1) = -S3dir(condflipZ,1);
        S3dir(condflipZ,2) = -S3dir(condflipZ,2);
        S3dir(condflipZ,3) = -S3dir(condflipZ,3);
        benchmarkvectorX = S3dir*0 + [1 0 0]; %Generic vector parallel to x-axis and pointing to the positive direction
        benchmarkvectorX = benchmarkvectorX'; 
        condflipX = dot(S3dir',benchmarkvectorX) < 0;
        condflipX = condflipX';
        S3dir(condflipX,1) = -S3dir(condflipX,1);
        S3dir(condflipX,2) = -S3dir(condflipX,2);
        S3dir(condflipX,3) = -S3dir(condflipX,3);
        benchmarkvectorY = S3dir*0 + [0 1 0]; %Generic vector parallel to y-axis and pointing to the positive direction
        benchmarkvectorY = benchmarkvectorY'; 
        condflipY = dot(S3dir',benchmarkvectorY) < 0;
        condflipY = condflipY';
        S3dir(condflipY,1) = -S3dir(condflipY,1);
        S3dir(condflipY,2) = -S3dir(condflipY,2);
        S3dir(condflipY,3) = -S3dir(condflipY,3);
    end

    S3Interp = scatteredInterpolant(Xstressgrid,Ystressgrid,...
        Zstressgrid,S3);
    S3xInterp = scatteredInterpolant(Xstressgrid,Ystressgrid,...
        Zstressgrid,S3dir(:,1));
    S3yInterp = scatteredInterpolant(Xstressgrid,Ystressgrid,...
        Zstressgrid,S3dir(:,2));
    S3zInterp = scatteredInterpolant(Xstressgrid,Ystressgrid,...
        Zstressgrid,S3dir(:,3));

    varargout{1} = S3Interp;
    varargout{2} = S3xInterp;
    varargout{3} = S3yInterp;
    varargout{4} = S3zInterp;

end

end

function TopoInterp = CreateTopoInterp(Flat,varargin)

if Flat 

    StartPoints = varargin{1};
    minStartPoint = min([min(StartPoints(:,1)) min(StartPoints(:,2))]);
    maxStartPoint = max([max(StartPoints(:,1)) max(StartPoints(:,2))]);
    Xsurf = linspace(minStartPoint,maxStartPoint,10);
    Ysurf = Xsurf;
    [Xsurf,Ysurf] = meshgrid(Xsurf,Ysurf);
    Zsurf = Xsurf*0;

else

    Xsurf = varargin{1};
    Ysurf = varargin{2};
    Zsurf = varargin{3};

end

TopoInterp = scatteredInterpolant(Xsurf(:),Ysurf(:),Zsurf(:));

end

function MDT = findMDT(P1,P2,P3)

minP12dist = min(((P1(:,1)-P2(:,1)).^2 + (P1(:,2) - P2(:,2)).^2 + (P1(:,3)-P2(:,3)).^2).^(1/2));
minP23dist = min(((P2(:,1)-P3(:,1)).^2 + (P2(:,2) - P3(:,2)).^2 + (P2(:,3)-P3(:,3)).^2).^(1/2));
minP31dist = min(((P3(:,1)-P1(:,1)).^2 + (P3(:,2) - P1(:,2)).^2 + (P3(:,3)-P1(:,3)).^2).^(1/2));

MDT = min([minP12dist,minP23dist,minP31dist]);

end