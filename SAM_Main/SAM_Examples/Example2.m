%  Calculates and plots the pathways of a set of SAM dikes starting 
%  from the edge of a sill-like reservoir below a caldera lying on a 
%  coastline and bounded by hills (see Complex-Coastline scenario from
%  Mantiloni et al., 2023). The stress field, due to the gravitational
%  loading/unloading associated with the topography and a superimposed 
%  tectonic stress, is imported from a .mat file (Example2_Data.mat).
%  SAM runs in "Interpolate" mode.
%
%  Colour conventions: SAM forward pathways are represented by red patches,
%  and the final step is marked by a yellow patch. Arrival points projected 
%  to the free surface are shown as by blue dots. SAM backtracked 
%  trajectories are represented by green and cyan dots when starting from
%  the true arrival point and from its projection to the free surface,
%  respectively.The points where they stop are also marked by blue dots.
% 
%  I appreciate any feedback or bug report. Reach out to me at 
%  lorenzo@gfz-potsdam.de
%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  References:
%
%  Mantiloni, L., Rivalta, E., & Davis, T. (2023). Mechanical Modeling of 
%  Pre-Eruptive Magma Propagation Scenarios at Calderas. Journal of 
%  Geophysical Research: Solid Earth, 128, 3. 
%  https://doi.org/10.1029/2022JB025956

%  Import the stress model, topography mesh and scenario parameters:
load 'Example2_Data.mat'

%  Variables included in "Example2_Data.mat":
%  _ ObsPoints: P x 6 array with x, y, z coordinates of observation points
%  _ Points, Triangles: together they describe the mesh of triangular
%    dislocations
%  _ StartPoints: x, y, z coordinates of 10 dike starting points
%  _ Sxx, Sxy, Sxz, Syy, Syz, Szz: components of the Cartesian stress
%    tensor evaluated at each observation point
%  _ TectS: Sxx, Syy, Sxy components of tectonic stress
%  _ Topo_Interp: interpolating function for the topography giving the
%    height of free surface at point (x, y)

%  Host rock parameters
Pr = 2850; %Host rock density (kg/m^3)

%  Setting dike starting points and parameters 
ImportStartPoint = 0; %Set to 0 to use the imported starting points
if ~ImportStartPoint
    N = 1; %Number of dikes
    th=0:2*pi/N:2*pi-2*pi/N;
    rad=3e3;
    StartX = -rad * cos(th);
    StartY = rad * sin(th);
    StartZ = StartX*0 - 6e3;
    StartPoints = [StartX' StartY' StartZ'];
else
    N = size(StartPoints,1);
end
c = 800;
n = 12;
Pm = 2300;

%  Setting Direction
Direction = 'F';

%  Creating the Stensor variable 
Stensor = [Sxx Syy Szz Sxy Sxz Syz];

%  Setting StressModel
StressModel.StressOption = 'Interpolate'; %Cambiare in Interpolate
StressModel.InterpOption = 'Cartesian';
StressModel.ObsPoints = ObsPoints;
StressModel.Stensor = Stensor;
StressModel.MDT = 800;

%  Setting options
%  Setting plotting and propagation/arrest options
display = 1; %Display SAM pathway
RigidSAM = 1; %Preventing SAM dikes from dipping/steering abruptly
displast = 1; %Display the last step of SAM pathways in yellow
options.behavior = [display,RigidSAM,displast]; 
options.adjustdir = 1;
%  Project SAM pathways to the free surface 
PathStep = 60;
options.ProjectToFreeSurface = [1,PathStep]; 

%  Plotting the boundary-element mesh representing the topography
figure
trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4)*5,'FaceAlpha',...
    (.2),'FaceColor', [0.5 0 0.9 ]);
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); axis equal
hold on

%  Running SAM:
[dike_path_tip,dike_path_dip,~] = SAM(Direction,StartPoints,c,n,Pr,Pm,...
    Topo_Interp,StressModel,options);

%  Retrieving the arrival poins for each SAM pathway
ArrivalPoints_ProjectedToFreeSurface = zeros(N,3);
ArrivalPoints_Actual = zeros(N,3);
for i=1:numel(dike_path_tip)
    dike = dike_path_tip{i};
    PathwayDifference = size(dike,1) - size(dike_path_dip{i},2);
    ArrivalPoints_ProjectedToFreeSurface(i,:) = dike(end,:);
    ArrivalPoints_Actual(i,:) = dike(end-PathwayDifference,:);
end

%  Backtracking each pathway from the respective arrival point 
Direction = 'B';
%  Here we employ the same radius c as in the forward model, but they may 
%  be different
cB = c; 
%  We stop SAM backtracked trajectories at a depth slightly larger than the 
%  starting depth of dikes
options.cutoff = [1e4,-6.1e3];

%  First backtracking SAM dikes from their actual arrival points
[dike_path_tip1,dike_path_dip1,dike_path_strike1,BSP1] = SAM(Direction,...
    ArrivalPoints_Actual,cB,n,Pr,Pm,Topo_Interp,StressModel,options,Pr);

%  Extract BSPs
BSPs = zeros(N,3);
for i=1:N
    BSPs(i,:) = BSP1{i};
end

%  Calculate the difference between BSPs and actual starting points
DeltaBSP1 = abs(BSPs - StartPoints);
meanDeltaBSP1 = (sum((DeltaBSP1(:,1).^2 + DeltaBSP1(:,2).^2 + ...
    DeltaBSP1(:,3).^2).^(1/2)))/N;
message = ['Starting points retrieved with a mean absolute error of ' ...
    num2str(meanDeltaBSP1) ' m'];
disp(message)

%  Then repeating the procedure with projected arrival points
Method = 1;
options.FromFreeSurface = [1 Method];
[dike_path_tip2,dike_path_dip2,dike_path_strike2,BSP2] = SAM(Direction,...
    ArrivalPoints_ProjectedToFreeSurface,cB,n,Pr,Pm,Topo_Interp,...
    StressModel,options,Pr);

%  Extract BSPs
BSPs = zeros(N,3);
for i=1:N
    BSPs(i,:) = BSP2{i};
end

%  Calculate the difference between BSPs and actual starting points
DeltaBSP2 = abs(BSPs - StartPoints);
meanDeltaBSP2 = (sum((DeltaBSP2(:,1).^2 + DeltaBSP2(:,2).^2 + ...
    DeltaBSP2(:,3).^2).^(1/2)))/N;
message = ['Starting points retrieved with a mean absolute error of ' ...
    num2str(meanDeltaBSP2) ' m'];
disp(message)

%%

%  Running another forward simulation with the same starting points, but 
%  also setting a common dike volume

V = 1e5; %  Dike volume (common, m^3)

%  Host rock parameters
Kc = 70e6; %  Rock fracture toughness (Pa*m^(1/2)
mu = 6e9; % Rigidity modulus (Pa)
nu = 0.25; % Poisson's ratio

%  Setting Direction
Direction = 'F';

figure
trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4)*5,'FaceAlpha',...
    (.2),'FaceColor', [0.5 0 0.9 ]);
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); axis equal
hold on

[dike_path_tip,dike_path_dip,~] = SAM(Direction,StartPoints,c,n,Pr,Pm,...
    Topo_Interp,StressModel,options,V,Kc,mu,nu);
