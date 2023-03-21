%  Calculates and plots the pathways of a set of SAM dikes starting 
%  from the edge of a sill-like reservoir below a caldera lying on a 
%  coastline and bounded by hills (see Complex-Coastline scenario from
%  Mantiloni et al., 2023). The stress field, due to the gravitational
%  loading/unloading associated with the topography and a superimposed 
%  tectonic stress, is imported from a .mat file (Example2_Data.mat).
%  SAM runs in "Numerical" mode: stresses are calculated through the 
%  Boundary-Element (BE) tool "Cut&Displace" (Davis et al., 2017; 2019) 
%  throughout the dike trajectory simulation. SAM needs the displacement 
%  discontinuities on the BEs, which are pre-calculated and imported from 
%  the .mat file "Example3_Complex_Coastline_Data".
%
%  Colour conventions: SAM forward pathways are represented by red patches,
%  and the final step is marked by a yellow patch. Arrival points projected 
%  to the free surface are shown as by blue dots. 
% 
%  I appreciate any feedback or bug report. Reach out to me at 
%  lorenzo@gfz-potsdam.de
%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  References:
%
%  Davis, T., Healy, D., Bubeck, A., & Walker, R. (2017). Stress 
%  concentrations around voids in three dimensions: The roots of failure. 
%  Journal of Structural Geology, 102, 193-207.
% 
%  Davis, T., Healy, D., & Rivalta, E. (2019). Slip on wavy frictional 
%  faults: Is the 3rd dimension a sticking point?. Journal of Structural  
%  Geology, 119, 33-49.
%
%  Mantiloni, L., Rivalta, E., & Davis, T. (2023). Mechanical Modeling of 
%  Pre-Eruptive Magma Propagation Scenarios at Calderas. Journal of 
%  Geophysical Research: Solid Earth, 128, 3. 
%  https://doi.org/10.1029/2022JB025956

%  Import the stress model, topography mesh and scenario parameters:
load 'Example3_Data.mat'

%  Variables included in "Example3_Data.mat":
%  _ Dds, Dn, Dss: k x 3 arrays with, respectively, dip-slip, normal and
%    strike-slip components of displacement discontinuity on each BE (k =
%    number of BEs)
%  _ mu, nu: rigidity modulus and Poisson's ratio of the host rock
%  _ P1, P2, P3: x, y, z coordinates of the three vertices of each 
%    triangular dislocation
%  _ Points, Triangles: together they describe the mesh of triangular
%    dislocations
%  _ Pr: host rock density (kg/m^3)
%  _ StartPoints: x, y, z coordinates of 10 dike starting points
%  _ TectS: Sxx, Syy, Sxy components of tectonic stress
%  _ Topo_Interp: interpolating function for the topography giving the
%    height of free surface at point (x, y)

%  Setting dike starting points and parameters: here, we simulate only one
%  dike, because running times are longer than those of the other SAM
%  modes. Change 'ImportStartPoint' or 'N' variables below to add more
%  dikes
ImportStartPoint = 0; %Set to 1 to use the imported starting points
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
N = size(StartPoints,1);
c = 800;
n = 12;
Pm = 2300;

%  Creating RockPar and DisplMat arrays (all physical and host rock
%  parameters are imported)
RockPar = [nu mu];
DisplMat = [Dds Dss Dn]; %Dip-slip, strike-slip and normal displacements of
                         %each Boundary Element

%  Setting Direction
Direction = 'F';

%  Setting StressModel
StressModel.StressOption = 'Numerical';
StressModel.RockPar = RockPar;
StressModel.P1 = P1;
StressModel.P2 = P2;
StressModel.P3 = P3;
StressModel.DisplMat = DisplMat;
StressModel.TectS = TectS;
StressModel.MDT = 800;

%  Setting options
%  Setting plotting and propagation/arrest options
display = 1; %Display SAM pathway
RigidSAM = 1; %Preventing SAM dikes from dipping/steering abruptly
displast = 1; %Display the last step of SAM pathways in yellow
options.behavior = [1,1,1];
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

%  Running SAM
[dike_path_tipB,dike_path_dipB,~] = SAM(Direction,StartPoints,c,n,Pr,Pm,...
    Topo_Interp,StressModel,options);

%  Retrieving the arrival poins for each SAM pathway
ArrivalPoints_ProjectedToFreeSurface = zeros(N,3);
ArrivalPoints_Actual = zeros(N,3);
for i=1:numel(dike_path_tipB)
    dike = dike_path_tipB{i};
    PathwayDifference = size(dike,1) - size(dike_path_dipB{i},2);
    ArrivalPoints_ProjectedToFreeSurface(i,:) = dike(end,:);
    ArrivalPoints_Actual(i,:) = dike(end-PathwayDifference,:);
end
