%  Calculates and plots the pathways of SAM dykes starting below 
%  an axisymmetric shallow topographic depression. Stress due to the
%  surface unloading associated to the depression is computed analytically 
%  as a distribution of normal forces (see function 
%  "HalfSpaceNormalForce.m" in folder "Shared_functions").
%
%  Colour conventions: SAM forward pathways are represented by red patches,
%  and the final step is marked by a yellow patch. Arrival points are shown
%  as blue dots. SAM backtracked trajectories are represented by green 
%  dots, and the points where they stop are also marked by blue dots.
% 
%  I appreciate any feedback or bug report. Reach out to me at 
%  lorenzo@gfz-potsdam.de
%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  References:
% 
%  Jaeger, J. C., Cook, N. G. W., & Zimmermann, R. W. (2007). Fundamentals 
%  of rock mechanics, edited by Chapman and Hall. New York, United States 
%  of America, 1-593.
%
%  Mantiloni, L., Rivalta, E., & Davis, T. (2023). Mechanical Modeling of 
%  Pre-Eruptive Magma Propagation Scenarios at Calderas. Journal of 
%  Geophysical Research: Solid Earth, 128, 3. 
%  https://doi.org/10.1029/2022JB025956

%  Physical and host rock parameters
g = 9.81;
nu = 0.25; %Poisson's ratio
Pr = 2850; %Host rock density (kg/m^3)

%  Unloading geometry and parameters
Depth = 3e3;
P = -Pr*g*Depth;
Radius = 1e3;

%  Number of point forces
Npf = 200;

%  Grid of point forces
Xloc = linspace(-Radius,Radius,Npf);
Yloc = Xloc;
[Xloc,Yloc] = meshgrid(Xloc,Yloc);
Xloc = Xloc(:);
Yloc = Yloc(:);
condOutofRadius = (Xloc.^2 + Yloc.^2) > Radius^2;
Xloc = Xloc(~condOutofRadius);
Yloc = Yloc(~condOutofRadius);

%  Setting dyke starting points and parameters
N = 4; %Number of dykes
th=0:2*pi/N:2*pi-2*pi/N;
rad=3e3;
StartX = -rad * cos(th);
StartY = rad * sin(th);
StartZ = StartX*0 - 6e3;
StartPoints = [StartX' StartY' StartZ'];
c = 800;
n = 12;
Pm = 2200; %Magma density

%  Setting Direction
Direction = 'F';

%  Setting StressModel
StressModel.StressOption = 'Analytical';
StressModel.Stressfun = @HalfSpaceNormalForce;
StressModel.Stressfun_Par = {P,Pr,nu,Xloc,Yloc}; 

%  Setting plotting and propagation/arrest options
display = 1; %Display SAM pathway
RigidSAM = 1; %Preventing SAM dykes from dipping/steering abruptly
displast = 1; %Display the last step of SAM pathways in yellow
options.behavior = [display,RigidSAM,displast]; 

%  Setting grid to plot the free surface
Xsurf = linspace(-2,2,40)*1e4; 
Ysurf = Xsurf; 
[Xsurf,Ysurf] = meshgrid(Xsurf,Ysurf); 
Zsurf = Xsurf*0;  

%  Assigning the free surface elevation (here, a flat surface at z=0)
Topo_Interp = 0; 

%  Plotting free surface and distribution of point forces
figure
hold on
scatter(Xloc,Yloc,'b')
quiver3(Xloc,Yloc,Xloc*0,Xloc*0,Xloc*0,Xloc*0-sign(P),25,'b')
surf(Xsurf,Ysurf,Zsurf,'facecolor','g','FaceAlpha',0.8,'LineStyle','none');
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
xlim([min(Xsurf(:)),max(Xsurf(:))]);
ylim([min(Ysurf(:)),max(Ysurf(:))]);
axis equal
zlim([-9e3,1e3])

%  Running SAM: 

[dyke_path_tip,dyke_path_dip,dyke_path_strike] = SAM(Direction,...
    StartPoints,c,n,Pr,Pm,Topo_Interp,StressModel,options);

% Retrieving the arrival poins for each SAM pathway:
ArrivalPoints = zeros(N,3);
for i=1:numel(dyke_path_tip)
    dyke = dyke_path_tip{i};
    ArrivalPoints(i,:) = dyke(end-1,:);
end

%  Backtracking each pathway from the respective arrival point 
Direction = 'B';
%  Here we employ the same radius c as in the forward model, but they may 
%  be different
cB = c; 
%  We set the vertical threshold slightly deeper than the starting depth 
%  of dykes
options.cutoff = [1e4,-6.01e3];

%  Run backtrack SAM 
[dyke_path_tip1,dyke_path_dip1,dyke_path_strike1,BSP1] = SAM(Direction,...
    ArrivalPoints,cB,n,Pr,Pm,Topo_Interp,StressModel,options,Pr);

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
