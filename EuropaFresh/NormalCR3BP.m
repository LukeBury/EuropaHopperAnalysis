clear
clc
close all
addpath('../ProjectBin')

% ------------------------------------------------------------------------
%%% Plot/Run Options; 0 = no, 1 = yes
% ------------------------------------------------------------------------
plotRotatingSystem  = 1;
plot3DEuropa        = 1;
plotJCPercentChange = 0; 
motion2DPlot        = 1;
plotJCContours      = 0;
latlonReachability  = 0;

% ------------------------------------------------------------------------
%%% System Parameters
% ------------------------------------------------------------------------
a = 671100; % semimajor axis (distance-normalizing value)
rad1_n = 69911/a; % radius of primary body (km normalized)
rad2_n = 1560.8/a; % radius of secondary body (km normalized)
% rad1_n = 6371/a; % Earth
% rad2_n = 1737/a; % Moon
w2 = 2.047200349303344e-05; % rad/s
w_n = [0, 0, w2/w2]; % angular velocity

%%% Gravitational Parameters
u1 = 126672520; % km^3 / s^2
u2 = 3203.413216; % km^3 / s^2
% u1 = 398600.4415; % Earth
% u2 = 4904.86; % Moon
% u1 = 10;
% u2 = 0.5;
u = u2/u1;

%%% Rotating frame cooridinates
rB1_BCR_n = [-u, 0, 0];
rB2_BCR_n = [1-u, 0, 0];

%%% Normalizing Factors (real / xN = normalized)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Never felt awesome about these
wNorm = 1/w2; % Angular Velocity (rad/s)
rNorm = a; % Position (km)
vNorm = rNorm/wNorm; % Velocity (km/s)
tNorm = wNorm;

% ------------------------------------------------------------------------
%%% Defining Particle State
% ------------------------------------------------------------------------
%%% Initial Particle Position
% Surface Position (latitude / longitude)
lat0 = -90; % deg (-90:90)
lon0 = -90; % deg (-180:180)
[rH0_ECEF_n] = latlon2surfECEF(lat0, lon0, rad2_n);
rH0_BCR_n = rH0_ECEF_n + rB2_BCR_n;

%%% Initial Partical Velocity
v_mag_n = 1/vNorm; % Velocity magnitude (km/s normalized)
% vH0rad_BCR_n = (rH0_ECEF_n/norm(rH0_ECEF_n))*v_mag_n; % radial velocity
az = 0;
el = 90;
[vH0_BCR_n] = angleVelocity(v_mag_n, az, el, lat0, lon0);

r10_n = norm(rB1_BCR_n - rH0_BCR_n);
r20_n = norm(rB2_BCR_n - rH0_BCR_n);

% ------------------------------------------------------------------------
%%% Propagating State
% ------------------------------------------------------------------------
%%% Setting normalized time vector
ti = 0*w2;
dt = 1*w2;
tf = 50000*w2;
time = ti:dt:tf;

%%% Choosing ode45 tolerance
tol = 1e-13;

%%% Setting integrator options
options = odeset('Events',@normalCR3BP_impactEvent,'RelTol',tol,'AbsTol',tol);

%%% Setting Initial State Vector (ECEF)
X0_n = [rH0_BCR_n, vH0_BCR_n]'; % km, km/s

%%% Propagating the State
[Times_n,States_n] = ode45(@normalCR3BP_Int,time,X0_n,options,u,rB1_BCR_n,rB2_BCR_n,rad2_n);


%%% To get real times
Times = Times_n.*tNorm;

%%% Assigning variables to state components
rH_BCR_n = States_n(:,1:3);
vH_BCR_n = States_n(:,4:6);
% [0 1 1]
if plotRotatingSystem == 1
    figure; hold all
    bodySurface2(rad1_n, rB1_BCR_n, [1 .5 0],1.5)
    bodySurface2(rad2_n, rB2_BCR_n, [0 1 1], 1.5)
    plot3(States_n(:,1),States_n(:,2),States_n(:,3))
    PlotBoi3('X','Y','Z',16)
    axis equal
end

% ------------------------------------------------------------------------
%%% 3D Europa w/ Trajectory Plot
% ------------------------------------------------------------------------
if plot3DEuropa == 1
    figure; hold all
    plot3(States_n(:,1),States_n(:,2),States_n(:,3),'linewidth',1.5)
    bodySurface3(rad2_n, rB2_BCR_n, [0 1 1])
    quiver3(rB2_BCR_n(1),rB2_BCR_n(2),rB2_BCR_n(3),-rB2_BCR_n(1), -rB2_BCR_n(2), -rB2_BCR_n(3),.004)
    PlotBoi3('X','Y','Z',16)
    axis equal
end



% ------------------------------------------------------------------------
%%% Jacobi Constants
% ------------------------------------------------------------------------
[JCs_n] = JacobiConstantCalculator(u,rH_BCR_n,vH_BCR_n);

if plotJCPercentChange == 1
    figure
    plot(Times,(JCs_n-JCs_n(1))*100./(JCs_n(1)),'-o','linewidth',1.5)
    PlotBoi2('Time, sec','JC Percent Change',16)
end

% ------------------------------------------------------------------------
%%% 2D Motion
% ------------------------------------------------------------------------
if motion2DPlot == 1
    xmin = min(States_n(:,1))-(v_mag_n*vNorm*100)/rNorm; % min minus 1 km
    xmax = max(States_n(:,1))+(v_mag_n*vNorm*100)/rNorm; % max plus 1 km
    ymin = min(States_n(:,2))-(v_mag_n*vNorm*100)/rNorm;
    ymax = max(States_n(:,2))+(v_mag_n*vNorm*100)/rNorm;
    
    x = linspace(xmin, xmax, 1000);
    y = linspace(ymin, ymax, 1000);
    [X, Y] = meshgrid(x,y);
    clear x y

    Z = zeros(size(X));
    for xk = 1:size(X,1)
        for yk = 1:size(X,2)
            r1 = norm(rB1_BCR_n - [X(xk,yk), Y(xk,yk), 0]);
            r2 = norm(rB2_BCR_n - [X(xk,yk), Y(xk,yk), 0]);
            Z(xk,yk) = X(xk,yk)^2 + Y(xk,yk)^2 + 2*((1-u)/r1 + u/r2);
        end
    end

    figure; hold all
    contour(X,Y,Z,[max(JCs_n), JCs_n(1)+.01],'linewidth',2,'ShowText','on')
    plot3(States_n(:,1),States_n(:,2),States_n(:,3))
    bodySurface2(rad2_n, rB2_BCR_n, [0 1 1], 2)
    plot3(States_n(1,1),States_n(1,2),States_n(1,3),'ro','markersize',10)
    plot3(States_n(end,1),States_n(end,2),States_n(end,3),'rx','markersize',10)
    PlotBoi3('X','Y','Z',10)
    axis equal
    xlim([xmin xmax])
    ylim([ymin ymax])
%     axis equal
end


% ------------------------------------------------------------------------
%%% Jacobi Constant Contours
% ------------------------------------------------------------------------
if plotJCContours == 1
    x = linspace((1-u)-1.2*rad2_n, (1-u)+1.2*rad2_n, 100);
    y = linspace(-1.2*rad2_n, 1.2*rad2_n, 100);
    [X, Y] = meshgrid(x,y);
    clear x y

    Z = zeros(size(X));
    for xk = 1:size(X,1)
        for yk = 1:size(X,2)
            r1 = norm(rB1_BCR_n - [X(xk,yk), Y(xk,yk), 0]);
            r2 = norm(rB2_BCR_n - [X(xk,yk), Y(xk,yk), 0]);
            Z(xk,yk) = X(xk,yk)^2 + Y(xk,yk)^2 + 2*((1-u)/r1 + u/r2);
        end
    end
    
    figure
    meshc(X,Y,Z)
    
    figure
    hold all
    contourVals = linspace(3.01,3.0212,10);
    contour(X,Y,Z,[contourVals, JCs_n(1)],'ShowText','on')
    
    PlotBoi3('X','Y','Z',16)
    bodySurface2(rad1_n, rB1_BCR_n, [1 .5 0], 2)
    bodySurface2(rad2_n, rB2_BCR_n, [0 1 1], 2)
    axis equal
    xlim([(1-u)-1.2*rad2_n, (1-u)+1.2*rad2_n]); ylim([-1.2*rad2_n, 1.2*rad2_n])    
end

% ------------------------------------------------------------------------
%%% Latitude-Longitude Reachability
% ------------------------------------------------------------------------
if latlonReachability == 1
    fprintf('*Insert smart things\n')
    
end



