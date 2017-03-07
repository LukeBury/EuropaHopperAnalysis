clear
clc
close all
addpath('../ProjectBin')

% ------------------------------------------------------------------------
%%% System Parameters
% ------------------------------------------------------------------------
a = 671100; % semimajor axis (distance-normalizing value)
rad1 = 69911/a; % radius of primary body
rad2 = 1560.8/a; % radius of secondary body
w = [0, 0, 1]; % angular velocity

%%% Gravitational Parameters
u1 = 126672520; % km^3 / s^2
u2 = 3203.413216; % km^3 / s^2
% ------------------------------------------------------------------------
%%% Defining Planetary System
% ------------------------------------------------------------------------
% massPrimary = 10;
% massSecondary = 2;
% 
% u = massSecondary/massPrimary;
u = u2/u1;

%%% Rotating frame cooridinates
rB1_BCR = [-u, 0, 0];
rB2_BCR = [1-u, 0, 0];
 
% ------------------------------------------------------------------------
%%% Defining Particle State
% ------------------------------------------------------------------------
rH0_BCR = [1-u, 1.1*rad2, 0];
% velPart0 = cross(w,posPart0);
% posPart0 = [1, .2, 0];
vH0_BCR = [0, 0, 0];

r10 = norm(rB1_BCR - rH0_BCR);
r20 = norm(rB2_BCR - rH0_BCR);

% ------------------------------------------------------------------------
%%% Propagating State
% ------------------------------------------------------------------------
%%% Setting time vector
ti = 0; % sec
dt = .001; % sec
tf = 2; % sec
time = ti:dt:tf; % sec

%%% Choosing ode45 tolerance
tol = 1e-13;

%%% Setting integrator options
options = odeset('Events',@normalCR3BP_impactEvent,'RelTol',tol,'AbsTol',tol);

%%% Setting Initial State Vector (ECEF)
X0 = [rH0_BCR, vH0_BCR]; % km, km/s

%%% Propagating the State
[Times,States] = ode45(@normalCR3BP_Int,time,X0,options,u,rB1_BCR,rB2_BCR,rad2);

%%% Assigning variables to state components
rH_BCR = States(:,1:3);
vH_BCR = States(:,4:6);

% ------------------------------------------------------------------------
%%% Jacobi Constants
% ------------------------------------------------------------------------
[JCs] = JacobiConstantCalculator(u,rH_BCR,vH_BCR,rB1_BCR,rB2_BCR);

% ------------------------------------------------------------------------
%%% Jacobi Constant Contours
% ------------------------------------------------------------------------

x = linspace(0.8*(1-u), 1.2*(1-u), 100);
y = linspace(-.1, .1, 100);
% x = linspace(-.1, 1.2*(1-u), 100);
% y = linspace(-u, u, 100);
z = zeros(length(x),length(y));
for xk = 1:length(x)
    for yk = 1:length(y)
        r1 = norm(rB1_BCR - [x(xk), y(yk), 0]);
        r2 = norm(rB2_BCR - [x(xk), y(yk), 0]);
        z(xk,yk) = x(xk)^2 + y(yk)^2 + 2*((1-u)/r1 + u/r2);
    end
end

figure
hold all
contourVals = linspace(3,3.02,10);
contour(x,y,z,[contourVals, 3.01966,3.1],'ShowText','on')
PlotBoi3('X','Y','Z',10)
x_temp = rad2 * cos(0:.01:2*pi) + rB2_BCR(1);
y_temp = rad2 * sin(0:.01:2*pi);
plot(x_temp, y_temp,'b');
x_temp = rad1 * cos(0:.01:2*pi) + rB1_BCR(1);
y_temp = rad1 * sin(0:.01:2*pi);
plot(x_temp, y_temp,'r');
axis equal

% 
% 
% figure
% plot(Times,(JCs-JCs(1))*100./(JCs(1)),'.','linewidth',1.5)
% PlotBoi2('Time, sec','JC Percent Change',16)
% 
% figure; hold all
% x = rad1 * cos(0:.01:2*pi) + rB1_BCR(1);
% y = rad1 * sin(0:.01:2*pi);
% plot(x, y,'r');
% x = rad2 * cos(0:.01:2*pi) + rB2_BCR(1);
% y = rad2 * sin(0:.01:2*pi);
% plot(x, y,'b');
% 
% plot3(States(:,1),States(:,2),States(:,3))
% PlotBoi3('x','y','z',15)