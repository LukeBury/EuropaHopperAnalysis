clear
clc
close all
addpath('../ProjectBin')

% ------------------------------------------------------------------------
%%% Plot Options
% ------------------------------------------------------------------------
rH_ECEF_Diff_Plots = 0; % 0:off   1:on
Jacobian_Plots = 1;     % 0:off   1:on

% ------------------------------------------------------------------------
%%% Setting Initial Values
% ------------------------------------------------------------------------
%%% Time Constraints
ti = 0;
% tf = 520000; % sec
tf = 5327*2;
dt = .1;
time = ti:dt:tf;

%%% Jupiter Parameters
uJ = 126672520; % km^3 / s^2

%%% Europa Parameters
RE = 1560.8; % km
aE = 671100; % km
uE = 3203.413216; % km^3 / s^2
nE = sqrt(uJ / (aE^3)); % rad/s
vE = sqrt(uJ/aE); % Circular Velocity, km/s
wE = [0; 0; nE]; % Rot. velocity (tidally locked)

%%% Intial Europa State
E_theta0 = 0*(pi/180); % Initial position of Europa about Jupiter from +x, rads
rE0_JCI = R3([aE, 0, 0],E_theta0); % km
vE0_JCI = R3([0, vE, 0],E_theta0); % km

%%% Initial Hopper State
% Surface Position (latitude / longitude)
lat1 = 0; % deg (-90:90)
lon1 = 0; % deg (-180:180)
[rH0_ECEF] = latlon2surfECEF(lat1, lon1, RE); % km
% rH0_ECEF = [RE+500, 0, 0];%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for orbit
rH0_ECI = R3(rH0_ECEF,E_theta0); % km
rH0_JCI = rH0_ECI + rE0_JCI; % km

%%% Radial Velocity (Europa relative)
v_mag = .015; % km/s
vH0_ECEF = (rH0_ECEF/norm(rH0_ECEF))*v_mag;
% vH0_ECEF = [0, sqrt(uE/(RE+500)),0] - cross(wE,rH0_ECI); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for orbit

%%% Creating ECI and JCI Initial Hopper Velocity
vH0_ECI = R3(vH0_ECEF,E_theta0) + cross(wE,rH0_ECI); % km/s
vH0_JCI = vH0_ECI + vE0_JCI; % km/s

% ------------------------------------------------------------------------
%%% Propagating the Inertial State with Numerical Integration
% ------------------------------------------------------------------------
%%% Setting integrator options
tol = 1E-10;
optionsI = odeset('Events',@impactEvent_I_exp,'RelTol',tol,'AbsTol',tol);

%%% Setting Initial State Vector (JCI)
X0_JCI = [rH0_JCI, vH0_JCI]; % km km/s km km/s

%%% Propagating the State
[TimesI,StatesI] = ode45(@EH_NumIntegrator_CR3BP_exp,time,X0_JCI,optionsI,RE,uE,uJ,nE,aE,E_theta0);  % make RK4 integrator (fixed step), ode113 (a couple others as well), Plot total energy

%%% Calculating Europa States
rE_JCI = zeros(length(TimesI),3);
vE_JCI = zeros(length(TimesI),3);
rotAngles = zeros(length(TimesI),1);
rBC_JCI = zeros(length(TimesI),3);
for k = 1:length(TimesI)
    rotAngles(k) = nE * TimesI(k); % rads
    rE_JCI(k,:) = R3(rE0_JCI,rotAngles(k)); % km, JCI Europa position vectors
    vE_JCI(k,:) = R3(vE0_JCI,rotAngles(k)); % km, JCI Europa position vectors
    rBC_JCI(k,:) = R3([aE, 0, 0].*(uE/uJ),rotAngles(k)); % km, JCI barycenter
end

%%% Creating positional and velocity matrices with inertial results
rH_JCI_I = StatesI(:,1:3); % km
vH_JCI_I = StatesI(:,4:6); % km/s

rH_ECI_I = rH_JCI_I(:,1:3) - rE_JCI; % km

rH_ECEF_I = zeros(length(TimesI),3);
for k = 1:length(TimesI)
    rH_ECEF_I(k,:) = R3(rH_ECI_I(k,:),-rotAngles(k)); % km
end

% ------------------------------------------------------------------------
%%% Propagating the Body-Frame State with Numerical Integration
% ------------------------------------------------------------------------
%%% Setting integrator options
optionsB = odeset('Events',@impactEvent_B_exp,'RelTol',tol,'AbsTol',tol);

%%% Setting Initial State Vector (ECEF)
X0_ECEF = [rH0_ECEF, vH0_ECEF]; % km, km/s

%%% Propagating the State
[TimesB,StatesB] = ode45(@EH_Body_NumIntegrator_CR3BP_exp,time,X0_ECEF,optionsB,RE,uE,uJ,nE,aE,E_theta0);

%%% Checking that solutions are same length
if length(TimesI) ~= length(TimesB)
    warning('Time discrepancy in simulations')
end

%%% Creating positional and velocity matrices with body results
rH_ECEF_B = StatesB(:,1:3); % km
vH_ECEF_B = StatesB(:,4:6); % km/s

rH_ECI_B = zeros(size(rH_ECEF_B));
rH_JCI_B = zeros(size(rH_ECI_B));
for k = 1:length(TimesB)
    rH_ECI_B(k,:) = R3(rH_ECEF_B(k,:),rotAngles(k)); % km
    rH_JCI_B(k,:) = rH_ECI_B(k,:) + rE_JCI(k,:); % km
end

% ------------------------------------------------------------------------
%%% Jacobi Constant Analysis
% ------------------------------------------------------------------------
JC_I = zeros(length(TimesI),1);
JC_B = zeros(length(TimesI),1);

%%% Initializing matrices to hold individual components of JC calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xs_I = zeros(length(TimesI),1);
xs_B = zeros(length(TimesI),1);
ys_I = zeros(length(TimesI),1);
ys_B = zeros(length(TimesI),1);
r1s_I = zeros(length(TimesI),1);
r1s_B = zeros(length(TimesI),1);
r2s_I = zeros(length(TimesI),1);
r2s_B = zeros(length(TimesI),1);
vs_I = zeros(length(TimesI),1);
vs_B = zeros(length(TimesI),1);
for k = 1:length(TimesI)
    %%% Non-Normalized Method (Inertial Results)
    rBC_I = [aE, 0, 0].*(uE/uJ); % Barycenter position (JC_Rot)
    x_I = rH_ECEF_I(k,1) + aE - rBC_I(1); % (BC_Rot)
    y_I = rH_ECEF_I(k,2); % (BC_Rot)
    r1_I = norm(rH_ECEF_I(k,:)+[aE, 0, 0]);  % (Jupiter Distance)
    r2_I = norm(rH_ECEF_I(k,:)); % (Europa Distance)
%     vH_BC_I = vH_JCI_I(k,:) - vE_JCI(k,:) - cross(wE', rH_ECI_I(k,:));
    vH_BC_I = R3(vH_JCI_I(k,:) - cross(wE',rBC_JCI(k,:)), -rotAngles(k)) - cross(wE', [x_I, y_I, 0]);
    
    % Storing JC Values
    JC_I(k) = (nE^2)*(x_I^2 + y_I^2) + 2*(uJ/r1_I + uE/r2_I) - norm(vH_BC_I)^2;
    
    %%% Non-Normalized Method (Body Results)
    rBC_B = [aE, 0, 0].*(uE/uJ); % Barycenter position (JC_Rot)
    x_B = rH_ECEF_B(k,1) + aE - rBC_B(1); % (BC_Rot)
    y_B = rH_ECEF_B(k,2); % (BC_Rot)
    r1_B = norm(rH_ECEF_B(k,:)+[aE, 0, 0]); % (Jupiter Distance)
    r2_B = norm(rH_ECEF_B(k,:)); % (Europa Distance)
    vH_BC_B = vH_ECEF_B(k,:); % (BC_Rot)
    
    % Storing JC Values
    JC_B(k) = (nE^2)*(x_B^2 + y_B^2) + 2*(uJ/r1_B + uE/r2_B) - norm(vH_BC_B)^2;
    
    %%% Grabbing JC components for comparison
    xs_I(k) = (nE^2)*(x_I^2);
    xs_B(k) = (nE^2)*(x_B^2);
    ys_I(k) = (nE^2)*(y_I^2);
    ys_B(k) = (nE^2)*(y_B^2);
    r1s_I(k) = 2*uJ/r1_I;
    r1s_B(k) = 2*uJ/r1_B;
    r2s_I(k) = 2*uE/r2_I;
    r2s_B(k) = 2*uE/r2_B;
    vs_I(k) = norm(vH_BC_I)^2;
    vs_B(k) = norm(vH_BC_B)^2;
    
    
end

%%% Clearing Temporary Variables
clear x_I y_I rBC_I vBC_I vJC_I r1_I r2_I
clear x_B y_B rBC_B vBC_B vJC_B r1_B r2_B

% ------------------------------------------------------------------------
%%% Plots
% ------------------------------------------------------------------------
%%% Plotting differences in rH_ECEF components
if rH_ECEF_Diff_Plots == 1
    figure
    subplot(3,1,1)
    plot(TimesI,(rH_ECEF_I(:,1) - rH_ECEF_B(:,1)),'m.','linewidth',1.5)
    title('rH\_ECEF\_I - rH\_ECEF\_B')
    PlotBoi2('','rX diff, km',16)
    subplot(3,1,2)
    plot(TimesI,(rH_ECEF_I(:,2) - rH_ECEF_B(:,2)),'m.','linewidth',1.5)
    PlotBoi2('','rY diff, km',16)
    subplot(3,1,3)
    plot(TimesI,(rH_ECEF_I(:,3) - rH_ECEF_B(:,3)),'m.','linewidth',1.5)
    PlotBoi2('Times, sec','rZ diff, km',16)
end

%%% Plotting Jacobian Constants
if Jacobian_Plots == 1
    figure; hold all
    plot(TimesI,JC_I,'r','linewidth',1.5)
    plot(TimesB,JC_B,'b','linewidth',1.5)
    legend('JC\_I','JC\_B')
    PlotBoi2('Time, sec','Jacobian Constant',16)
    
    figure; hold all
    plot(TimesI,(JC_I-JC_I(1))*100./(JC_I(1)),'linewidth',1.5)
    PlotBoi2('Time, sec','JC Percent Change',16)

    figure
    plot(TimesI,JC_I-JC_B,'m','linewidth',1.5)
    legend('JC\_I - JC\_B')
    PlotBoi2('Time, sec','Jacobian Constant Difference',16)

end

%%% Plotting differences of various components of JC calculation
% figure; hold all;
% plot(TimesI, xs_I - xs_B, 'm.')
% title('xs_I - xs_B')
% 
% figure; hold all;
% plot(TimesI, ys_I - ys_B, 'm.')
% title('ys_I - ys_B')
% 
% figure; hold all;
% plot(TimesI, r1s_I - r1s_B, 'm.')
% title('r1s_I-r1s_B')
% 
% figure; hold all;
% plot(TimesI, r2s_I - r2s_B, 'm.')
% title('r2s_I-r2s_B')
% 
% figure; hold all
% plot(TimesI,vs_I - vs_B,'m.','linewidth',1.5)
% title('Jacobian Velocity Differences')
% PlotBoi2('Times, sec','v^2(z) diff',16)


%%% Plotting ECI Motion
figure; hold all
plot3(rH_ECI_I(:,1),rH_ECI_I(:,2),rH_ECI_I(:,3),'r')
plot3(rH_ECI_B(:,1),rH_ECI_B(:,2),rH_ECI_B(:,3),'b')
%%% Plotting Europa equator
x = RE * cos(0:.01:2*pi);
y = RE * sin(0:.01:2*pi);
plot(x, y,'b');
title('ECI')
axis equal

%%% Plotting ECEF Motion
figure; hold all
plot3(rH_ECEF_I(:,1),rH_ECEF_I(:,2),rH_ECEF_I(:,3),'r')
plot3(rH_ECEF_B(:,1),rH_ECEF_B(:,2),rH_ECEF_B(:,3),'b')
%%% Plotting Europa equator
x = RE * cos(0:.01:2*pi);
y = RE * sin(0:.01:2*pi);
plot(x, y,'b');
title('ECEF')
axis equal








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Integrators   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------------------------------------------------------------
%%% Inertial Numerical Integrator
% ------------------------------------------------------------------------
function [ dY ] = EH_NumIntegrator_CR3BP_exp(t,Y,RE,uE,uJ,nE,aE,E_theta0)
dY = zeros(6,1);

%%% Creating Europa Position (JCI)
ta = nE*t; % rads
rE = R3([aE; 0; 0],E_theta0); % km
rE = R3(rE,ta); % km

%%% Unpack the Hopper state vector (JCI)
yH = Y(1:3); % Hopper Position, km
dyH = Y(4:6); % Hopper Velocity, km/s

%%% Europa-Centric hopper position (ECI)
rH = yH - rE; % km

%%% Hopper Dynamics
fJ = (-uJ/(norm(yH)^3))*yH; % Jupiter Pull, km/s^2
fE = (-uE/(norm(rH)^3))*rH; % Europa Pull, km/s^2
ddyH = fJ + fE; % Jupiter and Europa pull km/s^2

%%% Output the derivative of the state
dY(1:3) = dyH; % km/s
dY(4:6) = ddyH; % km/s^2
end


% ------------------------------------------------------------------------
%%% Body Frame Numerical Integrator
% ------------------------------------------------------------------------
function [ dY ] = EH_Body_NumIntegrator_CR3BP_exp(t,Y,RE,uE,uJ,nE,aE,E_theta0)
dY = zeros(6,1);

%%% Unpack the Hopper state vector (ECEF)
yH = Y(1:3); % Hopper Position, km
dyH = Y(4:6); % Hopper Velocity, km/s

%%% Creating Europa Position
ta = nE*t; % rads
rE_JCI = R3([aE; 0; 0],E_theta0); % km
rE_JCI = R3(rE_JCI,ta); % km

%%% Creating Hopper Position (ECI)
rH_ECI = R3(yH,ta); % km

%%% Creating Hopper Position (JCI)
rH_JCI = rH_ECI + rE_JCI; % km

%%% Determining Inertial Accelerations (JCI)
aH_JCI = (-uJ/(norm(rH_JCI)^3))*rH_JCI...
    + (-uE/(norm(rH_ECI)^3))*rH_ECI; % km/s^2
aH_ECEF = R3(aH_JCI, -ta); % km/s^2

%%% Determining Acceleration of Europa (JCI)
aE_JCI = (-uJ/(norm(rE_JCI)^3))*rE_JCI;
aE_ECEF = R3(aE_JCI, -ta); % km/s^2

%%% Determining Body Frame Hopper Acceleration (ECEF)
ddyH = aH_ECEF - aE_ECEF - 2*cross([0;0;nE],dyH) - cross([0;0;nE],cross([0;0;nE],yH)); % km/s^2 (Shaub, pg 19)

%%% Output the derivative of the state
dY(1:3) = dyH; % km/s
dY(4:6) = ddyH; % km/s^2
end


% ------------------------------------------------------------------------
%%% Inerital Frame Impact Event
% ------------------------------------------------------------------------
function [value, isterminal, direction] = impactEvent_I_exp(t,Y,RE,uE,uJ,nE,aE,E_theta0)
%%% Creating Europa Position (JCI)
ta = nE*t; % rads
rE = R3([aE; 0; 0],E_theta0); % km
rE = R3(rE,ta); % km

%%% Unpack the Hopper position vector (JCI)
yH = Y(1:3); % km

%%% Europa-Centric hopper position (ECI)
rH = yH - rE; % km

%%% Event function watching for when "value" = 0 (Hopper impacts Europa)
value = norm(rH) - RE;
isterminal = 1; % stops the integration
direction = -1; % negative direction only
end


% ------------------------------------------------------------------------
%%% Body Frame Impact Event
% ------------------------------------------------------------------------
function [value, isterminal, direction] = impactEvent_B_exp(t,Y,RE,uE,uJ,nE,aE,E_theta0)
%%% Unpack the Hopper position vector (ECEF)
yH = Y(1:3);

%%% Event function watching for when "value" = 0 (Hopper impacts Europa)
value = norm(yH) - RE;
isterminal = 1; % stops the integration
direction = -1; % negative direction only
end






