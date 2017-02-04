clear
clc
close all
addpath('ProjectBin')
% ------------------------------------------------------------------------
%%% Setting Initial Values
% ------------------------------------------------------------------------
%%% Time Constraints
ti = 0;
tf = 520000; % sec
dt = .1;
time = ti:dt:tf;

%%% Jupiter Parameters
uJ = 126672520; % km^3 / s^2

%%% Europa Parameters
E_radius = 1560.8; % km
E_a = 671100; % km
uE = 3203.413216; % km^3 / s^2
nE = sqrt(uJ / (E_a^3)); % rad/s
% Circular Velocity
vE = sqrt(uJ/E_a); % km/s
wE = [0; 0; nE]; % Rot. velocity (tidally locked)

%%% Intial Europa State
E_theta0 = 0*(pi/180); % Initial position of Europa about Jupiter from +x, rads
rE0_JCI = R3([E_a, 0, 0],E_theta0); % km
vE0_JCI = R3([0, vE, 0],E_theta0); % km

%%% Initial Hopper State
% Surface Position (latitude / longitude)
lat1 = 0; % deg (-90:90)
lon1 = -45; % deg (-180:180)
[rH0_ECEF] = latlon2surfECEF(lat1, lon1, E_radius); % km
rH0_ECI = R3(rH0_ECEF,E_theta0); % km
rH0_JCI = rH0_ECI + rE0_JCI; % km

%%% Radial Velocity (Europa relative)
v_mag = .015; % km/s
vH0_ECEF = (rH0_ECEF/norm(rH0_ECEF))*v_mag;

%%% Creating ECI and JCI Initial Hopper Velocity
vH0_ECI = R3(vH0_ECEF,E_theta0) + cross(wE,rH0_ECI); % km/s
vH0_JCI = vH0_ECI + vE0_JCI; % km/s


% ------------------------------------------------------------------------
%%% Propagating the Inertial State with Numerical Integration
% ------------------------------------------------------------------------
%%% Setting integrator options
tol = 1E-13;
optionsI = odeset('Events',@impactEvent_I_exp,'RelTol',tol,'AbsTol',tol);

%%% Setting Initial State Vector (JCI)
X0_JCI = [rH0_JCI vH0_JCI]; % km km/s km km/s

%%% Propagating the State
[TimesI,StatesI] = ode45(@EH_NumIntegrator_CR3BP_exp,time,X0_JCI,optionsI,E_radius,uE,uJ,nE,E_a);

%%% Calculating Europa States
rE_JCI = zeros(length(TimesI),3);
rotAngles = zeros(length(TimesI),1);
for k = 1:length(TimesI)
    rotAngles(k) = nE * TimesI(k); % rads
    rE_JCI(k,:) = R3(rE0_JCI,rotAngles(k)); % km, JCI Europa position vectors
end

%%% Creating all positional and velocity matrices
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
X0_ECEF = [rH0_ECEF vH0_ECEF]; % km, km/s

%%% Propagating the State
[TimesB,StatesB] = ode45(@EH_Body_NumIntegrator_CR3BP_exp,time,X0_ECEF,optionsB,E_radius,uE,uJ,nE,E_a);

%%% Creating all positional and velocity matrices
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

vs_I = zeros(length(TimesI),3);
vs_B = zeros(length(TimesI),3);
for k = 1:length(TimesI)
    %%% Non-Normalized Method (Inertial Results)
    rBC_I = [E_a, 0, 0].*(uE/uJ); % Barycenter position (JC_Rot)
    x_I = rH_ECEF_I(k,1) + E_a - rBC_I(1); % (BC_Rot)
    y_I = rH_ECEF_I(k,2); % (BC_Rot)
    r1_I = norm(rH_ECEF_I(k,:)+[E_a, 0, 0]);  % (Jupiter Distance)
    r2_I = norm(rH_ECEF_I(k,:)); % (Europa Distance)
    vH_BC_I = vH_JCI_I(k,:) - cross(wE', rE_JCI(k,:)); % (BC_Rot vH)
    
    % Storing Values
    JC_I(k) = (nE^2)*(x_I^2 + y_I^2) + 2*(uJ/r1_I + uE/r2_I) - norm(vH_BC_I)^2;
    vs_I(k,:) = vH_BC_I;
    
    
    %%% Non-Normalized Method (Body Results)
    rBC_B = [E_a, 0, 0].*(uE/uJ); % Barycenter position (JC_Rot)
    x_B = rH_ECEF_B(k,1) + E_a - rBC_B(1); % (BC_Rot)
    y_B = rH_ECEF_B(k,2); % (BC_Rot)
    r1_B = norm(rH_ECEF_B(k,:)+[E_a, 0, 0]); % (Jupiter Distance)
    r2_B = norm(rH_ECEF_B(k,:)); % (Europa Distance)
    vH_BC_B = vH_ECEF_B(k,:); % (BC_Rot)
    
    % Storing Values
    JC_B(k) = (nE^2)*(x_B^2 + y_B^2) + 2*(uJ/r1_B + uE/r2_B) - norm(vH_BC_B)^2;
    vs_B(k,:) = vH_BC_B;
    
    
%     %%% Plotting differences of various components of JC calculation
%     figure(5); hold all;
%     plot(TimesI(k), x_I-x_B, 'r.')
%     
%     figure(6); hold all;
%     plot(TimesI(k), y_I-y_B, 'r.')
%     
%     figure(7); hold all;
%     plot(TimesI(k), r1_I-r1_B, 'r.')
%     
%     figure(8); hold all;
%     plot(TimesI(k), r2_I-r2_B, 'r.')

end

%%% Clearing Temporary Variables
clear x_I y_I rBC_I vBC_I vJC_I r1_I r2_I
clear x_B y_B rBC_B vBC_B vJC_B r1_B r2_B

% ------------------------------------------------------------------------
%%% Plots
% ------------------------------------------------------------------------
%%% Plotting differences in rH_ECEF components
figure
subplot(3,1,1)
plot(TimesI,rH_ECEF_I(:,1) - rH_ECEF_B(:,1))
title('rH\_ECEF\_I - rH\_ECEF\_B')
subplot(3,1,2)
plot(TimesI,rH_ECEF_I(:,2) - rH_ECEF_B(:,2))
subplot(3,1,3)
plot(TimesI,rH_ECEF_I(:,3) - rH_ECEF_B(:,3))

%%% Plotting Jacobian Constants
figure
plot(TimesI,JC_I)
title('JC\_I')

figure
plot(TimesB,JC_B)
title('JC\_B')

% %%% Filling in titles of "difference" plots
% figure(5); title('x_I - x_B')
% figure(6); title('y_I - y_B')
% figure(7); title('r1_I-r1_B')
% figure(8); title('r2_I-r2_B')

%%% Plotting differences in all components of vJC terms
figure; hold all
plot(TimesI,vs_I(:,1) - vs_B(:,1),'linewidth',1.5)
plot(TimesI,vs_I(:,2) - vs_B(:,2),'linewidth',1.5)
plot(TimesI,vs_I(:,3) - vs_B(:,3),'linewidth',1.5)
legend('vx diff','vy diff','vz diff')
title('vJC\_I-vJC\_B')

















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
function [ dY ] = EH_NumIntegrator_CR3BP_exp(t,Y,E_radius,uE,uJ,nE,E_a)
dY = zeros(6,1);

%%% Creating Europa Position (JCI)
ta = nE*t; % rads
rE = [E_a; 0; 0]; % km
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
function [ dY ] = EH_Body_NumIntegrator_CR3BP_exp(t,Y,E_radius,uE,uJ,nE,E_a)
dY = zeros(6,1);

%%% Unpack the Hopper state vector (ECEF)
yH = Y(1:3); % Hopper Position, km
dyH = Y(4:6); % Hopper Velocity, km/s

%%% Creating Europa Position (JCI)
ta = nE*t; % rads
rE_JCI = [E_a; 0; 0]; % km
rE_JCI = R3(rE_JCI,ta); % km

%%% Creating Hopper Position (ECI)
rH_ECI = R3(yH,ta); % km

%%% Creating Hopper Position (JCI)
rH_JCI = rH_ECI + rE_JCI; % km

%%% Determining Inertial Accelerations (JCI)
aH_JCI = (-uJ/(norm(rH_JCI)^3))*rH_JCI...
    + (-uE/(norm(rH_ECI)^3))*rH_ECI; % km/s^2

%%% Determining Acceleration of Europa (JCI)
aE = (-uJ/(norm(rE_JCI)^3))*rE_JCI;

%%% Determining Body Frame Hopper Acceleration (ECEF)
ddyH = aH_JCI - aE - 2*cross([0;0;nE],dyH) - cross([0;0;nE],cross([0;0;nE],yH)); % km/s^2

%%% Output the derivative of the state
dY(1:3) = dyH; % km/s
dY(4:6) = ddyH; % km/s^2
end


% ------------------------------------------------------------------------
%%% Inerital Frame Impact Event
% ------------------------------------------------------------------------
function [value, isterminal, direction] = impactEvent_I_exp(t,Y,E_radius,uE,uJ,nE,E_a)
%%% Creating Europa Position
ta = nE*t; % rads
rE = [E_a; 0; 0]; % km
rE = R3(rE,ta); % km

%%% Unpack the Hopper position vector (JCI)
yH = Y(1:3); % km

%%% Europa-Centric hopper position (ECI)
rH = yH - rE; % km

%%% Event function watching for when "value" = 0 (Hopper impacts Europa)
value = norm(rH) - E_radius;
isterminal = 1; % stops the integration
direction = -1; % negative direction only

end


% ------------------------------------------------------------------------
%%% Body Frame Impact Event
% ------------------------------------------------------------------------
function [value, isterminal, direction] = impactEvent_B_exp(t,Y,E_radius,uE,uJ,nE,E_a)
%%% Unpack the Hopper position vector (ECEF)
yH = Y(1:3);

%%% Event function watching for when "value" = 0 (Hopper impacts Europa)
value = norm(yH) - E_radius;
isterminal = 1; % stops the integration
direction = -1; % negative direction only

end






