function [dVelocity, dKE, T_impact, AngChange, Traveled, azColor, vH01, vH02, vH03] = EH_Mult_CR3BP(lat1, lon1, vmag)
% ------------------------------------------------------------------------
%%% Setting Initial Values
% ------------------------------------------------------------------------

%%% Time Constraints
ti = 0;
tf = 50000;
time = ti:.1:tf; 

%%% Jupiter Parameters
uJ = 126672520; % km^3 / s^2
J_pos = [0, 0, 0]; % km
J_radius = 69911; % km

%%% Europa Parameters
E_radius = 1560.8; % km
E_a = 671100; % km
uE = 3203.413216; % km^3 / s^2
nE = sqrt(uJ / E_a^3); % rad/s
% Circular Velocity
vE = sqrt(uJ/E_a); % km/s
wE = [0; 0; nE]; % Rot. velocity (tidally locked)

%%% Intial Europa State (Jupiter-Centric)
rE0 = [E_a, 0, 0]; % km
vE0 = [0, vE, 0]; % km

%%% Initial Hopper State (Jupiter-Centric) 
% Surface Position (latitude / longitude)
[rH01] = latlon2surfECEF(lat1, lon1, E_radius); % km

%%% Radial Velocity (Europa relative)
vH01 = (rH01/norm(rH01))*vmag;

%%% Correcting Velocity
vH02 = vH01 + cross(wE,rH01); % Adding rotational velocity component
vH03 = vH02 + vE0; % Making velocity relative to Jupiter
%%% Correcting Position
rH02 = rH01 + rE0; % Making position relative to Jupiter
%%% Setting Initial State Vector
X0 = [rH02 vH03]; % km km/s km km/s

% ------------------------------------------------------------------------
%%% Propagating the State with Numerically Integration
% ------------------------------------------------------------------------
%%% Setting integrator accuracy
tol = 1E-12;
options = odeset('Events',@impactEvent_CR3BP,'RelTol',tol,'AbsTol',tol);

%%% Propagating the STate
[Times,States] = ode45(@EH_NumIntegrator_CR3BP,time,X0,options,E_radius,uE,uJ,nE,E_a);
T_impact = Times(end);

% ------------------------------------------------------------------------
%%% Calculating Europa states
% ------------------------------------------------------------------------
r_Europa = zeros(size(Times,1),3);
for k = 1:size(Times,1)
    ta = nE * Times(k); % rads
    r_Europa(k,:) = R3(rE0,ta); % km
end

%%% Final States
rEf = r_Europa(end,:); % km
taf = nE * Times(end); % rad
vEf = R3(vE0,taf); % km/s

% ------------------------------------------------------------------------
%%% Calculating Hopper Initial and Final Energies in ECI
% ------------------------------------------------------------------------

v0_m = norm(States(1,4:6) - vE0)*1000; % m/s
vf_m = norm(States(end,4:6) - vEf)*1000; % m/s
dVelocity = vf_m - v0_m; % m/s
KE0 = .5*v0_m^2;
KEf = .5*vf_m^2;
dKE = KEf - KE0; % J/kg

% ------------------------------------------------------------------------
%%% Calculating Change in Hopper Surface Angle & Distance Traveled
% ------------------------------------------------------------------------
Pos0 = States(1,1:3) - rE0; % ECI / ECEF
Posf = States(end,1:3) - rEf; % ECI
% Rotating Final Position
theta = -nE * Times(end);
Posf = R3(Posf,theta); % ECEF

%%% Angle change and surface distance traveled
AngChange = atan2(norm(cross(Pos0,Posf)),dot(Pos0,Posf))*180/pi; % deg
Traveled = (AngChange*pi/180)*E_radius*1000; % m

%%% Finding direction traveled
[lat2, lon2] = ECEF2latlon(Posf); % rad
lat2 = lat2 * 180/pi; % deg
lon2 = lon2 * 180/pi; % deg
az = azimuth(lat1,lon1,lat2,lon2); % deg
[azColor] = azimuthColorWheel(az); % color [c1 c2 c3]
if lat1 == 90
    azColor = [0 0 1]; % Can only go south...
elseif lat1 == -90
    azColor = [1 0 0]; % Can only go north...
end

if Times(end) == tf
    warning off backtrace
    warning('No Impact in Simulation\n')
    fprintf('Lat1 = %3.2f',lat1')
    fprintf('Lon1 = %3.2f',lon1')
end

end
