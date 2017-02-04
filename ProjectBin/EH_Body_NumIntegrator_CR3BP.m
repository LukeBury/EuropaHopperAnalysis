% ------------------------------------------------------------------------
%%% Body Frame Numerical Integrator
% ------------------------------------------------------------------------
function [ dY ] = EH_Body_NumIntegrator_CR3BP(t,Y,E_radius,uE,uJ,nE,E_a)
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
