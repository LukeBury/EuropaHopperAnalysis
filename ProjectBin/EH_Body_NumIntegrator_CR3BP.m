% ------------------------------------------------------------------------
%%% Europa-and-Jupiter Numerical Integrator
% ------------------------------------------------------------------------
function [ dY ] = EH_NumIntegrator_CR3BP(t,Y,E_radius,uE,uJ,nE,E_a)
dY = zeros(6,1);
% 
% %%% Creating Europa Position (JCI)
% ta = nE*t; % rads
% rE = [E_a; 0; 0]; % km
% rE = R3(rE,ta); % km
% 
% %%% Unpack the Hopper state vector (JCI)
% yH = Y(1:3); % Hopper Position, km
% dyH = Y(4:6); % Hopper Velocity, km/s
% 
% %%% Europa-Centric hopper position (ECI)
% rH = yH - rE; % km
% 
% %%% Hopper Dynamics
% fJ = (-uJ/(norm(yH)^3))*yH; % Jupiter Pull, km/s^2
% fE = (-uE/(norm(rH)^3))*rH; % Europa Pull, km/s^2
% ddyH = fJ + fE; % Jupiter and Europa pull km/s^2

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


