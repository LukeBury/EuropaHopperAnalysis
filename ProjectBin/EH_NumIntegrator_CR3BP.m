% ------------------------------------------------------------------------
%%% Inertial Numerical Integrator
% ------------------------------------------------------------------------
function [ dY ] = EH_NumIntegrator_CR3BP(t,Y,RE,uE,uJ,nE,aE,E_theta0)
dY = zeros(6,1);

%%% Creating Europa Position (JCI)
ta = nE*t; % rads
rE_JCI = R3([aE; 0; 0],E_theta0); % km
rE_JCI = R3(rE_JCI,ta); % km

%%% Unpack the Hopper state vector (JCI)
yH = Y(1:3); % Hopper Position, km
dyH = Y(4:6); % Hopper Velocity, km/s

%%% Europa-Centric hopper position (ECI)
rH = yH - rE_JCI; % km

%%% Hopper Dynamics
fJ = (-uJ/(norm(yH)^3))*yH; % Jupiter Pull, km/s^2
fE = (-uE/(norm(rH)^3))*rH; % Europa Pull, km/s^2
ddyH = fJ + fE; % Jupiter and Europa pull km/s^2

%%% Output the derivative of the state
dY(1:3) = dyH; % km/s
dY(4:6) = ddyH; % km/s^2
end


