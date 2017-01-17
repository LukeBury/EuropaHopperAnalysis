% ------------------------------------------------------------------------
%%% Europa-and-Jupiter Numerical Integrator
% ------------------------------------------------------------------------
function [ dY ] = EH_NumIntegrator_CR3BP(t,Y,E_radius,uE,uJ,nE,E_a)
dY = zeros(6,1);

%%% Creating Europa Position (JCI)
ta = nE*t; % rads
rE = [E_a; 0; 0]; % km
rE = R3(rE,ta);

%%% Unpack the Hopper state vector (JCI)
yH = Y(1:3);
dyH = Y(4:6);

%%% Europa-Centric hopper position (ECI)
rH = yH - rE; % km

%%% Hopper Dynamics
fJ = (-uJ/(norm(yH)^3))*yH;
fE = (-uE/(norm(rH)^3))*rH;
ddyH = fJ + fE; % Jupiter and Europa pull km/s^2

%%% Output the derivative of the state
dY(1:3) = dyH;
dY(4:6) = ddyH;
end


