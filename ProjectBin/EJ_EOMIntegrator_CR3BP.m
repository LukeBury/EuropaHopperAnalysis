% ------------------------------------------------------------------------
%%% Europa-and-Jupiter Numerical Integrator
% ------------------------------------------------------------------------
function [ dY ] = EJ_EOMIntegrator_CR3BP(t,Y,E_radius,uE,uJ,nE,E_a)
dY = zeros(6,1);

%%% Creating Europa Position (JCI)
ta = nE*t; % rads
rE = [E_a, 0, 0]; % km
rE = R3(rE,ta);

%%% Unpack the Hopper state vector (JCI)
yH = Y(1:3);
dyH = Y(4:6);

%%% Europa-Centric hopper position (ECI)
rH = yH - rE; % km

%%% Hopper Dynamics
ddyH = (-uJ/(norm(yH)^3))*yH + (-uE/(norm(rH)^3))*rH; % Jupiter and Europa pull km/s^2

%%% Output the derivative of the state
dY(1:3) = dyH;
dY(4:6) = ddyH;
end


