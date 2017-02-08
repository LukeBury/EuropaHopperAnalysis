% ------------------------------------------------------------------------
%%% Inertial Frame Impact Event
% ------------------------------------------------------------------------
function [value, isterminal, direction] = impactEvent_CR3BP(t,Y,RE,uE,uJ,nE,aE,E_theta0)
%%% Creating Europa Position
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