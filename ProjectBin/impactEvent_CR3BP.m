% ------------------------------------------------------------------------
%%% Inertial Frame Impact Event
% ------------------------------------------------------------------------
function [value, isterminal, direction] = impactEvent_CR3BP(t,Y,E_radius,uE,uJ,nE,E_a)
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