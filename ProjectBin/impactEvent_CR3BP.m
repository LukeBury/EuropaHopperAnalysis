% ------------------------------------------------------------------------
%%% Europa-and-Jupiter Numerical Integrator
% ------------------------------------------------------------------------
function [value, isterminal, direction] = impactEvent(t,Y,E_radius,uE,uJ,nE,E_a)
%%% Creating Europa Position
ta = nE*t; % rads
rE = [E_a; 0; 0]; % km
rE = R3(rE,ta); % km

%%% Unpack the Hopper position vector
yH = Y(1:3);

%%% Europa-Centric hopper position
rH = yH - rE; % km

%%% Event function watching for when "value" = 0 (Hopper impacts Europa)
value = norm(rH) - E_radius;
isterminal = 1; % stops the integration
direction = -1; % negative direction only

end