% ------------------------------------------------------------------------
%%% Body Frame Impact Event
% ------------------------------------------------------------------------
function [value, isterminal, direction] = impactEvent_Body_CR3BP(t,Y,RE,uE,uJ,nE,aE,E_theta0)
%%% Unpack the Hopper position vector (ECEF)
yH = Y(1:3);

%%% Event function watching for when "value" = 0 (Hopper impacts Europa)
value = norm(yH) - RE;
isterminal = 1; % stops the integration
direction = -1; % negative direction only

end