% ------------------------------------------------------------------------
%%% Body Frame Impact Event
% ------------------------------------------------------------------------
function [value, isterminal, direction] = impactEvent_Body_CR3BP(t,Y,E_radius,uE,uJ,nE,E_a)
%%% Unpack the Hopper position vector (ECEF)
yH = Y(1:3);

%%% Event function watching for when "value" = 0 (Hopper impacts Europa)
value = norm(yH) - E_radius;
isterminal = 1; % stops the integration
direction = -1; % negative direction only

end