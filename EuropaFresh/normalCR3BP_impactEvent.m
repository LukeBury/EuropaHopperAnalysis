function [value, isterminal, direction] = normalCR3BP_impactEvent(t,X,u,rB1,rB2,rad2)
%%% Event function watching for when "value" = 0
r2 = sqrt((X(1)-rB2(1))^2 + (X(2)-rB2(2))^2 + (X(3)-rB2(3))^2);

value = r2 - rad2; % When the surface is impacted
isterminal = 1; % stops the integration
direction = -1; % negative direction only
end