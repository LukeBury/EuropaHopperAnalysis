function [dX] = normalCR3BP_Int(t,X,u,rB1,rB2,rad2)
dX = zeros(6,1);

%%% Unpack the Hopper state vector
rP = X(1:3); % Hopper Position
vP = X(4:6); % Hopper Velocity

%%% Assigning Variables
x = rP(1); y = rP(2); z = rP(3);
dx = vP(1); dy = vP(2); dz = vP(3);

%%% Creating hopper distances to bodies
r1 = sqrt((rP(1)-rB1(1))^2 + (rP(2)-rB1(2))^2 + (rP(3)-rB1(3))^2);
r2 = sqrt((rP(1)-rB2(1))^2 + (rP(2)-rB2(2))^2 + (rP(3)-rB2(3))^2);

%%% Equations of Motion
ddx = 2*dy + x - (1-u)*(x+u)/(r1^3) - u*(x+u-1)/(r2^3);
ddy = -2*dx + y -((1-u)/(r1^3) + u/(r2^3))*y;
ddz = -((1-u)/(r1^3) + u/(r2^3))*z;

%%% Storing EOMs
ddxP = [ddx; ddy; ddz];

%%% Output the derivative of the state
dX(1:3) = vP; % km/s
dX(4:6) = ddxP; % km/s^2
end

