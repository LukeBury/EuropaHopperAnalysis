%%% bodySurface3: Plots a 3D spherical body
%%% Inputs: 
% 1) Radius of body
% 2) [x,y] position of body
% 3) [n,n,n] color assignment
function bodySurface2(body_rad, body_pos, color, linewidth)
x = body_rad * cos(0:.001:2*pi) + body_pos(1);
y = body_rad * sin(0:.001:2*pi) + body_pos(2);
plot(x, y, 'color', color, 'linewidth', linewidth);


end

