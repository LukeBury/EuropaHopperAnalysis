%%% bodySurface3: Plots a 3D spherical body
%%% Inputs: 
% 1) Radius of body
% 2) [x,y,z] position of body
% 3) [n,n,n] color assignment
function bodySurface3(body_rad, Epos, color)
[xE,yE,zE] = sphere;
surf(xE*body_rad+Epos(1),yE*body_rad+Epos(2),zE*body_rad+Epos(3))
colormap(color)
view(3);
end