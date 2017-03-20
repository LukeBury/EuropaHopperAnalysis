%%% bodySurface3: Plots a 3D spherical body
%%% Inputs: 
% 1) Radius of body
% 2) [x,y,z] position of body
% 3) [n,n,n] color assignment
function bodySurface3(body_rad, body_pos, color)
[x,y,z] = sphere;
surf(x*body_rad+body_pos(1),y*body_rad+body_pos(2),z*body_rad+body_pos(3))
colormap(color)
view(3);
end