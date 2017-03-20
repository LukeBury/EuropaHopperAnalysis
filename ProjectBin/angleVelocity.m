%%% Inputs:
% 1) Magnitude of initial velocity
% 2) Azimuth of desired velocity (deg)
% 3) Elevation of desired velocity (deg)
% 4) Initial hopper latitude (deg)
% 5) Initial hopper longitude (deg)
%%% Outputs:
% 1) Initial velocity vector in body-fixed frame
function [ v0_bodyFixed ] = angleVelocity(v0, az, el, lat0, lon0)
fprintf('Make sure angleVelocity works for poles\n')
%%% Converting to radians
az = az*pi/180; % rad
el = el*pi/180; % rad
lat0 = lat0*pi/180; % rad
lon0 = lon0*pi/180; % rad

%%% Starting pointing straight South (Because of SEZ frame simplicity)
v0_sez = [v0, 0, 0]; 

%%% Rotating to elevation
v0_sez = R2(v0_sez, -el); 

%%% Rotating to azimuth
v0_sez = R3(v0_sez, -(az+pi));

%%% Rotating to Body Fixed frame
v0_bodyFixed = R3(R2(v0_sez, (pi/2-lat0)), lon0);
end

