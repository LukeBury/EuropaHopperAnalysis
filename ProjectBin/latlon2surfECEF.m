%%% Given latitude (deg), longitude (deg), and body radius (km)
%%% Returns surface coordinates in ECEF
function [rECEF] = latlon2surfECEF(lat, lon, rad)
%%% Convert to radians
lat = lat * pi/180;
lon = lon * pi/180;

%%% Assign x coordinate
x = rad * cos(lat) * cos(lon);

%%% Assign y coordinate
y = rad * cos(lat) * sin(lon);

%%% Assign z coordinate
z = sin(lat) * rad;

rECEF = [x, y, z];
end