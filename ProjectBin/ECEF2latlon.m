%%% Given ECEF position vector, returns latitdue and longitude in RADIANS
function [lat, lon] = ECEF2latlon(r)
%%% Define axes
x = [1 0 0]';
z = [0 0 1];

%%% Latitude
lat = pi/2 - acos(dot(z,r)/norm(r));

%%% Longitude
lon = atan(r(2)/r(1));

%%% Quadrant corrections
% -x, +y
if r(1) < 0 && r(2) >= 0
    lon = lon + pi;
% -x, -y
elseif r(1) < 0 && r(2) < 0
    lon = lon - pi;
end
end