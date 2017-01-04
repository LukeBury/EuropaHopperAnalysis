%%% Inputs
% 1) Azimuth (0-360 degrees)
%%% Outputs
% 1) Color code matrix [c1 c2 c3]
function [newColor] = azimuthColorWheel(az)
%%% Define Colors
r = [1 0 0]; % Red, North
g = [0 1 0]; % Green, East
y = [1 1 0]; % Yellow, West
b = [0 0 1]; % Blue, South

if (0 <= az) && (az < 90)
    newColor = ((az-0)/90)*(g-r) + r;
    
elseif (90 <= az) && (az < 180)
    newColor = ((az-90)/90)*(b-g) + g;
    
elseif (180 <= az) && (az < 270)
    newColor = ((az-180)/90)*(y-b) + b;
    
elseif (270 <= az) && (az < 360)
    newColor = ((az-270)/90)*(r-y) + y;
    
else
    warning('Bad Azimuth')
end

end