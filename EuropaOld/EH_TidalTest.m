clear
clc
close all
addpath('ProjectBin')

% ------------------------------------------------------------------------
%%% Setting Initial Values
% ------------------------------------------------------------------------
%%% Jupiter Parameters
uJ = 126672520; % km^3 / s^2
J_pos = [0, 0, 0]; % km
J_radius = 69911; % km

%%% Europa Parameters
E_radius = 1560.8; % km
E_a = 671100; % km
uE = 3203.413216; % km^3 / s^2

%%% Intial Europa State (Jupiter-Centric)
rEJ = [E_a, 0, 0]; % km

%%% Creating hopper states around equator
figure
hold all
s = 500000000;
k = 1;
for lon = -180:5:180
    [rHE] = latlon2surfECEF(0, lon, E_radius); % km
    rHJ = rHE + rEJ;
    
    aHJ = (-uJ/(norm(rHJ)^3))*rHJ; % Hopper --> Jupiter, km/s^2
    aEJ = (-uJ/(norm(rEJ)^3))*rEJ; % Europa --> Jupiter, km/s^2
    
    aT = aHJ - aEJ;
    k = k + 1;
   
    quiver3(rHE(1),rHE(2),rHE(3),aT(1)*s,aT(2)*s,aT(3)*s, 'k','linewidth',2,'maxheadsize',10);
end
th = 0:.01:2*pi;
x = E_radius * cos(th);
y = E_radius * sin(th);
plot(x, y,'b','linewidth',2);
title('Tidal Acceleration Directions')
PlotBoi2('X','Y',16)



%%% Magnitude grows as altitude increases
figure
hold all
for i = 1:.01:1.2
    [rHE] = latlon2surfECEF(0, -45, E_radius); % km
    rHE = rHE*i;
    rHJ = rHE + rEJ;

    aHJ = (-uJ/(norm(rHJ)^3))*rHJ; % Hopper --> Jupiter, km/s^2
    aEJ = (-uJ/(norm(rEJ)^3))*rEJ; % Europa --> Jupiter, km/s^2

    aT = aHJ - aEJ;
    norm(aT)
    quiver3(rHE(1),rHE(2),rHE(3),aT(1)*s,aT(2)*s,aT(3)*s, 'k','linewidth',2,'maxheadsize',10);
end
title('Tidal Acceleration Magnitude Increasing w/ Altitude')
PlotBoi2('X','Y',16)