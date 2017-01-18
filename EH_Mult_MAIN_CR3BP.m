clear
clc
% close all
addpath('ProjectBin')
tic

%%% Plot results?
plotsOn = 1; % 0 = no, 1 = yes

%%% Set initial parameters
vmag = 0.015; % km/s, magnitude of velocity of radial hop
ddeg = 10; % deg, degree step size
latRange = [-90, 90]; % deg, range of latitudes to analyze 
lonRange = [-180; 180]; % deg, range of longitudes to analyze
%0.015
%%% Check compatibility of parameters
if rem((latRange(2)-latRange(1)),ddeg) ~= 0 || rem((lonRange(2)-lonRange(1)),ddeg) ~= 0
    warning off backtrace
    warning('delta-Degree value does not pair with latitude or longitude range')
    return
end

%%% Initialize matrices
np = size(latRange(1):ddeg:latRange(2),2) * size(lonRange(1):ddeg:lonRange(2),2);
results = zeros(np - 2*(size(lonRange(1):ddeg:lonRange(2),2)-1),6);

Energy = zeros(size(latRange(1):ddeg:latRange(2),2),size(lonRange(1):ddeg:lonRange(2),2)); % J/kg
dVels = zeros(size(latRange(1):ddeg:latRange(2),2),size(lonRange(1):ddeg:lonRange(2),2)); % m/s
TOF = zeros(size(latRange(1):ddeg:latRange(2),2),size(lonRange(1):ddeg:lonRange(2),2)); % sec
dAng = zeros(size(latRange(1):ddeg:latRange(2),2),size(lonRange(1):ddeg:lonRange(2),2)); % deg
dTraveled = zeros(size(latRange(1):ddeg:latRange(2),2),size(lonRange(1):ddeg:lonRange(2),2)); % m
azCols = zeros(size(latRange(1):ddeg:latRange(2),2),size(lonRange(1):ddeg:lonRange(2),2),3); % Color Codes [c1 c2 c3]

vH01s = zeros(size(latRange(1):ddeg:latRange(2),2),size(lonRange(1):ddeg:lonRange(2),2),3); % Velocity [vx vy vz]
vH02s = zeros(size(latRange(1):ddeg:latRange(2),2),size(lonRange(1):ddeg:lonRange(2),2),3); % Velocity [vx vy vz]
vH03s = zeros(size(latRange(1):ddeg:latRange(2),2),size(lonRange(1):ddeg:lonRange(2),2),3); % Velocity [vx vy vz]

%%% Loop through positions
kk = 1; % index
lat_i = 1; % index
comp = 1; % for assessing completeness of simulation
for rlat = latRange(1):ddeg:latRange(2)
    lon_i = 1; % index
    for rlon = lonRange(1):ddeg:lonRange(2)
        %%% Run Simulation
        [dVelocity, dEnergy, T_impact, AngChange, Traveled, azColor, vH01, vH02, vH03] =...
            EH_Mult_CR3BP(rlat, rlon, vmag);
        
        %%% Assign results
        results(kk,1) = rlat; % deg (input)
        results(kk,2) = rlon; % deg (input)
        results(kk,3) = vmag; % km/s (input)
        
        results(kk,4) = dVelocity; % m/s (output)
        results(kk,5) = dEnergy; % kJ/kg (output)
        results(kk,6) = T_impact; % sec (output)
        kk = kk + 1;
        
        Energy(lat_i,lon_i) = dEnergy; % J/kg
        TOF(lat_i,lon_i) = T_impact; % sec
        dVels(lat_i,lon_i) = dVelocity; % m/s
        dAng(lat_i,lon_i) = AngChange; % deg
        dTraveled(lat_i,lon_i) = Traveled; % m
        azCols(lat_i,lon_i,:) = azColor; % Color Codes [c1 c2 c3]
        vH02;
        vH01s(lat_i,lon_i,:) = vH01; % Hopper Velocity, km/s
        vH02s(lat_i,lon_i,:) = vH02; % Hopper Velocity, km/s
        vH03s(lat_i,lon_i,:) = vH03; % Hopper Velocity, km/s

        %%% Catching pole singularities
        if rlat == -90 || rlat == 90
            Energy(lat_i,:) = dEnergy;
            TOF(lat_i,:) = T_impact;
            dVels(lat_i,:) = dVelocity;
            dAng(lat_i,:) = AngChange;
            dTraveled(lat_i,:) = Traveled; % m
            breaklon_i = 1;
            for breaklon = lonRange(1):ddeg:lonRange(2)
                azCols(lat_i,breaklon_i,:) = azColor; % Color Codes [c1 c2 c3]
                breaklon_i = breaklon_i + 1;
            end
            break
        end
        lon_i = lon_i + 1;
    end
    lat_i = lat_i  + 1;
    clc
    fprintf('Complete: %3.2f%%\n',(comp/size(latRange(1):ddeg:latRange(2),2))*100)
    fprintf('Time:     %3.1f s\n',toc)
    fprintf('Latitude: %2.2f°\n',rlat)
    comp = comp + 1;
end
clear comp

if plotsOn == 1
    %%% Plotting delta-Energy Surface
    figure
    hold all
    surf(lonRange(1):ddeg:lonRange(2),latRange(1):ddeg:latRange(2),Energy)
    % surf(-180:ddeg:180,-90:ddeg:90,Energy)
    PlotBoi3('Longitude [°]','Latitude [°]','\DeltaSpecific Kinetic Energy [J/kg]',16);
    view(-15,10)

    %%% Plotting TOF Surface
    figure
    hold all
    surf(lonRange(1):ddeg:lonRange(2),latRange(1):ddeg:latRange(2),TOF)
    PlotBoi3('Longitude [°]','Latitude [°]','TOF [sec]',16);
    view(-15,10)

    %%% Plotting Changes in Surface Angle
    figure
    hold all
    surf(lonRange(1):ddeg:lonRange(2),latRange(1):ddeg:latRange(2),dAng)
    PlotBoi3('Longitude [°]','Latitude [°]','\Delta° Along Surface',16);
    view(-15,10)
    
    %%% Plotting Directions of Travel (3D plot)
    figure
    hold all
    lat_i = 1; % index
    E_radius = 1560.8; % km
    for rlat = latRange(1):ddeg:latRange(2)
        lon_i = 1; % index
        for rlon = lonRange(1):ddeg:lonRange(2)
            plot3(rlon,rlat,dAng(lat_i,lon_i)*(pi/180)*E_radius,'.','color',[azCols(lat_i,lon_i,1),azCols(lat_i,lon_i,2),azCols(lat_i,lon_i,3)],'markersize',40)
            lon_i = lon_i + 1;
        end
        lat_i = lat_i + 1;
    end
    PlotBoi3('Longitude [°]','Latitude [°]','Distance, km',16);
    dim = [.3 .6 .3 .3];
    str = 'Red = North | Green = East | Blue = South | Yellow = West';
    annotation('textbox',dim,'String',str,'FitBoxToText','on')
    view(-15,10)
    
    %%% Plotting Directions of Travel (Sphere)
    figure
    hold all
    lat_i = 1; % index
    for rlat = latRange(1):ddeg:latRange(2)
        lon_i = 1; % index
        for rlon = lonRange(1):ddeg:lonRange(2)
            [rSphere] = latlon2surfECEF(rlat, rlon, 1);
            rSphere = rSphere.*72;
            plot3(rSphere(1),rSphere(2),rSphere(3),'.','color',[azCols(lat_i,lon_i,1),azCols(lat_i,lon_i,2),azCols(lat_i,lon_i,3)],'markersize',50)
            lon_i = lon_i + 1;
        end
        lat_i = lat_i + 1;
    end
    quiver3(0,0,0,-150,0,0,...
        'linewidth',2,'color',[1 .5 0]); % Plotting Europa-Jupiter Vector
    axis equal
    PlotBoi3('X','Y','Z',16);
    
%     %%% Plotting directions of vH02 velocities (includes rotational effect)
%     figure
%     hold on
%     lon_i = 1;
%     Equator_i = (size(latRange(1):ddeg:latRange(2),2)+1)/2;
%     for rlon = lonRange(1):ddeg:lonRange(2)  
%         quiver3(0,0,rlon,vH02s(Equator_i,lon_i,1),vH02s(Equator_i,lon_i,2),vH02s(Equator_i,lon_i,3),'linewidth',2)
%         lon_i = lon_i + 1;
%     end
%     PlotBoi3('X','Y','Equatorial Longitude',16)
    
end


