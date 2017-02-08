clear
clc
close all
addpath('ProjectBin')

% ------------------------------------------------------------------------
%%% Plot Switches
% ------------------------------------------------------------------------
% Plot Europa ECEF?
ECEFplot = 0; % no = 0, yes = 1
scale1 = 8; % Plot Scale (limits = scale1 x E_radius)

% Plot Europa ECI?
ECIplot = 0; % no = 0, yes = 1

% Plot Jupiter System Inertial?
JCIplot = 0; % no = 0, yes = 1

% Plot ECI distance vs Time?
rECIplot = 0; % no = 0, yes = 1

% Plot East-ness vs Time?
EquatorialEastPlot = 0; % no = 0, yes = 1

% Plot Numerical vs Analytical Results?
compareAnalyticalNumericalPlots = 0; % no = 0, yes = 1

% Run movie?
runECEFMovie = 0; % no = 0, yes = 1
runECIMovie = 0; % no = 0, yes = 1
scale2 = 1.4;
framespeed = 1; % Higher is faster
% ------------------------------------------------------------------------
%%% Setting Initial Values
% ------------------------------------------------------------------------
%%% Time Constraints
ti = 0;
tf = 520000; % sec
dt = .1;
time = ti:dt:tf;

%%% Jupiter Parameters
uJ = 126672520; % km^3 / s^2
rJ0 = [0, 0, 0]; % km
RJ = 69911; % km

%%% Europa Parameters
RE = 1560.8; % km
aE = 671100; % km
uE = 3203.413216; % km^3 / s^2
nE = sqrt(uJ / (aE^3)); % rad/s
% Circular Velocity
vE = sqrt(uJ/aE); % km/s
wE = [0; 0; nE]; % Rot. velocity (tidally locked)

%%% Intial Europa State
E_theta0 = 0*(pi/180); % Initial position of Europa about Jupiter from +x, rads
rE0_JCI = R3([aE, 0, 0],E_theta0); % km
vE0_JCI = R3([0, vE, 0],E_theta0); % km/s

%%% Initial Hopper State
% Surface Position (latitude / longitude)
lat1 = 0; % deg (-90:90)
lon1 = -45; % deg (-180:180)
[rH0_ECEF] = latlon2surfECEF(lat1, lon1, RE); % km
rH0_ECI = R3(rH0_ECEF,E_theta0); % km
rH0_JCI = rH0_ECI + rE0_JCI; % km

%%% Radial Velocity (Europa relative)
v_mag = .015; % km/s
vH0_ECEF = (rH0_ECEF/norm(rH0_ECEF))*v_mag;
% 
% %%% Vector Velocity (Europa relative)
% vH0_ECEF = [.2, -vE/7.9, .8]; % km/s (-45,0)
% % vH0_ECEF = [-1.7, -.49, .6]; % km/s (-180,0)
% % vH0_ECEF = [-1.7, -.2, .6]; % km/s (-180,0)
% % vH0_ECEF = [.2, -1.6, 1]; % km/s (-180,0)

%%% Creating ECI and JCI Initial Hopper Velocity
vH0_ECI = R3(vH0_ECEF,E_theta0) + cross(wE,rH0_ECI); % km/s
vH0_JCI = vH0_ECI + vE0_JCI; % km/s

% ------------------------------------------------------------------------
%%% Propagating the State with Numerical Integration
% ------------------------------------------------------------------------
%%% Setting integrator options
tol = 1E-13;
optionsI = odeset('Events',@impactEvent_CR3BP,'RelTol',tol,'AbsTol',tol);

%%% Setting Initial State Vector (JCI)
X0_JCI = [rH0_JCI vH0_JCI]; % km km/s

%%% Propagating the State
[Times,StatesI] = ode45(@EH_NumIntegrator_CR3BP,time,X0_JCI,optionsI,RE,uE,uJ,nE,aE,E_theta0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Post-Integration Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------------------------------------------------
%%% Calculating Europa states
% ------------------------------------------------------------------------
rE_JCI = zeros(length(Times),3);
vE_JCI = zeros(length(Times),3);
rotAngles = zeros(length(Times),1);
for k = 1:length(Times)
    rotAngles(k) = nE * Times(k); % rads
    rE_JCI(k,:) = R3(rE0_JCI,rotAngles(k)); % km, JCI Europa position vectors
    vE_JCI(k,:) = R3(vE0_JCI,rotAngles(k)); % km/s, JCI Europa velocity vectors
end
%%% Assigning final states of Europa (for simplicity)
rEf_JCI = rE_JCI(end,:); % km
vEf_JCI = R3(vE0_JCI,rotAngles(end)); % km/s

% ------------------------------------------------------------------------
%%% Calculating rECI positions
% ------------------------------------------------------------------------
rH_ECI = StatesI(:,1:3) - rE_JCI;

% ------------------------------------------------------------------------
%%% Calculating rECEF positions
% ------------------------------------------------------------------------
rH_ECEF = zeros(length(Times),3);
for k = 1:length(Times)
    rH_ECEF(k,:) = R3(rH_ECI(k,:),-rotAngles(k));
end

% ------------------------------------------------------------------------
%%% Calculating Hopper Initial and Final Energies in ECI
% ------------------------------------------------------------------------
v0_m = norm(vH0_ECI)*1000; % m/s, ECI
vf_m = norm(StatesI(end,4:6) - vEf_JCI)*1000; % m/s, ECI
dVelocity = vf_m - v0_m; % m/s
KE0 = .5*v0_m^2;
KEf = .5*vf_m^2;
dKE = KEf - KE0; % J/kg .... m^2/s^2

% ------------------------------------------------------------------------
%%% Calculating Change in Hopper Surface Angle & Distance Traveled
% ------------------------------------------------------------------------
%%% Angle change and surface distance traveled in ECEF
AngChange = atan2(norm(cross(rH_ECEF(1,:),rH_ECEF(end,:))),dot(rH_ECEF(1,:),rH_ECEF(end,:)))*180/pi; % deg
Traveled = (AngChange*pi/180)*RE*1000; % m

%%% Finding direction traveled
[lat2, lon2] = ECEF2latlon(rH_ECEF(end,:)); % rad
lat2 = lat2 * 180/pi; % deg
lon2 = lon2 * 180/pi; % deg
az = azimuth(lat1,lon1,lat2,lon2); % deg
[azColor] = azimuthColorWheel(az); % color [c1 c2 c3]
if lat1 == 90
    azColor = [0 0 1]; % Can only go south...
elseif lat1 == -90
    azColor = [1 0 0]; % Can only go north...
end

% ------------------------------------------------------------------------
%%% Calculating ECEF Hopper Accelerations
% ------------------------------------------------------------------------
%%% Initializing ECI and ECEF acceleration matrices
aH_JCI = zeros(length(Times),3); % km/s^2
aH_ECEF_I = zeros(length(Times),3); % km/s^2

vH_B_tempI = zeros(length(Times),3);
t4_tempI = zeros(length(Times),3);
t3_tempI = zeros(length(Times),3);
aH_JCI_tempI = zeros(length(Times),3);
aE_tempI = zeros(length(Times),3);

%%% Calculating ECI and ECEF Hopper accelerations at each time step
for k = 1:length(Times)        
    %%% Determining Inertial Accelerations in Body Frame
    aH_JCI(k,:) = (-uJ/(norm(StatesI(k,1:3))^3))*StatesI(k,1:3)...
        + (-uE/(norm(rH_ECI(k,:))^3))*rH_ECI(k,:); % km/s^2
    aH_JCI_tempI(k,:) = aH_JCI(k,:);
    
    %%% Determining Inertial Velocities in Body Frame
    vH_B = StatesI(k,4:6) - vE_JCI(k,:) - cross(wE,rH_ECEF(k,:)); % km/s
    vH_B_tempI(k,:) = vH_B;
    
    %%% Determining Acceleration of Europa in JCI
    aE_tempI = (-uJ/(norm(rE_JCI(k,:))^3))*rE_JCI(k,:);
    
    t3_tempI = 2*cross(wE',vH_B);
    t4_tempI = cross(wE',cross(wE',rH_ECEF(k,:)));
    %%% Determining Body Frame Acceleration
    aH_ECEF_I(k,:) = aH_JCI(k,:) - aE_tempI - 2*cross(wE',vH_B) - cross(wE',cross(wE',rH_ECEF(k,:))); % km/s^2
    
end
clear vH_B
fprintf('Clear some of these vars after testing\n')
vAnalytical = cumtrapz(Times, aH_ECEF_I(:,:)) + vH0_ECEF;
rAnalytical = cumtrapz(Times, vAnalytical(:,:)) + rH0_ECEF;

% ------------------------------------------------------------------------
%%% Calculating Initial Tidal Accelerations
% ------------------------------------------------------------------------
%%% Finding Acceleration Vectors
aHJ0_JCI = (-uJ/(norm(rH0_JCI)^3))*rH0_JCI; % Hopper --> Jupiter, km/s^2
aEJ0_JCI = (-uJ/(norm(rE0_JCI)^3))*rE0_JCI; % Europa --> Jupiter, km/s^2

% %%% Tidal Acceleration
% aT0 = aHJ0_JCI - aEJ0_JCI; % km/s^2
% 
% %%% Theoretical Distance Traveled for Tangential Tidal Acceleration
% d_TanTidal = .5*norm(aT0)*Times(end)*Times(end)*1000;

% ------------------------------------------------------------------------
%%% East Calculations
% ------------------------------------------------------------------------
% relPos= zeros(length(Times),3);
% EastAnalyticalPos = zeros(length(Times),1);
% EastAnalyticalVel = zeros(length(Times),1);
% EastAnalyticalAcc = zeros(length(Times),1);
% aTs = zeros(length(Times),3);
% EastAT = zeros(length(Times),1);
% EastUVec = zeros(length(Times),3);
% 
% liftoffUVec = (rH_ECEF(1,:)./norm(rH_ECEF(1,:))); % Uvec to liftoff spot in ECEF
% EastUVec(1,:) = cross(liftoffUVec,[0,0,-1]); % unit vector pointing local east from liftoff spot, km
% 
% for k = 1:length(Times)
%     %%% Finding Eastern Unit Vector
% %     EastUVec(k,:) = cross(rECEF_Hopper(k,:)./norm(rECEF_Hopper(k,:)), [0,0,-1]);
%     EastUVec(k,:) = cross(rH_ECEF(1,:)./norm(rH_ECEF(1,:)), [0,0,-1]);
%     
%     %%% Numerical East
%     th = nE*Times(k); % How far Europa has rotated, rad
%     relPos(k,:) = rH_ECEF(k,:) - rH_ECEF(1,:); % relative position of hopper to starting point, km
%     
%     %%% Analytical East
%     EastAnalyticalPos(k) = dot(EastUVec(k,:),rAnalytical(k,:));
%     EastAnalyticalVel(k) = dot(EastUVec(k,:),vAnalytical(k,:));
%     EastAnalyticalAcc(k) = dot(EastUVec(k,:),aH_ECEF(k,:));
%     
%     %%% Tidal East
%     aHJ0_JCI = (-uJ/(norm(States(k,1:3))^3))*States(k,1:3); % Hopper --> Jupiter, km/s^2
%     aEJ0_JCI = (-uJ/(norm(rE_JCI(k,:))^3))*rE_JCI(k,:); % Europa --> Jupiter, km/s^2
%     aTs(k,:) = R3(aHJ0_JCI - aEJ0_JCI,-th); % Tidal Accelerations in ECI frame
%     EastAT(k) = dot(EastUVec(k,:),aTs(k,:));
%     
% end
% 
% %%% Finding Numerical Acceleration (from integration)
% relVel = [diff(rH_ECEF(:,1)) diff(rH_ECEF(:,2)) diff(rH_ECEF(:,3))].*(1/dt); % km/s, ECEF
% % relVel = [diff(relPos(:,1)) diff(relPos(:,2)) diff(relPos(:,3))].*(1/dt); % km/s, ECEF
% relAcc = [diff(relVel(:,1)) diff(relVel(:,2)) diff(relVel(:,3))].*(1/dt); % km/s^2, ECEF
% 
% EastNumericalPos = zeros(length(Times),1);
% EastNumericalVel = zeros(length(Times)-1,1);
% EastNumericalAcc = zeros(length(Times)-2,1);
% 
% for k = 1:length(Times)
%     EastNumericalPos(k) = dot(EastUVec(k,:),relPos(k,:)); % component of relative position in the East direction
%     if k < length(Times)
%         EastNumericalVel(k) = dot(EastUVec(k,:),relVel(k,:));
%     end
%     if k < (length(Times)-1)
%         EastNumericalAcc(k) = dot(EastUVec(k,:),relAcc(k,:));
%     end
% end

% ------------------------------------------------------------------------
%%% Jacobi Constant (Numerical)
% ------------------------------------------------------------------------
JC_N = zeros(length(Times),1);
% rBC = rE0.*(uE/uJ);
mu = uE/(uE + uJ);
P = 2*pi*sqrt((aE^3)/uJ);
vNorm = P/(aE*2*pi);

for k = 1:length(Times)
    th = Times(k)*nE; % rads
    
    %%% Normalized method
%     rBC = R3(rE0.*(uE/uJ), th); % Barycenter position (JCI)
%     vBC = cross(wE',rBC); % Barycenter velocity (JCI) 
%     
%     x = (rECEF_Hopper(k,1)+rE0(1)-rBC(1))./rE0(1); % Hopper x-position in rotating barycentric frame
%     y = (rECEF_Hopper(k,2))./rE0(1); % Hopper y-position in rotating barycentric frame
%     r1 = norm(States(k,1:3))./rE0(1); % Distance from hopper to Jupiter
%     r2 = norm(rECEF_Hopper(k,:))./rE0(1); % Distance from hopper to Europa
%     r = [x, y, 0]; % Hopper position in rotating barycentric frame
%     vJC = (States(k,4:6) - vBC).*vNorm - cross([0,0,1],r); % Hopper velocity in rotating barycentric frame
% 
%     JC(k) = (x)^2 + (y)^2 + 2*(1-mu)/r1 + 2*mu/r2 - norm(vJC)^2; % Jacobian constant
    
    %%% Non-Normalized Method
    rBC = R3(rE0_JCI.*(uE/uJ), th); % Barycenter position (JCI)
    vBC = cross(wE',rBC); % Barycenter velocity (JCI) 
    x = rH_ECEF(k,1)+rE0_JCI(1)-rBC(1);
    y = rH_ECEF(k,2);
    r = [x, y, 0];
%     r = R3(States(k,1:3)-rBC,-th)
%     x = r(1);
%     y = r(2);
    r1 = norm(StatesI(k,1:3));
    r2 = norm(rH_ECEF(k,:));
    
    vJC = StatesI(k,4:6) - vBC - cross(wE',r);
    
    JC_N(k) = (nE^2)*(x^2 + y^2) + 2*(uJ/r1 + uE/r2) - norm(vJC)^2;

end


clear x y rBC vBC vJC mu r r1 r2 P 

% figure
% plot(Times,JC_N)
% PlotBoi2('Times, sec','Jacobian Constant',16)
% figure
% plot(Times,(JC_N-JC_N(1))*100./(JC_N(1)))
% PlotBoi2('Times, sec','Jacobian Constant %Change',16)


% ------------------------------------------------------------------------
%%% Running Numerical Integration of Analytical Body Acceleration
% ------------------------------------------------------------------------
%%% Setting integrator options
optionsB = odeset('Events',@impactEvent_Body_CR3BP,'RelTol',tol,'AbsTol',tol);

%%% Setting Initial State Vector (ECEF)
X0_ECEF = [rH0_ECEF vH0_ECEF]; % km, km/s

%%% Propagating the State
[TimesB,StatesB] = ode45(@EH_Body_NumIntegrator_CR3BP,time,X0_ECEF,optionsB,RE,uE,uJ,nE,aE,E_theta0);

% rE_JCI_tempB = zeros(length(Times),3);
% rH_ECI_tempB = zeros(length(Times),3);
% rH_JCI_tempB = zeros(length(Times),3);
t3_tempB = zeros(length(Times),3);
t4_tempB = zeros(length(Times),3);
aH_JCI_tempB = zeros(length(Times),3);
aE_tempB = zeros(length(Times),3);

aH_ECEF_B = zeros(length(Times),3);
for k = 1:length(TimesB)
    %%% Unpack the Hopper state vector (ECEF)
    rH_ECEF_tempB = StatesB(k,1:3); % Hopper Position, km
    vH_ECEF_tempB = StatesB(k,4:6); % Hopper Velocity, km/s
    
    %%% Creating Europa Position (JCI)
    rE_JCI_tempB = [aE,0 0]; % km
    rE_JCI_tempB = R3(rE_JCI_tempB,rotAngles(k)); % km
    
    %%% Creating Hopper Position (ECI)
    rH_ECI_tempB = R3(rH_ECEF_tempB,rotAngles(k)); % km
    
    %%% Creating Hopper Position (JCI)
    rH_JCI_tempB = rH_ECI_tempB + rE_JCI_tempB; % km
    
    %%% Determining Inertial Accelerations (JCI)
    aH_JCI_tempB(k,:) = (-uJ/(norm(rH_JCI_tempB)^3))*rH_JCI_tempB...
        + (-uE/(norm(rH_ECI_tempB)^3))*rH_ECI_tempB; % km/s^2
    
    %%% Determining Acceleration of Europa (JCI)
    aE_tempB(k,:) = (-uJ/(norm(rE_JCI_tempB)^3))*rE_JCI_tempB;
    
    t3_tempB(k,:) = 2*cross(wE',vH_ECEF_tempB);
    t4_tempB(k,:) = cross(wE',cross(wE',rH_ECEF_tempB));
    %%% Storing Body-Frame Acceleration (ECEF)
    aH_ECEF_B(k,:) = aH_JCI_tempB(k,:) - aE_tempB(k,:) - 2*cross(wE',vH_ECEF_tempB) - cross(wE',cross(wE',rH_ECEF_tempB)); % km/s^2
end
clear rH_ECEF_tempB vH_ECEF_tempb rE_JCI_temp rH_ECI_temp rH_JCI_temp aH_JCI_temp aE_temp

figure
subplot(3,1,1)
plot(Times,aH_JCI_tempI(:,1)-aH_JCI_tempB(:,1))
title('aH-JCI-temp Diff')
subplot(3,1,2)
plot(Times,aH_JCI_tempI(:,2)-aH_JCI_tempB(:,2))
subplot(3,1,3)
plot(Times,aH_JCI_tempI(:,3)-aH_JCI_tempB(:,3))

figure
subplot(3,1,1)
plot(Times,aE_tempI(:,1)-aE_tempB(:,1))
title('aE-temp Diff')
subplot(3,1,2)
plot(Times,aE_tempI(:,2)-aE_tempB(:,2))
subplot(3,1,3)
plot(Times,aE_tempI(:,3)-aE_tempB(:,3))

figure
subplot(3,1,1)
plot(Times,t3_tempI(:,1)-t3_tempB(:,1))
title('t3 Diff')
subplot(3,1,2)
plot(Times,t3_tempI(:,2)-t3_tempB(:,2))
subplot(3,1,3)
plot(Times,t3_tempI(:,3)-t3_tempB(:,3))

figure
subplot(3,1,1)
plot(Times,t4_tempI(:,1)-t4_tempB(:,1))
title('t4 Diff')
subplot(3,1,2)
plot(Times,t4_tempI(:,2)-t4_tempB(:,2))
subplot(3,1,3)
plot(Times,t4_tempI(:,3)-t4_tempB(:,3))

figure
subplot(3,1,1)
plot(Times,vH_B_tempI(:,1) - StatesB(:,4))
title('vH-B Diff')
subplot(3,1,2)
plot(Times,vH_B_tempI(:,2) - StatesB(:,5))
subplot(3,1,3)
plot(Times,vH_B_tempI(:,3) - StatesB(:,6))

% figure
% subplot(3,1,1)
% plot(Times,rH_ECEF(:,1) - StatesB(:,1))
% title('Pos Diff')
% subplot(3,1,2)
% plot(Times,rH_ECEF(:,2) - StatesB(:,2))
% subplot(3,1,3)
% plot(Times,rH_ECEF(:,3) - StatesB(:,3))
% 
% figure
% subplot(3,1,1)
% plot(Times,aH_ECEF_I(:,1) - aH_ECEF_B(:,1))
% title('Acc Diff')
% subplot(3,1,2)
% plot(Times,aH_ECEF_I(:,2) - aH_ECEF_B(:,2))
% subplot(3,1,3)
% plot(Times,aH_ECEF_I(:,3) - aH_ECEF_B(:,3))


% ------------------------------------------------------------------------
%%% Jacobi Constant (Analytical)
% ------------------------------------------------------------------------
JC_B = zeros(length(Times),1);
mu = uE/(uE + uJ);
P = 2*pi*sqrt((aE^3)/uJ);
vNorm = P/(aE*2*pi);

for k = 1:length(Times)
    th = Times(k)*nE; % rads
    
%     %%% Normalized method
%     rBC = R3(rE0.*(uE/uJ), th); % Barycenter position (JCI)
%     vBC = cross(wE',rBC); % Barycenter velocity (JCI) 
%     
%     x = (rECEF_Hopper(k,1)+rE0(1)-rBC(1))./rE0(1); % Hopper x-position in rotating barycentric frame
%     y = (rECEF_Hopper(k,2))./rE0(1); % Hopper y-position in rotating barycentric frame
%     r1 = norm(States(k,1:3))./rE0(1); % Distance from hopper to Jupiter
%     r2 = norm(rECEF_Hopper(k,:))./rE0(1); % Distance from hopper to Europa
%     r = [x, y, 0]; % Hopper position in rotating barycentric frame
%     vJC = (States(k,4:6) - vBC).*vNorm - cross([0,0,1],r); % Hopper velocity in rotating barycentric frame
% 
%     JC(k) = (x)^2 + (y)^2 + 2*(1-mu)/r1 + 2*mu/r2 - norm(vJC)^2; % Jacobian constant
    
    %%% Non-Normalized Method
    rBC = R3(rE0_JCI.*(uE/uJ), th); % Barycenter position (JCI)
    vBC = cross(wE',rBC); % Barycenter velocity (JCI) 
    x = StatesB(k,1)+rE0_JCI(1)-rBC(1);
    y = StatesB(k,2);
    r1 = norm(StatesB(k,1:3)+rE0_JCI);
    r2 = norm(StatesB(k,1:3));
    
    vJC = StatesB(k,4:6);
    
    JC_B(k) = (nE^2)*(x^2 + y^2) + 2*(uJ/r1 + uE/r2) - norm(vJC)^2;
     
%     figure(3)
%     hold all
%     plot(Times(k),(nE^2)*(x^2),'ro','markersize',6)
%     title('x')
%     
%     figure(4)
%     hold all
%     plot(Times(k),(nE^2)*(y^2),'ro','markersize',6)
%     title('y')
%     
%     figure(5)
%     hold all
%     plot(Times(k),2*uJ/r1,'ro','markersize',6)
%     title('r1')
%     
%     figure(6)
%     hold all
%     plot(Times(k),2*uE/r2,'ro','markersize',6)
%     title('r2')
%     
%     figure(7)
%     hold all
%     plot(Times(k),-norm(vJC)^2,'ro','markersize',6)
%     title('vJC')
end
clear x y rBC vBC vJC mu r r1 r2 P 

figure
plot(Times,JC_B,'linewidth',2)
PlotBoi2('Times, sec','Jacobian Constant',16)
title('Analytical JC')
figure
plot(Times,(JC_B-JC_B(1))*100./(JC_B(1)),'linewidth',2)
PlotBoi2('Times, sec','Jacobian Constant %Change',16)
title('Analytical JC %Change')


%%% Attempt at plotting contour lines
% x = linspace(-2*J_radius, 1.2*E_a, 1000);
% y = linspace(-2*E_radius, 2*E_radius, 1000);
% [X, Y] = meshgrid(x,y);
% JC_Test = 5.703572862448475e+02;
% rBC_Test = rE0.*(uE/uJ);
% 
% r1 = norm(StatesB(k,1:3)+rE0);
% r2 = norm(StatesB(k,1:3));
% 
% 
% Z = (nE^2).*(X.^2 + Y.^2) + 2*(uJ./norm([X,Y,0]+rBC_Test) + uE./norm([X,Y,0]-rE0));
% % (nE^2)*(x^2 + y^2) + 2*(uJ/r1 + uE/r2) - norm(vJC)^2;
% 
% % A_Test = JC_Test/(nE^2) + (2/(nE^2))*(uJ/norm(States(1,1:3)) + uE/norm(StatesB(1,1:3)));
% 
% 
% figure
% contour(X,Y,Z)

rBC_Test = rE0_JCI.*(uE/uJ);
x = linspace(aE - 1.3*RE, aE + 1.3*RE, 100);
y = linspace(-1.3*RE, 1.3*RE, 100);
z = zeros(length(x),length(x));
r1 = zeros(length(x),length(x));
r2 = zeros(length(x),length(x));
for xk = 1:length(x)
    for yk = 1:length(x)
        r1(xk,yk) = norm([x(xk), y(yk), 0] + rBC_Test);
        r2(xk,yk) = norm([x(xk), y(yk), 0] - rE0_JCI + rBC_Test);
        z(xk,yk) = (nE^2)*(x(xk)^2 + y(yk)^2) + 2*(uJ/r1(xk,yk) + uE/r2(xk,yk));
    end
end

% figure
% hold all
% % surf(x,y,z)
% contour(x,y,z,[567 568 569 570 570.35728 571 874],'ShowText','on')
% PlotBoi3('X','Y','Z',10)
% %%% Plotting Europa equator
% th = 0:.01:2*pi;
% xE = E_radius * cos(th);
% yE = E_radius * sin(th);
% plot(xE + E_a, yE,'b','linewidth',2);





































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------------------------------------------------
%%% Plotting the Europa-Centered-Europa-Fixed Frame (ECEF)
% ------------------------------------------------------------------------
if ECEFplot == 1
    figure
    hold all
    axis equal
    trackWidth = 2;
    lineWidth = 3;
    PlotBoi3('X, km','Y, km','Z, km',16)

    %%% Plotting Hopper motion
    plot3(rH_ECEF(:,1),rH_ECEF(:,2),rH_ECEF(:,3),'m','linewidth',trackWidth);
    
    %%% Plotting Europa equator
    th = 0:.01:2*pi;
    x = RE * cos(th);
    y = RE * sin(th);
    plot(x, y,'b','linewidth',lineWidth);

    %%% Coloring in Europa (2D)
    n = 5000;
    THETA=linspace(0,2*pi,n);
    RHO=ones(1,n)*(RE);
    [X,Y] = pol2cart(THETA,RHO);
    fill(X, Y, 'c');

    %%% Plotting Europa Surface (3D)
    bodySurface3(RE, [0 0 0], [0 1 1]);

    %%% Plotting Europa-Jupiter Vector
    quiver3(0,0,0,-2*RE,0,0,...
        'linewidth',2,'color',[1 .5 0]);

    %%% Plotting Intial Velocity
    sc1 = 1000; % scalar
    quiver3(rH0_ECEF(1),rH0_ECEF(2),rH0_ECEF(3),...
        sc1*(vH0_ECEF(1)),sc1*(vH0_ECEF(2)),sc1*(vH0_ECEF(3)),...
        'linewidth',2,'color',[0 0 0])

    xlim([-RE*scale1 RE*scale1])
    ylim([-RE*scale1 RE*scale1])
    zlim([-RE*scale1 RE*scale1])
    
    %%% Frame
    title('ECEF')
end

% ------------------------------------------------------------------------
%%% Plotting the Europa Inertial Frame (ECI)
% ------------------------------------------------------------------------
if ECIplot == 1
    figure
    hold all
    axis equal
    trackWidth = 2;
    lineWidth = 3;
    PlotBoi3('X, km','Y, km','Z, km',16)

    %%% Plotting Hopper motion
    plot3(rH_ECI(:,1),rH_ECI(:,2),rH_ECI(:,3),'m','linewidth',trackWidth);

    %%% Plotting Europa equator
    th = 0:.01:2*pi;
    x = RE * cos(th) + rE0_JCI(1);
    y = RE * sin(th) + rE0_JCI(2);
    plot(x, y,'b','linewidth',lineWidth);

    %%% Coloring in Europa (2D)
    n = 5000;
    THETA=linspace(0,2*pi,n);
    RHO=ones(1,n)*(RE);
    [X,Y] = pol2cart(THETA,RHO);
    fill(X + rE0_JCI(1),Y+ rE0_JCI(2),'c');

    %%% Plotting Europa Surface (3D)
    bodySurface3(RE, [0 0 0], [0 1 1]);

    %%% Plotting Initial Europa-Jupiter Vector
    quiver3(0,0,0,-2*RE,0,0,...
        'linewidth',2,'color',[.7 .7 .7])
    
    %%% Plotting Final Europa-Jupiter Vector
    E2Jhat_f = -rE_JCI(end,:)./norm(rE_JCI(end,:));
    jpECI = quiver3(0, 0, 0, 2*RE*E2Jhat_f(end,1), 2*RE*E2Jhat_f(end,2), 2*RE*E2Jhat_f(end,3),...
        'linewidth',2,'color',[1 .5 0]);
    legend([jpECI], 'Jupiter-Pointing (Final)')
    
    %%% Plotting Intial ECI Velocity
    sc1 = 1000; % scalar
    quiver3(rH0_ECEF(1),rH0_ECEF(2),rH0_ECEF(3),...
        sc1*(vH0_ECI(1)),sc1*(vH0_ECI(2)),sc1*(vH0_ECI(3)),...
        'linewidth',2,'color',[0 0 0])
    
    %%% To Focus on Europa
    xlim([-RE*scale1 RE*scale1])
    ylim([-RE*scale1 RE*scale1])
    zlim([-RE*scale1 RE*scale1])
%     xlim([min(rECI_Hopper(:,1))-.1 max(rECI_Hopper(:,1))+.1])
%     ylim([min(rECI_Hopper(:,2))-.1 max(rECI_Hopper(:,2))+.1])
%     zlim([-.1 .1])
    
    %%% Frame
    title('ECI')
    
end

% ------------------------------------------------------------------------
%%% Plotting the Jupiter System (JCI)
% ------------------------------------------------------------------------
lineWidth = 3;
trackWidth = 2;
if JCIplot == 1
    figure
    hold all
    grid on
    axis equal
    PlotBoi3('X, km','Y, km','Z, km',16)

    %%% Plotting Hopper motion
    % p3 = plot3(r_Hopper(:,1),r_Hopper(:,2),r_Hopper(:,3),'m','linewidth',trackWidth,'markersize',5);
    plot3(StatesI(:,1),StatesI(:,2),StatesI(:,3),'m','linewidth',trackWidth,'markersize',5);

    %%% Plotting Europa motion
    plot3(rE_JCI(:,1),rE_JCI(:,2),rE_JCI(:,3),'--r','linewidth',trackWidth);

    %%% Plotting Europa equator
    th = 0:.01:2*pi;
    x = RE * cos(th) + rEf_JCI(1);
    y = RE * sin(th) + rEf_JCI(2);
    plot(x, y,'b','linewidth',lineWidth);

    %%% Coloring in Europa (2D)
    n = 5000;
    THETA=linspace(0,2*pi,n);
    RHO=ones(1,n)*(RE);
    [X,Y] = pol2cart(THETA,RHO);
    fill(X + rEf_JCI(1),Y+ rEf_JCI(2),'c');

    %%% Plotting Europa-Jupiter Vector
    EJ_hat = -[rEf_JCI(1),rEf_JCI(2),rEf_JCI(3)]/...
        norm(rEf_JCI);
    scale = 2*RE;
    quiver3(rEf_JCI(1),rEf_JCI(2),rEf_JCI(3),...
        scale*EJ_hat(1),scale*EJ_hat(2),scale*EJ_hat(3),...
        'linewidth',2,'color',[1 .5 0])

    %%% Plotting Jupiter Equator
    th = 0:.001:2*pi;
    x = RJ * cos(th);
    y = RJ * sin(th);
    plot(x, y,'k','linewidth',lineWidth);

    %%% Coloring in Jupiter (2D)
    n = 5000;
    THETA=linspace(0,2*pi,n);
    RHO=ones(1,n)*(RJ);
    [X,Y] = pol2cart(THETA,RHO);
    fill(X,Y,[1 .5 0]);
    
    %%% Frame
    title('JCI')
end

% ------------------------------------------------------------------------
%%% ECI Distance vs Time plot
% ------------------------------------------------------------------------
if rECIplot == 1
figure
hold all
plot(Times,ones(size(Times)).*RE,'--b','linewidth',1.5)
plot(Times,arrayfun(@(x) norm(rH_ECI(x,:)), 1:length(Times))','m','linewidth',trackWidth)
PlotBoi2('Time, sec','Distance to Europa Center, km',16)
legend('Europa Mean Radius')
end

% ------------------------------------------------------------------------
%%% East-ness vs Time
% ------------------------------------------------------------------------
if EquatorialEastPlot == 1
%%% Plotting Numerical Eastern Position, Velocity, and Acceleration
figure
subplot(3,1,1)
title('Numerically Calculated Eastern Pos/Vel/Acc')
hold all
plot(Times, EastNumericalPos,'.')
plot([0 Times(end)],[0 0],'--r')
PlotBoi2('', 'East Pos, km', 14)
subplot(3,1,2)
hold all
plot(Times(1:end-2), EastNumericalVel(1:end-1),'.')
plot([0 Times(end)],[0 0],'--r')
PlotBoi2('', 'East Vel, km/s', 14)
subplot(3,1,3)
hold all
plot(Times(1:end-3), EastNumericalAcc(1:end-1),'.')
plot([0 Times(end)],[0 0],'--r')
PlotBoi2('Time, sec', 'East Acc, km/s^2', 14)

%%% Plotting Analytical Acceleration
figure
subplot(3,1,1);hold on
title('Analytically Calculated Eastern Pos/Vel/Acc')
plot(Times, EastAnalyticalPos,'.')
plot(Times(1),EastAnalyticalPos(1),'o','markersize',10)
plot([0 Times(end)],[0 0],'--r')
PlotBoi2('', 'East Pos, km', 14)
subplot(3,1,2);hold on
plot(Times, EastAnalyticalVel,'.')
plot(Times(1),EastAnalyticalVel(1),'o','markersize',10)
plot([0 Times(end)],[0 0],'--r')
PlotBoi2('', 'East Vel, km/s', 14)
subplot(3,1,3);hold on
plot(Times, EastAnalyticalAcc,'.')
plot(Times(1),EastAnalyticalAcc(1),'o','markersize',10)
plot([0 Times(end)],[0 0],'--r')
PlotBoi2('Time, sec','Eastern Acc, km/s^2',14)

% %%% Plotting Tidal Acceleration
% figure
% plot(Times, EastAT)
% title('Eastern Tidal Acceleration')
% PlotBoi2('Times, sec','Eastern Tidal Acceleration, km/s^2',16)

end

% ------------------------------------------------------------------------
%%% Movies
% ------------------------------------------------------------------------
%%% ECEF Movie
if runECEFMovie == 1
    figure
    for i = 1:length(Times) %size(States,1)
        if rem(i,framespeed) == 0 || i == length(Times)
            clf
            hold all
            %%% Plotting Europa motion
            trackWidth = 2;
            plot3(rH_ECEF(i,1),rH_ECEF(i,2),rH_ECEF(i,3),'rX','linewidth',trackWidth,'markersize',10);
            PlotBoi3('X, km','Y, km','Z, km',16)
            grid on
            axis square
            xlim([-RE*scale2 RE*scale2])
            ylim([-RE*scale2 RE*scale2])
            zlim([-RE*scale2 RE*scale2])
%             xlim([min(rECEF_Hopper(:,1)) max(rECEF_Hopper(:,1))])
%             ylim([min(rECEF_Hopper(:,2)) max(rECEF_Hopper(:,2))])
%             zlim([-.01 .01])

            %%% Plotting Past
            plot3(rH_ECEF(1:i-1,1),rH_ECEF(1:i-1,2),rH_ECEF(1:i-1,3),'-m')

            %%% Europa Equator 
            th = 0:.01:2*pi;
            x = RE * cos(th);
            y = RE * sin(th);
            p1 = plot(x, y,'b','linewidth',lineWidth);

            %%% Plotting Europa Surface (3D)
            bodySurface3(RE, [0 0 0], [0 1 1]);


            %%% Plotting Initial Europa-Jupiter Vector
            quiver3(0,0,0,-2*RE,0,0,...
                'linewidth',2,'color',[.7 .7 .7])
            
            %%% View(azimuth,elevation)
%             view(-80,60)
            
            %%% Draw frame
            drawnow limitrate
        end
    end
end

%%% ECI Movie
if runECIMovie == 1
    figure
    for i = 1:length(Times) %size(States,1)
        if rem(i,framespeed) == 0 || i == length(Times)
            clf
            hold all
            %%% Plotting Europa motion
            trackWidth = 2;
            plot3(rH_ECI(i,1),rH_ECI(i,2),rH_ECI(i,3),'rX','linewidth',trackWidth,'markersize',10);
            PlotBoi3('X, km','Y, km','Z, km',16)
            grid on
            axis square
            xlim([-RE*scale2 RE*scale2])
            ylim([-RE*scale2 RE*scale2])
            zlim([-RE*scale2 RE*scale2])

            %%% Plotting Past
            plot3(rH_ECI(1:i-1,1),rH_ECI(1:i-1,2),rH_ECI(1:i-1,3),'-m')

            %%% Europa Equator 
            th = 0:.01:2*pi;
            x = RE * cos(th);
            y = RE * sin(th);
            p1 = plot(x, y,'b','linewidth',lineWidth);

            %%% Plotting Europa Surface (3D)
            bodySurface3(RE, [0 0 0], [0 1 1]);


            %%% Plotting Initial Europa-Jupiter Vector
            quiver3(0,0,0,-2*RE,0,0,...
                'linewidth',2,'color',[.7 .7 .7])

            %%% Plotting Final Europa-Jupiter Vector
            E2Jhat_f = -rE_JCI(i,:)./norm(rE_JCI(i,:));
            quiver3(0, 0, 0, 2*RE*E2Jhat_f(1), 2*RE*E2Jhat_f(2), 2*RE*E2Jhat_f(3),...
                'linewidth',2,'color',[1 .5 0]);
            
            %%% View(azimuth,elevation)
            view(5,30)
            
            %%% Draw frame
            drawnow limitrate
        end
    end
end

% ------------------------------------------------------------------------
%%% Comparing Analytical and Numerical Results (kind of)
% ------------------------------------------------------------------------
if compareAnalyticalNumericalPlots == 1
figure
subplot(3,1,1)
hold all
plot3(rH_ECEF(:,1),rH_ECEF(:,2),rH_ECEF(:,3))
plot3(rH_ECEF(1,1),rH_ECEF(1,2),rH_ECEF(1,3),'o','markersize',10)
plot3(rH_ECEF(end,1),rH_ECEF(end,2),rH_ECEF(end,3),'+','markersize',10)
PlotBoi3('X','Y','Z',12)
title('Numerical Body Position')
view(0,90); axis equal;

subplot(3,1,2)
plot3(relVel(1:end-1,1),relVel(1:end-1,2),relVel(1:end-1,3),'.')
title('Numerical Body Velocity')
PlotBoi3('X','Y','Z',12)
view(0,90); axis equal;

subplot(3,1,3)
plot3(relAcc(1:end-1,1),relAcc(1:end-1,2),relAcc(1:end-1,3),'.')
title('Numerical Body Acceleration')
PlotBoi3('X','Y','Z',12)
view(0,90); axis equal; 

figure
subplot(3,1,1)
hold all
plot3(rAnalytical(:,1),rAnalytical(:,2),rAnalytical(:,3))
plot3(rAnalytical(1,1),rAnalytical(1,2),rAnalytical(1,3),'o','markersize',10)
plot3(rAnalytical(end,1),rAnalytical(end,2),rAnalytical(end,3),'+','markersize',10)

title('Analytical Body Position')
PlotBoi3('X','Y','Z',12)
view(0,90); axis equal; 

subplot(3,1,2)
plot3(vAnalytical(:,1),vAnalytical(:,2),vAnalytical(:,3))
title('Analytical Body Velocity')
PlotBoi3('X','Y','Z',12)
view(0,90); axis equal; 

subplot(3,1,3)
plot3(aH_ECEF_I(:,1),aH_ECEF_I(:,2),aH_ECEF_I(:,3))
title('Analytical Body Acceleration')
PlotBoi3('X','Y','Z',12)
view(0,90); axis equal; 
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('===== NUMERICAL RESULTS =====\n')
% ------------------------------------------------------------------------
%%% Printing Impact Time
% ------------------------------------------------------------------------
if Times(end) ~= tf
    fprintf('Time of Flight:    %1.4f sec\n', Times(end))
elseif Times(end) == tf
    warning off backtrace
    warning('No Impact in Simulation')
else
    fprintf('***Impact Time Error***\n\n')
end

% ------------------------------------------------------------------------
%%% Printing Energy Comparison
% ------------------------------------------------------------------------
fprintf('Initial Velocity:  %4.3f m/s\n', v0_m)
fprintf('Final Velocity:    %4.3f m/s\n', vf_m)
fprintf('Delta-V:           %4.3f m/s\n', dVelocity)
fprintf('Delta-E:           %4.3f J/kg\n', dKE)
fprintf('Angle Change:      %3.4f°\n', AngChange)
fprintf('Distance Traveled: %6.3f m\n', Traveled)
fprintf('Azimuth:           %3.2f°\n',az)

% ------------------------------------------------------------------------
%%% Stating Assumptions
% ------------------------------------------------------------------------
fprintf('\n\n\n***************************\n')
fprintf('ASSUMPTIONS:\n')
fprintf('-Europa and Jupiter are point masses at centers\n')
fprintf('-Europa is spherical\n')
fprintf('-Europa in circular orbit\n')
fprintf('-Jacobian Constant calculation assumes z = 0 (Equatorial)\n')
fprintf('-Plenty of others...\n')
