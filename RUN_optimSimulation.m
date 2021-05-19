clc,clear
deg2rad=pi/180;                 % degree to radian

%%%%%%%%%%%%%%%%%%%%%%%% Optimization Innputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lb_lamdas=-1;               % Lower boundary of Adjoint Variables
ub_lamdas=1;                % Upper boundary of Adjoint Variables
d_lamdas=3;                 % Number of Adjoint Variables
lb_nonIt3=0.5;              % Lower boundary of non-Thrust Dynamic Process
ub_nonIt3=3;                % Upper boundary of non-Thrust Dynamic Process
d_nonIt3=1;                 % Number of non-Thrust Dynamic Process
lb_It4=0.01;                % Lower boundary of Kick Stage Dynamic Process
ub_It4=0.405;               % Upper boundary of Kick Stage Dynamic Process
d_It4=1;                    % Number of Kick Stage Dynamic Process
lb_gamma=82*deg2rad;        % Lower boundary of first Flight Path Angle
ub_gamma=90*deg2rad;        % Upper boundary of first Flight Path Angle    
d_gamma=1;                  % Number of first Flight Path Angle                 

ssize=50;                   % Swarm Size
w=1;                        % Inertial coefficient
c1=1;                       % Cognitive coefficient
c2=1;                       % social coefficient
iter_number=1000;           % İteration Number

[sbestpos,sbestval,objit] = PSO(lb_lamdas,ub_lamdas,d_lamdas,lb_nonIt3,...
    ub_nonIt3,d_nonIt3,lb_It4,ub_It4,d_It4,ub_gamma,lb_gamma,d_gamma,...
    ssize,w,c1,c2,iter_number);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% BEST SOLUTION RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DU2m=6378.136e3;                 % canoncial distance unit to meter
DU2km=6378.136;                  % canoncial distance unit to kilometer
TU2s=806.8;                      % canoncial time unit to time
nu_E=1;                        % Earth gravitational parameter (DU^3/TU^2)
R_E=1;                         % Earth Radius (DU)
omega_E=7.292115e-5*TU2s;      % Earth angular velocity (rad/TU)


m1=11300;                                % First stage mass
mp1=10340;                               % first stage propellant mass
ms1=m1-mp1;                              % First stage structural mass
m2=2900;                                 % second stage mass
mp2=2640;                                % second stage propellant mass
ms2=m2-mp2;                              % second stage structural mass
ms3=130;                                 % kick stage structural mass
payload=150;                             % Default Payload mass
T1=192100/DU2m*TU2s^2;                   % first stage thrust (N)
Isp1=2972/DU2m*TU2s;                     % Spesific Thrust (N*s/kg)
TSFC1=T1/Isp1;                           % Fuel consumption (kg/s)
It1=mp1/TSFC1;                           % Ignition time  (s) 

T2=22200/DU2m*TU2s^2;                    % first stage thrust 
Isp2=3266/DU2m*TU2s;                     % Spesific Thrust (N*s/kg)
TSFC2=T2/Isp2;                           % Fuel consumption (kg/s)
It2=mp2/TSFC2;                           % Ignition time  (s) 
nonIt3=sbestpos(4);
It4=sbestpos(5);
options = odeset('RelTol',1e-2,'AbsTol',1e-2);


r1=R_E;                      % first stage first distance from Earth center
v1=10/DU2m*TU2s;             % first stage first velocity after Lift-off
gamma1=sbestpos(6);          % first stage first path angle
beta1=0;                     % first stage first latitude value
 
%%%%%%%%%%%%%%%%%%% first value of lagrange multipiers %%%%%%%%%%%%%%%%%%%
lamda_r1=sbestpos(1);                           
lamda_v1=sbestpos(2);
lamda_gamma1=sbestpos(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mass1=14500;                   % Lift-off mass

[t1,y1] = ode45(@stage1,[0 It1],[r1;v1;gamma1;beta1;lamda_r1;lamda_v1;...
    lamda_gamma1;mass1],options);     % differential process first stage

r2=y1(end,1);              % second stage first distance from Earth center 
v2=y1(end,2);              % second stage first velocity value
gamma2=y1(end,3);          % second stage first path angle
beta2=y1(end,4);           % second stage first latitude value


lamda_r2=y1(end,5);       
lamda_v2=y1(end,6);
lamda_gamma2=y1(end,7);

mass2=y1(end,8)-ms1;

[t2,y2] = ode45(@stage2,[It1 It1+It2],...
    [r2;v2;gamma2;beta2;lamda_r2;lamda_v2;...
    lamda_gamma2;mass2],options);     % differential process second stage


r3=y2(end,1);          % nonThrust Process first distance from Earth center 
v3=y2(end,2);          % nonThrust Process first velocity value
gamma3=y2(end,3);      % nonThrust Process first path angle
beta3=y2(end,4);       % nonThrust Process first latitude value


lamda_r3=y2(end,5);       
lamda_v3=y2(end,6);
lamda_gamma3=y2(end,7);

mass3=y2(end,8)-ms2;


[t3,y3] = ode45(@nonThrust,[It1+It2 It1+It2+nonIt3],...
    [r3;v3;gamma3;beta3;lamda_r3;lamda_v3;...
    lamda_gamma3;mass3],options);     % differential process second stage



r4=y3(end,1);          % kick Stage first distance from Earth center 
v4=y3(end,2);          % kick Stage first velocity value
gamma4=y3(end,3);      % kick Stage first path angle
beta4=y3(end,4);       % kick Stage first latitude value


lamda_r4=y3(end,5);       
lamda_v4=y3(end,6);
lamda_gamma4=y3(end,7);

mass4=y3(end,8);


[t4,y4] = ode45(@kickStage,[It1+It2+nonIt3 It1+It2+nonIt3+It4],...
    [r4;v4;gamma4;beta4;lamda_r4;lamda_v4;...
    lamda_gamma4;mass4],options);     % differential process second stage

r_final=y4(end,1);          
v_final=y4(end,2);         
gamma_final=y4(end,3);      


mass_final=y4(end,8)-ms3;


figure
plot(t1*TU2s,DU2km.*[y1(:,1)-R_E*ones(height(t1),1)],...
   t2*TU2s,DU2km.*[y2(:,1)-R_E*ones(height(t2),1)],...
   t3*TU2s,DU2km.*[y3(:,1)-R_E*ones(height(t3),1)],...
   t4*TU2s,DU2km.*[y4(:,1)-R_E*ones(height(t4),1)]);
ylabel('irtifa(km)');
xlabel('zaman(s)');


figure
plot(y1(:,4)./deg2rad,DU2km.*[y1(:,1)-R_E*ones(height(t1),1)],...
   y2(:,4)./deg2rad,DU2km.*[y2(:,1)-R_E*ones(height(t2),1)],...
   y3(:,4)./deg2rad,DU2km.*[y3(:,1)-R_E*ones(height(t3),1)],...
   y4(:,4)./deg2rad,DU2km.*[y4(:,1)-R_E*ones(height(t4),1)]);

ylabel('Altitude(km)');
xlabel('Longitude(degree)');

figure

semilogy(objit)
ylabel('Objective Function Value');
xlabel('Iteration Number');
grid on