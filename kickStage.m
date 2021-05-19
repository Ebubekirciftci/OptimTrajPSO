function dy3dt= kickStage(t,y3) 

r=y3(1);                                       % Distance from Earth Center
v=y3(2);                                       % velocity
gamma=y3(3);                                   % Flight Path angle 
beta=y3(4);

%%%%%%%%%%%%%%%%%%%%%%%% Lagrange Multipiers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lamda_r=y3(5);
lamda_v=y3(6);
lamda_gamma=y3(7); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=y3(8);

DU2m=6378.136e3;                % canoncial distance unit to meter
TU2s=806.8;                     % canoncial time unit to time
deg2rad=pi/180;                 % degree to radian


T=120/DU2m*TU2s^2;              % kick Stage thrust (N)
Isp=1961/DU2m*TU2s;             % TSFC (N*s/kg)
R_E=1;                          % Earth Radius (m)
nu_E=1;                         % Earth gravitational parameter (DU^3/TU^2)
omega_E=7.292115e-5*TU2s;       % Earth angular velocity (rad/s)
h=r-R_E;                       % Altitude (DU)  

if h<0
    dmdt=-T/Isp;    
    drdt=0;
    dvdt=0;
    dgammadt=0;
    dbetadt=0;
    dlamda_rdt=0;
    dlamda_vdt=0;
    dlamda_gammadt=0;
else

    alpha=atan(lamda_gamma/(lamda_v*v)); % control angle equation regulated...
    %as d(Hamiltonian)/d(alpha)=0, for aerodynamic parameters equal to zero

    if alpha>10*deg2rad
        alpha=10*deg2rad;
    elseif alpha<-10*deg2rad
        alpha=-10*deg2rad;
    end


    % Aerodynmamic forces and derivatives equal to zero, because kick stage ...
    %ignition start at exoatmospheric altitude

    L=0;       % Lift force
    D=0;       % Drag force
    dLdr=0;    % deriv. of Lift with respect to distance from Earth center 
    dDdr=0;    % deriv. of Drag with respect to distance from Earth center
    dLdv=0;    % deriv. of Lift with respect to velocity
    dDdv=0;    % deriv. of Drag with respect to velocity



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% DIFERANTIAL EQUATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    dmdt=-T/Isp;    

    drdt=v*sin(gamma);
    dvdt=T/m*cos(alpha)-nu_E/r^2*sin(gamma)-...
        D/m+omega_E^2*r*sin(gamma);
    dgammadt=T/(m*v)*sin(alpha)+(v/r-nu_E/(r^2*v))*...
        cos(gamma)+L/(m*v)+2*omega_E+omega_E^2*r/v*...
        cos(gamma);
    dbetadt=v*cos(gamma)/r;


    dlamda_rdt=-(-lamda_gamma*v*cos(gamma)/r^2+...
        2*lamda_gamma*nu_E*cos(gamma)/(r^3*v)+...
        lamda_gamma*omega_E^2*cos(gamma)/v+...
        2*lamda_v*nu_E*sin(gamma)/r^3+...
        lamda_v*omega_E^2*sin(gamma)+...
        lamda_gamma*dLdr/(v*m)-lamda_v*dDdr/m);

    dlamda_vdt=-(lamda_r*sin(gamma)-...
        lamda_gamma*T*sin(alpha)/(m*v^2)+...
        lamda_gamma*cos(gamma)/r+lamda_gamma*...
        nu_E*cos(gamma)/(r^2*v^2)-lamda_gamma*...
        omega_E^2*r*cos(gamma)/v^2+lamda_gamma*...
        (dLdv/(v*m)-L/(m*v^2))-lamda_v*dDdv/m);

    dlamda_gammadt=-(lamda_r*v*cos(gamma)-...
        lamda_gamma*v*sin(gamma)/r+lamda_gamma*...
        nu_E*sin(gamma)/(r^2*v)-lamda_gamma*...
        omega_E^2*r*sin(gamma)/v-lamda_v*nu_E*cos(gamma)/r^2+...
        lamda_v*omega_E^2*r*cos(gamma));
end

dy3dt=[drdt;dvdt;dgammadt;dbetadt;dlamda_rdt;...
    dlamda_vdt;dlamda_gammadt;dmdt];
end