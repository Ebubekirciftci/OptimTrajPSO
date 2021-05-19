function dy2dt= stage2(t,y2) 

r=y2(1);                                     % Distance from Earth Center
v=y2(2);                                     % velocity
gamma=y2(3);                                 % Flight Path angle 
beta=y2(4);
%%%%%%%%%%%%%%%%%%%%%%%% Lagrange Multipiers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lamda_r=y2(5);
lamda_v=y2(6);
lamda_gamma=y2(7); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=y2(8);                        % mass    

DU2m=6378.136e3;                % canoncial distance unit to meter
TU2s=806.8;                     % canoncial time unit to time
deg2rad=pi/180;                 % degree to radian


T=22200/DU2m*TU2s^2;            % Stage2 thrust (N)
Isp=3266/DU2m*TU2s;             % TSFC (N*s/kg)
R_E=1;                          % Earth Radius (m)
nu_E=1;                         % Earth gravitational parameter (DU^3/TU^2)
omega_E=7.292115e-5*TU2s;       % Earth angular velocity (rad/s)
h=r-R_E;                        % Altitude (m)  

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
    h_aero=h*DU2m;              % canoncial distance unit to meter for altitude
    v_aero=v*DU2m/TU2s;         % canoncial velocity to m/sec for velocity
    
    
    
    
    alpha_deg=-10:0.5:10;                     % possible angle of attack values
    alpha=deg2rad.*alpha_deg;
    
    
    % For loop writed for calculating Aerodinamic forces and derivatives...
    % for every possible angle of attack value
    
    for i=1:length(alpha)
        
        [L(i),D(i),dLdr(i),dDdr(i),dLdv(i),dDdv(i),dLdalpha(i),dDdalpha(i)]...
            = AeroCalculation(alpha_deg(i),v_aero,h_aero);
        L(i)=L(i)/DU2m*TU2s^2;
        D(i)=D(i)/DU2m*TU2s^2;
        dLdr(i)=dLdr(i)*TU2s^2;
        dDdr(i)=dDdr(i)*TU2s^2;
        dLdv(i)=dLdv(i)*TU2s;
        dDdv(i)=dDdv(i)*TU2s;
        dLdalpha(i)=dLdalpha(i)/DU2m*TU2s^2/deg2rad;
        dDdalpha(i)=dDdalpha(i)/DU2m*TU2s^2/deg2rad;
        
        dHdalpha(i)=abs(lamda_gamma*T/(m*v)*cos(alpha(i))-lamda_gamma*...
            1/(m*v)*dLdalpha(i)-lamda_v*T/m*sin(alpha(i))-...
            lamda_v/m*dDdalpha(i));
        
    end
    
    %Then, alpha value that caused to minimum hamiltonian function value is...
    %selected and Aerodinamic forces and derivatives are selected for this...
    %alpha value.
    
    
    [val, idx] = min(dHdalpha);
    alpha=alpha(idx);
    L=L(idx);       % Lift force
    D=D(idx);       % Drag force
    dLdr=dLdr(idx); % deriv. of Lift with respect to distance from Earth center
    dDdr=dDdr(idx); % deriv. of Drag with respect to distance from Earth center
    dLdv=dLdv(idx); % deriv. of Lift with respect to velocity
    dDdv=dDdv(idx); % deriv. of Drag with respect to velocity
    
    
    
    
    
    
    
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

dy2dt=[drdt;dvdt;dgammadt;dbetadt;dlamda_rdt;...
    dlamda_vdt;dlamda_gammadt;dmdt];
end