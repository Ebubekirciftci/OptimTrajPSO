function [L,D,dLdr,dDdr,dLdv,dDdv,dLdalpha,dDdalpha] = AeroCalculation(alpha,v,h)


%%%%%%%%%%%%%%%%%%%%%5%%% Atmosphere Modeling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho=1.225*exp(-h/9300);                    % Air density (kg/m^3)
drhodr=-1.225*exp(-h/9300);                
diameter=1.2;
A=1/4*diameter^2*pi;

if h>43000                                 % Air Pressure (Pa)
    Pa=0;
    dPadr=0;
else
    Pa=101325*(1-2.25577e-5*h)^5.25588;
    dPadr=-2.25577e-5*5.25588*101325*(1-2.25577e-5*h)^4.25588;
end

if h>10000
    a=294.9;
    dadr=0;
else
    a=sqrt(1.4*Pa/rho);
    dadr=sqrt(1.4)*0.5*(Pa/rho)^(-0.5)*...
        (dPadr*rho-drhodr*Pa)/rho^2;
        
end

M=v/a;                                        % Mach Number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% Aerodinamic Forces and Derivatives Calculations %%%%%%%%%%%%%%

if M>=10 || h>100000
    C_L=0;
    C_D=0;
    dC_LdM=0;
    dC_DdM=0;
    dC_Ldalpha=0;
else 
    if alpha>=0
        dC_Ldalpha=-7.43866786811709e-05*M^4+0.00189384537172436*M^3-0.0172488689806825*M^2+0.0628647820832881*M+0.0434061151786836;

        
        dC_LdalphadM=-4*7.43866786811709e-05*M^3+3*0.00189384537172436*M^2-2*0.0172488689806825*M+0.0628647820832881;
        
        C_L=dC_Ldalpha*alpha+(3.35087702133498e-05*M^5-0.000875486311164452*M^4+0.00755257870750183*M^3-0.0198014630489914*M^2-0.0236641997330676*M-0.0245102713238283);
        C_D=7.10569954010629e-05*M^5-0.00267889219860829*M^4+0.0373289362068858*M^3-0.231690183062741*M^2+0.586154598949159*M+0.0709482667411191;
        
        
        
        dC_LdM=dC_LdalphadM*alpha+(5*3.35087702133498e-05*M^4-4*0.000875486311164452*M^3+3*0.00755257870750183*M^2-2*0.0198014630489914*M-0.0236641997330676);
        dC_DdM=5*7.10569954010629e-05*M^4-4*0.00267889219860829*M^3+3*0.0373289362068858*M^2-2*0.231690183062741*M+0.586154598949159;
        
    else
        dC_Ldalpha=-(-7.43866786811709e-05*M^4+0.00189384537172436*M^3-0.0172488689806825*M^2+0.0628647820832881*M+0.0434061151786836);

        dC_LdalphadM=-(-4*7.43866786811709e-05*M^3+3*0.00189384537172436*M^2-2*0.0172488689806825*M+0.0628647820832881);
        
        
        
        C_L=-(dC_Ldalpha*alpha+(3.35087702133498e-05*M^5-0.000875486311164452*M^4+0.00755257870750183*M^3-0.0198014630489914*M^2-0.0236641997330676*M-0.0245102713238283));
        
        C_D=7.10569954010629e-05*M^5-0.00267889219860829*M^4+0.0373289362068858*M^3-0.231690183062741*M^2+0.586154598949159*M+0.0709482667411191;
        
        dC_LdM=-(dC_LdalphadM*alpha+(5*3.35087702133498e-05*M^4-4*0.000875486311164452*M^3+3*0.00755257870750183*M^2-2*0.0198014630489914*M-0.0236641997330676));
        dC_DdM=5*7.10569954010629e-05*M^4-4*0.00267889219860829*M^3+3*0.0373289362068858*M^2-2*0.231690183062741*M+0.586154598949159;
    end
    
end


dC_Ldr=-dC_LdM*(v/a^2)*dadr;
dC_Ddr=-dC_DdM*(v/a^2)*dadr;

L=1/2*rho*v^2*C_L*A;
D=1/2*rho*v^2*C_D*A;
dLdr=1/2*C_L*A*v^2*drhodr-1/2*rho*A*v^2*dC_Ldr;
dDdr=1/2*C_D*A*v^2*drhodr-1/2*rho*A*v^2*dC_Ddr;
dLdv=v*rho*A*C_L+1/2*rho*v^2*A/a*dC_LdM;
dDdv=v*rho*A*C_D+1/2*rho*v^2*A/a*dC_DdM;
dLdalpha=1/2*rho*v^2*dC_Ldalpha*A;
dDdalpha=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
