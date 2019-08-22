%gas parameters
    %names
    name_I='2/3*H2+1/3*O2';
    name_II='He';
    %ratios of specific heat
    gamma_I=1.4016; %slow material, first phase
    gamma_II=1.667; %fast material, second phase
    %molecular masses:
    mu_I=2/3*2+1/3*16;
    mu_II=4;
    %Universal gas constant (kg.mol^-1.K^-1)
    R_u=8.314;
    %Phase I gas constant
    R_I=R_u/(mu_I*10^-3);
 
%Experimental conditions:
    %Angle of incidence:
        %Cond Remy
        %omega_deg=25;
        
        %Pour cond RCM
        omega_deg=20;
        
    omega_rad=omega_deg*pi/180;
    %Incident wave pressure jump
        %Cond Remy
            %xi=18.5;
        
        %Cond RCM
            xi=10.4;
    %inital temperature:
        %cond Remy
        T0_I=298; %Phase I
        T0_II=1138; %Phase II
        
        %cond RCM
        %T0_I=1000;
        %T0_II=2500;
    %initial pressure:
        %Cond Remy
        %P0=1e5
        
        %Cond RCM
        P0=1e5;
    %initial specific volume:
    Vm0_I=R_I*T0_I/P0;
    %distance of shock to interface at t=0 (m):
    experiment_dim=1e-2;