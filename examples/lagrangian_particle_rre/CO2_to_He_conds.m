%gas parameters
    %names
    name_I='CO2';
    name_II='He';
    %ratios of specific heat
    gamma_I=1.288; %slow material, first phase
    gamma_II=1.667; %fast material, second phase
    %molecular masses:
    mu_I=44.01;
    mu_II=4;
    %Universal gas constant (kg.mol^-1.K^-1)
    R_u=8.314;
    %Phase I gas constant
    R_I=R_u/(mu_I*10^-3);
 
%Experimental conditions:
    %Angle of incidence:
    omega_deg=20;
    omega_rad=omega_deg*pi/180;
    %Incident wave pressure jump
    xi=30;
    %inital temperature:
    T0_I=300; %Phase I
    T0_II=T0_I; %Phase II
    %initial pressure:
    P0=1e5;
    %initial specific volume:
    Vm0_I=R_I*T0_I/P0;
    %distance of shock to interface at t=0 (m):
    experiment_dim=1e-2;
