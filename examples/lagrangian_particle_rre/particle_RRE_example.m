H2_O2_to_He_conds;

%experiment dimentions
    exp_dims=1e-2;
%Initial particle position when encountering shock wave (m)
    x=.5e-2; %x coordinate
    y=3e-2; %y coordinate
    %distance to point of incidence, when shock reaches particle
        r_0=y-x/tan(omega_rad); %has to be positive!
        if r_0<0
            "This code only follows particles in phase I!"
        else
        
%When in zone 1 between incident shock and expansion fan
    Msh=sqrt(xiToSqMach(xi,gamma_I,pi/2));
    Mi=Msh/sin(omega_rad);
    delta_1=atan(sqrt(tanDefSq(xi,Mi,gamma_I)));
    mu_1=asin(1/Mi); %velocity vector-expansion fan angle at beginning...
        %... of expansion fan
    t_1=1/sqrt(gamma_I*R_I*T0_I)*r_0*sin(omega_rad-delta_1+mu_1); %time spent...
        %between the incident shock and the expansion
    temp_jump=xiToTempJ(xi,gamma_I); %temperature jump of incident shock
    spec_vol_jump=((gamma_I-1)*xi+(gamma_I+1))/((gamma_I+1)*xi+(gamma_I-1));
    r_1=r_0/(Mi*sin(omega_rad-delta_1)); %radius when arriving at expansion


%When in zone 2, behind the expansion fan
    Mr=sqrt(postShockMachSq(xi,Mi,gamma_I));%Mach when arriving at reflection

    %Computing expansion deflection and pressure jump
    Mt=sqrt((gamma_I*T0_I*mu_II)/(gamma_II*T0_II*mu_I))*Mi;
    npts=200; %number of points
    xis=logspace(0,log(xi)/log(10),npts); %pressure jumps when going through...
        %... expansion fan
    deltas_r=ones(1,npts)*delta_1; %initializing array of expansion deflections
    deltas_t=zeros(1,npts); %initializing array of transmitted wave deflections
    for i=1:npts
        M_post_exp=sqrt(2/(gamma_I-1)*((1+(gamma_I-1)/2*Mr^2)*(xis(i)/xi)^...
            ((1-gamma_I)/gamma_I)-1)); %post-expansion Mach at given pressure jump
        deltas_r(i)= deltas_r(i) + prandtlMeyer(M_post_exp,gamma_I)-...
            prandtlMeyer(Mr,gamma_I); %computing deltas_r

        deltas_t(i)=atan(sqrt(tanDefSq(xis(i),Mt,gamma_II))); %computing...
            %... transmitted wave deflection
    end

    def_diff=abs(deltas_r-deltas_t);
    [m,k]=min(def_diff);
    delta_2=deltas_r(k)-delta_1; %deviation at end of expansion
    xi_t=xis(k); %pressure jump after expansion
    M_post_exp=sqrt(2/(gamma_I-1)*((1+(gamma_I-1)/2*Mr^2)*(xi_t/xi)^...
            ((1-gamma_I)/gamma_I)-1)); %Mach after expansion
    exp_temp_jump=(1+(gamma_I-1)/2*Mr^2)...
        /(1+(gamma_I-1)/2*M_post_exp^2); %temperature jump after expansion
    T_post_exp=T0_I*temp_jump*exp_temp_jump; %temperature after expansion
    exp_spec_vol_jump=((1+(gamma_I-1)/2*Mr^2)...
        /(1+(gamma_I-1)/2*M_post_exp^2))^(1/(1-gamma_I)); %spec. vol. jump...
            %... after expansion
    Vm_post_exp=Vm0_I*spec_vol_jump*exp_spec_vol_jump;
    
%In expansion fan
    no_dev_pts=1000; %number of points to solve differential equation
    ddelta=delta_2/(no_dev_pts); %deviation step

    T1=T0_I*temp_jump; %temperature before expansion fan

    u1=Mr*sqrt(gamma_I*R_I*T1); %speed behind expansion fan
    u1nexp=u1/Mr; %speed normal to the beginning of the expansion fan
    u1texp=sqrt(u1^2-u1nexp^2); %speed tangent to the beginning of expansion fan
    tsexp=zeros(1,no_dev_pts); %initializing time array
    r=r_0*Mr*sin(omega_rad-delta_1); %initializing radius i.e. distance to...
        %... point of incidence

    Msqexp=Mr^2; %initializing square of Mach in expansion fan
    Texps=zeros(1,no_dev_pts); %initializing array of temperatures in expansion fan
    Pexps=zeros(1,no_dev_pts); %initializing array of pressures in expansion fan
    Vmexps=zeros(1,no_dev_pts); %initializing array of spec. vols. in expansion fan
    Pexp=P0*xi; %pressure at current position
    Texp=T1; %temperature at current position
    Vmexp=Vm0_I*spec_vol_jump; %spec. vol. at current position

    for i=1:no_dev_pts
        dMsqexp=2*Msqexp*(1+(gamma_I-1)/2*Msqexp)/sqrt(Msqexp-1)*ddelta; %square ...
            %... Mach variation

        %Computing pressure and temperature
        dPexp=-Pexp*(gamma_I/2)/(1+(gamma_I-1)/2*Msqexp)*dMsqexp; %pressure variation
        Pexp=Pexp+dPexp; %new pressure
        dTexp=-Texp*((gamma_I-1)/2)/(1+(gamma_I-1)/2*Msqexp)*dMsqexp; %temperature variation
        Texp=Texp+dTexp; %new temperature
        dVmexp=Vmexp/2*1/(1+(gamma_I-1)/2*Msqexp)*dMsqexp; %spec. vol. variaton
        Vmexp=Vmexp+dVmexp;
        Pexps(i)=Pexp; %adding pressure to array
        Texps(i)=Texp; %adding temperature to array
        Vmexps(i)=Vmexp; %adding spec. vol. to array

        %geometry calculations
        dtheta=ddelta+1/(2*sqrt(Msqexp^2-Msqexp))*dMsqexp; %angle variation...
            %... in cylindrical coordinates
        dt=r/u1nexp*dtheta; %time interval
        dr=u1texp*dt; %radius variation
        du1nexp=u1nexp/2*1/(1+(gamma_I-1)/2*Msqexp)*dMsqexp; %variation of...
            %... speed normal to small disturbance 

        r=r+dr; %new radius
        tsexp(i)=dt; %adding time to array
        u1nexp=u1nexp+du1nexp; %new normal speed

        %changing Msqexp
        Msqexp=Msqexp+dMsqexp;
    end
    tsexp; %time spent in expansio fan

    %summing time intervals
    for i=2:no_dev_pts
        tsexp(i)=tsexp(i-1)+tsexp(i);
    end

% figure
% plot([0,tsexp],[P0*xi,Pexps])
% P0*xi_t


refr_dur=experiment_dim/(Msh*sqrt(gamma_I*R_I*T0_I)); %time needed for...
    %... incident wave to reach interface
ts=zeros(1,no_dev_pts+6);
t=refr_dur;
ts(2)=t;
t=t+x/(Msh*sqrt(gamma_I*R_I*T0_I));
ts(3:4)=[t,t];
t=t+t_1;
ts(5)=t;
ts(6:5+no_dev_pts)=t+tsexp;
t=ts(5+no_dev_pts)+refr_dur;
ts(6+no_dev_pts)=t;
ts_plot=ts*1e6; %in micros-secs

Ps=[P0,P0,P0,P0*xi,P0*xi,Pexps(1:end-1),P0*xi_t,P0*xi_t]; %presure array for plot
%figure
%hold on
subplot(3,1,1)
Ps_plot=Ps/10^5;
plot(ts_plot,Ps_plot) %plotting pressure of particle
title('Pressure evolution')
xlabel("time (micros)")
ylabel("Pressure (bar)")

%figure
subplot(3,1,2)
Ts=[T0_I,T0_I,T0_I,T0_I*temp_jump,T0_I*temp_jump,Texps(1:end-1),...
    T_post_exp,T_post_exp];
plot(ts_plot,Ts) %plotting temperature of particle
title('Temp evolution (K)')
xlabel("time (micros)")
ylabel("Temperature (K)")


%figure
subplot(3,1,3)
Vms=[Vm0_I,Vm0_I,Vm0_I,Vm0_I*spec_vol_jump,Vm0_I*spec_vol_jump,Vmexps(1:end-1),...
    Vm_post_exp,Vm_post_exp];
Vms_norm=Vms/Vm0_I;
plot(ts_plot,Vms_norm) %plotting temperature of particle
title('Spec. vol. evolution')
xlabel("time (micros)")
ylabel("Normalized specific volume")
ylim([0,1.1])

%exporting normalized specific volume
% VTIM_mat=[ts',Vms_norm'];
% 
% dlmwrite('VTIM.dat',VTIM_mat,'delimiter',' ','precision','%1.6e')

% %exporting temperature and pressure
%file_name='P_rre_omg_20_xi_10p4_x_8_y_23.csv';
% fid = fopen( file_name , 'wt' );
% fprintf(fid,'Time(sec), Pressure(Pa)\n');
% pressure_file_mat=[ts;Ps]
% fprintf(file_name,'%1.6e, %1.6e\n',pressure_file_mat);
% fclose(fid);
% pressure_file_mat=[ts',Ps'];
% pressure_file_cell=num2cell(pressure_file_mat);
% T=cell2table(pressure_file_cell);
% {'Time(sec)' 'Pressure(Pa)'}
% T.Properties.VariableNames={'Time' 'Pressure'};
% writetable(T,file_name)
end