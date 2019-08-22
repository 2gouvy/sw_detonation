%%Computing boundaries between different refraction systems
    %%...for Ar->He slow-fast refraction  in...
    %%... interface angle- incident pressure jump plane.
    %%A similar graph can be found in Abd-El-Fattah's and L. F. ...
    %%... Henderson 1978 paper Shock Waves at slow-fast gas interface, fig
    %%.. 13.
    %%The computed boundaries are compared with the article's data and we
    %%... can observe an overall good agreement, except for the FPR-TNR...
    %%... boundary.

limits_computation_time=cputime;%initialising timer for computation time

%gas parameters
    %gas names
    name_I='2/3*H2+1/3*O2';
    name_II='He';
    %ratios of specific heat
    gamma_I=1.4016; %slow material, first phase
    gamma_II=1.667; %fast material, second phase
    %molecular masses:
    mu_I=2/3*2+1/3*16;
    mu_II=4;
    %phase temperatures (K):
    T_0_I=298;
    T_0_II=1138;
    temp_ratio=T_0_I/T_0_II;

%Computing RRE<->... limit with graphical method:
    %Principle: Doing a dichotomy on "isSfStrongRRE" function to find omega...
        %at given chi.
        %"isSfStrongRRE" uses graphical method to determine if...
        %..."RRE->..." transition was made
    nchis_RRE=150;nys_RRE=100;%nchis_RRE: no. points on chi-axis
        %nys_RRE: no. points on xi axis when computing polar
    chis_RRE=linspace(0,1,nchis_RRE);%chi coordinates
    omegas_RRE=zeros(1,nchis_RRE);%pre-allocated omega degree coordinates
    rre_bound_res=.01;%omega coordinate precision, dichotomy termination threshold

    %executing dichotomy for each chi value
    for i=2:nchis_RRE-1 %not solving for chis_RRE(1)=0 and chis_RRE(end)=1 limits
        omega_inf_rre=0;omega_sup_rre=90;
        while omega_sup_rre-omega_inf_rre>rre_bound_res
            next_omega=(omega_sup_rre+omega_inf_rre)/2;
            is_sf_RRE=isSfStrongRRE(1/chis_RRE(i),pi/180*next_omega,gamma_I,...
                gamma_II,mu_I,mu_II,nys_RRE,temp_ratio);
            if is_sf_RRE
                omega_inf_rre=next_omega;
            else
                omega_sup_rre=next_omega;
            end
        end
        omegas_RRE(i)=omega_inf_rre;
    end
 
%Computing RRR<->BPR limit in omega-xi plane:
    nchis=100;nys=400;%nchis: number of points along chi axis;...
        %nys: nb of points along xi axis when computing polars

    %Principle: For each value of chi, the omega coordinate is computing by...
        %...doing a dichotomy with the isBeforeSfRRRToBPR, which determines if the...
        %... RRR->BPR transition has taken place at given chi and omega
    rrr_bpr_boundary_precision=.01; %omega precision, dichtomy termination...
        %...threshold

    chis_RRR_BPR=linspace(.52,1,nchis);%chi coordinates
    omegas_RRR_BPR=zeros(1,nchis);%pre-allocated omega coordinates
    for i=1:nchis
        omega_inf=0;omega_sup=90;%initial dichotomy bounds
        while omega_sup-omega_inf>rrr_bpr_boundary_precision %dichotomy...
            %...termination criteria

            %executing dichotomy, with lower bound having true value and...
                %... upper bound having false value
            next_omega=(omega_sup+omega_inf)/2;
            if isBeforeSfRRRToBPR(1/chis_RRR_BPR(i),pi/180*next_omega,gamma_I,...
                    gamma_II,mu_I,mu_II,nys,temp_ratio)
                omega_inf=next_omega;
            else
                omega_sup=next_omega;
            end
        end
        omegas_RRR_BPR(i)=omega_inf; %adding omega coordinate in degrees to...
            %...preallocated array
    end

%Computing BPR<->FNR boundary
    %Principle: Computes omegas for which Mt=1 at different values of chi
    ny_BPR_FNR=100;%number of points on chi axis
    chis_BPR_FNR=linspace(0,1,ny_BPR_FNR);%chi coordinates
    omegas_BPR_FNR=zeros(1,ny_BPR_FNR);%pre-allocated omega coordinates
    %calculating omega coordinates
    for i=1:ny_BPR_FNR
        omegas_BPR_FNR(i)=bPRFNROmega(1/chis_BPR_FNR(i),gamma_I,gamma_II...
            ,mu_I,mu_II,temp_ratio);
    end

%Computes TNR<->LSR limit
    %Principle: Solves the Mt=1 equation
    npts=400; %number of points on chi axis
    chis_TNR_LSR=linspace(0,1,npts);
    omegas_TNR_LSR=zeros(1,npts);
    syms omega
    for i=2:npts
        %solving Mt=1 equation at specific omega
        xi=1/chis_TNR_LSR(i);
        eqn = postShockMachSq(xi,sqrt(xiToSqMach(xi,gamma_I,omega)),gamma_I)==1;
        a=[vpasolve(eqn,omega,[0,pi/2])];
        if ~(isempty(a))
            omegas_TNR_LSR(i)=180/pi*a;
        end
    end

%Computing FPR<->TNR boundary
%Principle: Graphically resolves the beta1=beta2 problem, as described...
    %... in Abd-El-Fattah, L. F. Henderson: Shock waves at a slow-fast...
    %... gas interface (1978), p88, figure 6b. The boundary is computed...
    %... with a dichotomy on isBeforeFPRToTNR at constant chi
txt="Computing FPR-TNR boundary, this may take a while";
txt
nchis_tnr=200; %number of points on chi axis
nys_tnr=30; %number of points in computed polars
chis_FPR_TNR=linspace(1e-5,1,nchis_tnr); %chi values
omegas_FPR_TNR=zeros(1,nchis_tnr);%pre-allocated omega array
Msj=0; %j-wave Mach initialization
fpr_tnr_boundary_precision=.01; %dichotomy termination threshold
for i=1:nchis_tnr
    solve_Msj_eq=true;
    omega_inf_tnr=0;omega_sup_tnr=90;
    while omega_sup_tnr-omega_inf_tnr>fpr_tnr_boundary_precision
        next_omega=(omega_sup_tnr+omega_inf_tnr)/2;
        [is_before_trans,Msj]=isBeforeFPRToTNR(1/chis_FPR_TNR(i),pi/180*next_omega,...
             gamma_I,gamma_II,mu_I,mu_II,nys_tnr,solve_Msj_eq,Msj,temp_ratio);
        if solve_Msj_eq %Msj is only computed at the beginning...
                %... of the dichotomy
            solve_Msj_eq=false;
        end
        if is_before_trans
            omega_inf_tnr=next_omega;
        else
            omega_sup_tnr=next_omega;
        end
    end
    omegas_FPR_TNR(i)=omega_inf_tnr;
end

limits_computation_time=cputime-limits_computation_time %stopping timer


%Ploting limits like article
figure
hold on
plot(omegas_RRE(2:end-1),chis_RRE(2:end-1)) %ploting RRE<->... limit
hold on
plot(omegas_RRR_BPR,chis_RRR_BPR) %ploting RRR<->BPR limit
hold on
plot(omegas_BPR_FNR,chis_BPR_FNR) %ploting BPR<->FNR lim
hold on
plot(omegas_FPR_TNR(10:end),chis_FPR_TNR(10:end)) %ploting FPR<->TNR lim
hold on
plot(omegas_TNR_LSR(2:end),chis_TNR_LSR(2:end))  %ploting TNR<->LSR lim
legends={"RRE->... (graphical resolution)",...
    "RRR<->BPR","BPR<->FNR","FPR<->TNR",...
    "TNR<->LSR"};

legend(legends,'Location','eastoutside')
xlabel("$\omega_i$ (deg)",'interpreter','latex')
ylabel("\chi")
title(['Computed boundaries for slow-fast ' name_I '->' name_II ' refraction'])
xlim([15,90])
ylim([0,1])
hold off


%Plotting limits with Mach instead of shock strength
incident_Mach_shock_RRE=arrayfun(@xiToMach,1./chis_RRE(2:end-1));
incident_Mach_shock_RRR_BPR=arrayfun(@xiToMach,1./chis_RRR_BPR)
incident_Mach_shock_BPR_FNR=arrayfun(@xiToMach,1./chis_BPR_FNR)
incident_Mach_shock_FPR_TNR=arrayfun(@xiToMach,1./chis_FPR_TNR)
incident_Mach_shock_TNR_LSR=arrayfun(@xiToMach,1./chis_TNR_LSR(2:end))
figure
hold on
plot(omegas_RRE(2:end-1),incident_Mach_shock_RRE) %ploting RRE<->... limit
hold on
plot(omegas_RRR_BPR,incident_Mach_shock_RRR_BPR) %ploting RRR<->BPR limit
hold on
plot(omegas_BPR_FNR,incident_Mach_shock_BPR_FNR) %ploting BPR<->FNR lim
hold on
plot(omegas_FPR_TNR,incident_Mach_shock_FPR_TNR) %ploting FPR<->TNR lim
hold on
plot(omegas_TNR_LSR(2:end),incident_Mach_shock_TNR_LSR)  %ploting TNR<->LSR lim

legend(legends,'Location','eastoutside')
xlabel("\omega (deg)")
ylabel("Msh")
title(['Computed boundaries for slow-fast ' name_I '->' name_II ' refraction'])
xlim([15,90])
%ylim([0,1])
hold off

function val=xiToMach(xi)
    val=sqrt(xiToSqMach(xi,gamma_I,pi/2));
end