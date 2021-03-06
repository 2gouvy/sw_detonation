%%Computing boundaries between different refraction systems
    %%...for CO2->CH4 slow-fast refraction  in...
    %%... interface angle- incident pressure jump plane.
    %%A similar graph can be found in Abd-El-Fattah's and L. F. ...
    %%... Henderson 1978 paper Shock Waves at slow-fast gas interface, fig
    %%.. 13.
    %%The computed boundaries are compared with the article's data and we
    %%... can observe an overall good agreement, except for the FPR-TNR...
    %%... boundary.

limits_computation_time=cputime;%initialising timer for computation time

%gas parameters
    %ratios of specific heat
    gamma_CO2=1.288; %slow material, first phase
    gamma_CH4=1.303; %fast material, second phase
    %molecular masses:
    mu_CO2=44.01;
    mu_CH4=16.04; %non-contaminated CH4
    %mu_CH4=18.84; %contaminated CH4

 
%Computing analytical resolution of RRE<->RRR boundary:
    %Solves analytical RRE<->RRR equation at different given omega values
    npts=100;%number of points along chi axis
    omegas_RRE_RRR=linspace(31,37,npts); %omega coordinates in degrees
    omegas_RRE_RRR_rad=pi/180*omegas_RRE_RRR; %omega coordinates in radian
    xis_RRE_RRR=zeros(1,npts); %prealocated xi values

    %solving xi-equation for each value of omega
    syms xi
    for i=1:npts
        omega_rad=omegas_RRE_RRR_rad(i);
        eqn = rrrToRreBoundary(xi,...
            omega_rad,gamma_CO2,gamma_CH4,mu_CO2,mu_CH4)==0;
            a=[vpasolve(eqn,xi,[1,Inf])];
        if ~(isempty(a))
            xis_RRE_RRR(i)=a;
        end
    end
    chis_RRE_RRR=1./xis_RRE_RRR;

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
            is_sf_RRE=isSfStrongRRE(1/chis_RRE(i),pi/180*next_omega,gamma_CO2,...
                gamma_CH4,mu_CO2,mu_CH4,nys_RRE);
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

    chis_RRR_BPR=linspace(.42,1,nchis);%chi coordinates
    omegas_RRR_BPR=zeros(1,nchis);%pre-allocated omega coordinates
    for i=1:nchis
        omega_inf=0;omega_sup=90;%initial dichotomy bounds
        while omega_sup-omega_inf>rrr_bpr_boundary_precision %dichotomy...
            %...termination criteria

            %executing dichotomy, with lower bound having true value and...
                %... upper bound having false value
            next_omega=(omega_sup+omega_inf)/2;
            if isBeforeSfRRRToBPR(1/chis_RRR_BPR(i),pi/180*next_omega,gamma_CO2,...
                    gamma_CH4,mu_CO2,mu_CH4,nys)
                omega_inf=next_omega;
            else
                omega_sup=next_omega;
            end
        end
        omegas_RRR_BPR(i)=omega_inf; %adding omega coordinate in degrees to...
            %...preallocated array
    end


%Computing BPR<->NFR boundary
    %Principle: Computes omegas for which Mt=1 at different values of chi
    ny_BPR_NFR=100;%number of points on chi axis
    chis_BPR_NFR=linspace(0,1,ny_BPR_NFR);%chi coordinates
    omegas_BPR_NFR=zeros(1,ny_BPR_NFR);%pre-allocated omega coordinates
    %calculating omega coordinates
    for i=1:ny_BPR_NFR
        omegas_BPR_NFR(i)=bPRFNROmega(1/chis_BPR_NFR(i),gamma_CO2,gamma_CH4...
            ,mu_CO2,mu_CH4);
    end

%Computes TNR<->LSR limit
    %Principle: Solves the Mr=1 equation
    npts=400; %number of points on chi axis
    chis_TNR_LSR=linspace(0,1,npts);
    omegas_TNR_LSR=zeros(1,npts);
    syms omega
    for i=2:npts
        %solving Mr=1 equation at specific omega
        xi=1/chis_TNR_LSR(i);
        eqn = postShockMachSq(xi,sqrt(xiToSqMach(xi,gamma_CO2,omega)),gamma_CO2)==1;
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
nchis_tnr=100; %number of points on chi axis
nys_tnr=30; %number of points in computed polars
chis_FPR_TNR=linspace(.34,1,nchis_tnr); %chi values
omegas_FPR_TNR=zeros(1,nchis_tnr);%pre-allocated omega array
Msj=0; %j-wave Mach initialization
fpr_tnr_boundary_precision=.01; %dichotomy terminatin threshold
for i=1:nchis_tnr
    solve_Msj_eq=true;
    omega_inf_tnr=0;omega_sup_tnr=90;
    while omega_sup_tnr-omega_inf_tnr>fpr_tnr_boundary_precision
        next_omega=(omega_sup_tnr+omega_inf_tnr)/2;
        [is_before_trans,Msj]=isBeforeFPRToTNR(1/chis_FPR_TNR(i),pi/180*next_omega,...
             gamma_CO2,gamma_CH4,mu_CO2,mu_CH4,nys_tnr,solve_Msj_eq,Msj);
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

%Recovering On shock wave refraction at slow-fast gas interface (1978)...
    %... by Abd-El-Fattah, L.F. Henderson fig 13 article data in text file
fid=fopen('article_syst_lims.txt','r');
read_content=false;
no_graph_lines=6;
points=zeros(2*no_graph_lines,50);line_index=0;
no_read_lines=0;
for graph_line_no=1:no_graph_lines
    no_read_points=0;
    file_line=fgetl(fid);
    while ~contains(file_line,"Line")
        file_line=fgetl(fid);
    end
    if read_content
        no_read_lines=no_read_lines+1;
        line_index=line_index+1;
        file_line=fgetl(fid);
        while ~strcmp('',file_line) && ischar(file_line)
            point=str2num(file_line);
            no_read_points=no_read_points+1;
            points(2*no_read_lines-1:2*no_read_lines,no_read_points+1)=...
                point';
            file_line=fgetl(fid);
        end
        points(2*no_read_lines-1,1)=no_read_points;
    else
        read_content=true;
    end
end

%Ploting limits
figure
hold on
plot(omegas_RRE_RRR,chis_RRE_RRR,'m--') %ploting RRE<->RRR limit
hold on
plot(omegas_RRE(2:end-1),chis_RRE(2:end-1),'m') %ploting RRE<->... limit
hold on
plot(omegas_RRR_BPR,chis_RRR_BPR,'r') %ploting RRR<->BPR limit
hold on
plot(omegas_BPR_NFR,chis_BPR_NFR,'g') %ploting BPR<->NFR lim
hold on
plot(omegas_FPR_TNR,chis_FPR_TNR,'c') %ploting FPR<->TNR lim
hold on
plot(omegas_TNR_LSR(2:end),chis_TNR_LSR(2:end),'b')  %ploting TNR<->LSR lim
cols=['m';'r';'g';'c';'b'];
for i=1:no_read_lines
    hold on
    no_read_points=points(2*i-1,1);
    plot(points(2*i-1,2:no_read_points+1),points(2*i,2:no_read_points+1),[cols(i) 'o']);
    %ploting article data
end

vw_art_ang=[27,32.06,33.27,34.49,38,46];
n_vw_art_ang=length(vw_art_ang);
scatter(vw_art_ang,.78*ones(1,n_vw_art_ang),'<')
w_art_ang=[50.5,55];
n_w_art_ang=length(w_art_ang);
scatter(w_art_ang,.53*ones(1,n_w_art_ang),'^')
s_art_ang=[30,46,58];
n_s_art_ang=length(s_art_ang);
scatter(s_art_ang,.18*ones(1,n_s_art_ang),'>')
ltxt=sprintf(...
"RRE->Intromission->RRR->Shock critical\n->BPR->FPR from Nourgaliev et al. fem sims");
legends={"RRE<->RRR (analytical resolution)",...
    "RRE->... (graphical resolution)",...
    "RRR<->BPR","BPR<->FNR","FPR<->TNR",...
    "TNR<->LSR",...
    "RRE<->... article boundary points",...
    "RRR<->BPR article boundary points",...
    "BPR<->FNR article boundary points",...
    "FPR<->TNR article boundary points",...
    "TNR<->LSR article boundary points",ltxt,...
    "TRR->TNR from Nourgaliev et al. sims",...
    "RRE->BPR->TMR from Nourgaliev et al. fems sims"};

legend(legends,'Location','eastoutside')
xlabel("$\omega_i$ (deg)",'interpreter','latex')
ylabel("\chi")
title("Computed boundaries for CO2->CH4 refraction")
xlim([25,90])
ylim([0,1])
hold off
