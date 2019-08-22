npts=40;

a_H2O_low=[0.41986352E+01,-0.20364017E-02,0.65203416E-05,...
    -0.54879269E-08,0.17719680E-11,-0.30293726E+05,...
    -0.84900901E+00]; %200-1000K fit
a_H2O_high=[0.26770389E+01,0.29731816E-02,-0.77376889E-06,...
    0.94433514E-10,-0.42689991E-14,-0.29885894E+05,...
    0.68825500E+01]; %1000-6000K fit
a_H2_low=[2.34433112E+00,7.98052075E-03,...
    -1.94781510E-05,2.01572094E-08,...
    -7.37611761E-12,-9.17935173E+02,...
    6.83010238E-01]; %200-1000K fit
a_H2_high=[2.93286575E+00,8.26608026E-04,-1.46402364E-07,...
    1.54100414E-11,-6.88804800E-16,-8.13065581E+02,...
    -1.02432865E+00]; %1000K-6000K
a_O2_low=[3.78245636E+00,-2.99673416E-03,...
    9.84730201E-06,-9.68129509E-09,...
    3.24372837E-12,-1.06394356E+03,...
    3.65767573E+00]; %200-1000K fit
a_O2_high=[3.66096065E+00,6.56365811E-04,-1.41149627E-07,...
    2.05797935E-11,-1.29913436E-15,-1.21597718E+03,...
    3.41536279E+00]; %1000K-6000K
xi_mws=linspace(0,2/3,npts);
f_vals=arrayfun(@f,xi_mws);
plot(xi_mws,f_vals);
press_rats=linspace(0,1,npts);
temp_rats=linspace(298/6000,1,npts);
xis_mw=zeros(npts-1,npts-1);
for j=2:npts
    T=298/temp_rats(j);
    if T>1000
        a_H2=a_H2_high;a_O2=a_O2_high;a_H2O=a_H2O_high;
    else
        a_H2=a_H2_low;a_O2=a_O2_low;a_H2O=a_H2O_low;
    end
    H2_free_energy_term=standardEnthalpyPoly(T,a_H2)...
        -standardEntropyPoly(T,a_H2);
    O2_free_energy_term=standardEnthalpyPoly(T,a_O2)...
        -standardEntropyPoly(T,a_O2);
    H2O_free_energy_term=standardEnthalpyPoly(T,a_H2O)...
        -standardEntropyPoly(T,a_H2O);
    reaction_free_energy_term=H2O_free_energy_term-.5*O2_free_energy_term...
        -H2_free_energy_term;
    K0=exp(-reaction_free_energy_term);
    for i=2:npts
        '2d - plot computing indexes:'
        [i,j]
        syms xi_mw
        eqn = f(xi_mw)*sqrt(press_rats(i)) == K0;
        a=[vpasolve(eqn,xi_mw,[0,2/3])];
        if ~isempty(a)
            'Found solution';
            xis_mw(i-1,end-(j-1)+1)=a(1);
        end
    end
end
figure
x=[press_rats(2),press_rats(end)];
y=[temp_rats(2),temp_rats(end)];
imagesc(x,y,xis_mw)
colorbar
set(gca,'YDir','normal')
title({'$\xi_{mw}(\frac{T_{amb}}{T},\frac{P_{atm}}{P})$'},'interpreter','latex')
xlabel({'$\frac{P_{atm}}{P}$'},'interpreter','latex')
ylabel({'$\frac{T_{amb}}{T}$'},'interpreter','latex')

%plotting temperature jump for different values of chi
npts_temp=1000;
chis=linspace(0,1,npts_temp);
shock_temp_jumps=arrayfun(@shockTempJump,1./chis(2:end));
hold on
plot(chis(2:end),1./shock_temp_jumps,'r')
hold off

K0s=zeros(1,npts);
for j=1:npts
    T=298/temp_rats(j);
    if T>1000
        a_H2=a_H2_high;a_O2=a_O2_high;a_H2O=a_H2O_high;
    else
        a_H2=a_H2_low;a_O2=a_O2_low;a_H2O=a_H2O_low;
    end
    H2_free_energy_term=standardEnthalpyPoly(T,a_H2)...
        -standardEntropyPoly(T,a_H2);
    O2_free_energy_term=standardEnthalpyPoly(T,a_O2)...
        -standardEntropyPoly(T,a_O2);
    H2O_free_energy_term=standardEnthalpyPoly(T,a_H2O)...
        -standardEntropyPoly(T,a_H2O);
    reaction_free_energy_term=H2O_free_energy_term-.5*O2_free_energy_term...
        -H2_free_energy_term;
    K0=exp(-reaction_free_energy_term);
    K0s(j)=K0;
end
figure
plot(298./temp_rats,K0s)
ylabel({'$K^o(T)$'},'interpreter','latex')
xlabel({'$T$ $(K)$'},'interpreter','latex')
title({'$K^o(T)$'},'interpreter','latex')

% plotting xi_wm behind incident shock
npts_incident_shock=100;
chis=linspace(0,1,npts_incident_shock);
K0s=zeros(1,npts_incident_shock-1);
syms xi_wm
for i=2:npts_incident_shock
    i;
    xi=1/chis(i);
    T=298*shockTempJump(xi);
    if T>1000
        a_H2=a_H2_high;a_O2=a_O2_high;a_H2O=a_H2O_high;
    else
        a_H2=a_H2_low;a_O2=a_O2_low;a_H2O=a_H2O_low;
    end
    H2_free_energy_term=standardEnthalpyPoly(T,a_H2)...
        -standardEntropyPoly(T,a_H2);
    O2_free_energy_term=standardEnthalpyPoly(T,a_O2)...
        -standardEntropyPoly(T,a_O2);
    H2O_free_energy_term=standardEnthalpyPoly(T,a_H2O)...
        -standardEntropyPoly(T,a_H2O);
    reaction_free_energy_term=H2O_free_energy_term-.5*O2_free_energy_term...
        -H2_free_energy_term;
    K0=exp(-reaction_free_energy_term);
    K0s(i-1)=K0;
end
figure
plot(chis(2:end),K0s)
xlabel({'Incident shock $\chi$'},'interpreter','latex')
ylabel({'$K^o(T)$'},'interpreter','latex')
title({'$K^o(T)$ for incident shock'},'interpreter','latex')

function val=f(xi_mw)
    val=(xi_mw*(1-.5*xi_mw)^.5)/((2/3-xi_mw)*(1/3-1/2*xi_mw)^.5);
end

function val=standardEnthalpyPoly(T,a)
    val=0;
    for k=1:5
        val=val+(a(k)/k)*T^(k-1);
    end
    val=val+a(6)/T;
end

function val=standardEntropyPoly(T,a)
    val=a(1)*log(T);
    for k=2:5
        val=val+a(k)*T^(k-1);
    end
    val=val+a(7);
end

function shock_temp_jump=shockTempJump(xi)
    gamma=1.4016;
    shock_temp_jump=xi*((gamma-1)*xi+(gamma+1))/((gamma+1)*xi+(gamma-1));
end