%%Simulating a fast-slow refraction of a shock going from CO2 to CH4
    %%Similar ploting can be founf in Nourgaliew et al.: Shock wave...
        %%...shock wave refraction patterns at interfaces (2005), Fig 8 b)

%gas parameters
    %ratios of specific heat
    gamma_CO2=1.288; %slow material, first phase
    gamma_CH4=1.303; %fast material, second phase
    %molecular masses:
    mu_CO2=44.01;
    mu_CH4=16.04;

%incident shock properties
chi=.78; %(pre-shock pressure)/(post-shock pressure)
xi=1/chi; %corresponds to usual p2/p1 for stationary shock
Msh=sqrt(xiToSqMach(xi,gamma_CO2,pi/2)); %corresponding...
    %... Mach for incident shock
beta=33.27*pi/180; %shock-interface angle in rad

%Calculating transmitted shock Mach
Mi=Msh/sin(beta);
Mt=sqrt((gamma_CO2*mu_CH4)/(gamma_CH4*mu_CO2))*Mi;

%plot incident shock polar
ny=2000;
plotPolar(Mi,gamma_CO2,ny);
%plot transmitted shock polar
plotPolar(Mt,gamma_CH4,ny);

%plot incident wave pressure jump
deltai=atan(sqrt(tanDefSq(xi,Mi,gamma_CO2)));
semilogy([-1.1*deltai,1.1*deltai],[xi,xi],'--')

%plot reflected shock polar
hold on
Mr=sqrt(postShockMachSq(xi,Mi,gamma_CO2));
plotPolar(Mr,gamma_CO2,ny,2,xi,deltai)

%adding legend and setting window limits
legend('Incident CO2 polar','Transmited CH4 polar',...
    'incident \xi','y axis','Reflected shock polar')
title("Shock polars for CO2->CH4 refraction (s->f) with \beta=33.27 and \chi=.78")
xlabel('\delta (rad)')
ylabel('\xi')
xi_lim_CH4=xiLim(Mt,gamma_CH4);
ylim([1,1.1*xi_lim_CH4])
hold off