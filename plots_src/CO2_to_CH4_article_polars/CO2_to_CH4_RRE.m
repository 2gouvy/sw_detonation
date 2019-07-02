%%CO2->CH4 wave refraction with chi=.78...
    %%... and beta=31°. RRE system.
%%The same polars were plotted in the following article:...
    %%Abd-el-Fattah, Henderson (1978): Shock waves at...
    %%... slow-fast gas interface

%gas parameters
    %ratios of specific heat
    gamma_CO2=1.288; %slow material, first phase
    gamma_CH4=1.303; %fast material, second phase
    %molecular masses:
    mu_CO2=44.01;
    %mu_CH4=16.04; %for pure CH4
    mu_CH4=18.84; %for contaminated CH4

%incident shock properties
chi=.78; %(pre-shock pressure)/(post-shock pressure)
xi=1/chi; %corresponds to usual p2/p1 for stationary shock
Msh=sqrt(xiToSqMach(xi,gamma_CO2,pi/2)); %corresponding...
    %... Mach for incident shock
beta=31*pi/180; %shock-interface angle in rad

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

%adding legend and setting window limits
legend('Incident CO2 polar','Transmited CH4 polar',...
    'incident \xi','y axis')
title("Shock polars for CO2->CH4 refraction (s->f) with \beta=31° and \chi=.78")
xlabel('\delta (rad)')
ylabel('\xi')
xi_lim_CH4=xiLim(Mt,gamma_CH4);
ylim([1,1.1*xi_lim_CH4])
hold off