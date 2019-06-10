%%Simulating a fast-slow refraction of a shock going from CO2 to CH4

%gas parameters
gamma_CO2=1.288; %slow material, first phase
gamma_CH4=1.303; %fast material, second phase
%molecular masses:
    mu_CO2=44.01;
    mu_CH4=16.04;

%incident shock properties
Xi=.78; %(pre-shock pressure)/(post-shock pressure)
xi=1/Xi; %corresponds to usual p2/p1 for stationary shock
Msh=sqrt(xiToSqMach(xi,gamma_CO2,pi/2)); %corresponding...
    %... Mach for incident shock
beta=33.27*pi/180;
Mi=Msh/sin(beta);
Mt=sqrt((gamma_CO2*mu_CH4)/(gamma_CH4*mu_CO2))*Mi;

ny=2000;

%plot incident shock polar
plotPolar(Mi,gamma_CO2,ny);
%plot transmitted shock polar
plotPolar(Mt,gamma_CH4,ny);
xi_lim_CH4=xiLim(Mt,gamma_CH4);
delta_max=deltaMax(Mt,gamma_CH4);
semilogy([-1.*delta_max,1.1*delta_max],[xi,xi],'--')
hold on
deltai=atan(sqrt(tanDefSq(xi,Mi,gamma_CO2)));
Mr=sqrt(postShockMachSq(xi,Mi,gamma_CO2));
plotPolar(Mr,gamma_CO2,ny,2,xi,deltai)
legend('Incident CO2 polar','Transmited CH4 polar',...
    'incident \xi','y axis','Reflected shock polar')
title("Shock polars")
xlabel('\delta (rad)')
ylabel('\xi')
ylim([1,1.1*xi_lim_CH4])
hold off