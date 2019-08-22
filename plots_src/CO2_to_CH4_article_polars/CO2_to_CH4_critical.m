%%CO2->CH4 wave refraction with chi=.78...
    %%... and angle of incidence omega=34.4885�.
%%The same polars were plotted in log-scale in article:...
    %%Nourgaliev, Sushchikh, Dinh, Theofanous (2005)

%gas parameters
    %names
    name_I='CO2';
    name_II='CH4';
    %ratios of specific heat
    gamma_I=1.288; %slow material, first phase
    gamma_II=1.303; %fast material, second phase
    %molecular masses:
    mu_I=44.01;
    mu_II=16.04;

%incident shock properties
figure
chi=.78; %(pre-shock pressure)/(post-shock pressure)
xi=1/chi; %corresponds to usual p2/p1 for stationary shock
Msh=sqrt(xiToSqMach(xi,gamma_I,pi/2)); %corresponding...
    %... Mach for incident shock
omega_deg=34.4885; %shock-interface angle in deg
omega_rad=omega_deg*pi/180; %shock-interface angle in rad

%Calculating transmitted shock Mach
Mi=Msh/sin(omega_rad);
Mt=sqrt((gamma_I*mu_II)/(gamma_II*mu_I))*Mi;

%plot incident shock polar
ny=2000;
plotPolar(Mi,gamma_I,ny);
%plot transmitted shock polar
plotPolar(Mt,gamma_II,ny);

%plot incident wave pressure jump
deltai=atan(sqrt(tanDefSq(xi,Mi,gamma_I)));
semilogy(deltai*180/pi*[-1.1,1.1],[xi,xi],'--')

%plot reflected shock polar
hold on
Mr=sqrt(postShockMachSq(xi,Mi,gamma_I));
plotPolar(Mr,gamma_I,ny,2,xi,deltai);

%adding legend and setting window limits
legend('Incident CO2 polar','Transmited CH4 polar',...
    'incident \xi','y axis','Reflected shock polar')
title(['Polars for ' name_I '->'...
         name_II ' refraction with \chi=' num2str(chi) ' and \omega='...
        num2str(omega_deg) ' deg'])
xlabel('\delta (deg)')
ylabel('\xi')
xi_lim_II=xiLim(Mt,gamma_II);
ylim([1,1.1*xi_lim_II])
delta_max=deltaMax(Mt,gamma_II);
xlim(delta_max*180/pi*1.3*[-1,1])
set(gca,'yscale','log')
hold off
