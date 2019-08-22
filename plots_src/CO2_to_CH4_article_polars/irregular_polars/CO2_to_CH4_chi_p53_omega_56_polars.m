
%gas parameters
    %names
        name_I='CO2';
        name_II='CH4';
    %ratios of constant heat
        gamma_I=1.288;
        gamma_II=1.301;
    %molecular masses
        mu_I=44.01;
        mu_II=18.84;

%Plotting parameters
ny=200;
%Params of figure 6a of slow-fast gas interface
chi=.53;
omega_deg=56;
xi_i=1/chi;
omega_rad=pi/180*omega_deg;

Msi=sqrt(xiToSqMach(xi_i,gamma_I,pi/2)); %i shock Mach
M1i=Msi/sin(omega_rad); %i free-stream Mach
%Computing Mst:
c=(gamma_II+1)/(gamma_I+1)*sqrt((gamma_I*mu_II)/(gamma_II*mu_I))...
    *(Msi^2-1)/Msi;
Mst=1/2*(c+sqrt(c^2+4)); %t shock Mach
xi_t=((1-gamma_II)+2*gamma_II*Mst^2)/(1+gamma_II); %transmitted shock...
    %... pressure jump
%computing Msj, j shock Mach:
syms Msj
eqn= ((1-gamma_I)+2*gamma_I*Msj^2)/(1+gamma_I)*...
    (1-(gamma_II-1)/(gamma_I+1)*sqrt((gamma_I*mu_II)/(gamma_II*mu_I))*...
    (Msj-1/Msj))^(-2*gamma_II/(gamma_II-1))==xi_t;
a=[vpasolve(eqn,Msj,[1,Inf])];
if ~(isempty(a)) %checking if slution was found
    figure
    hold on
    plotPolar(M1i,gamma_I,ny); %plotting i polar
    hold on
    plotPolar(Mst,gamma_II,ny); %plotting t polar
    hold on
    %plotting reflected polar:
    Mr=sqrt(postShockMachSq(xi_i,M1i,gamma_I)); %reflected shock free-...
        %... stream Mach
    i_dev=atan(sqrt(tanDefSq(xi_i,M1i,gamma_I))); %i deviation
    plot([-1.1,1.1]*i_dev*180/pi,xi_i*[1,1],'--');
    plotPolar(Mr,gamma_I,ny,2,xi_i,i_dev); %plotting r polar
    hold on
    Msj=a(1); %j shock Mach
    xi_j=((1-gamma_I)+2*gamma_I*Msj^2)/(1+gamma_I); %j pressure jump
    M1k=sqrt(postShockMachSq(xi_j,M1i,gamma_I)); %k shock Mach
    j_dev=atan(sqrt(tanDefSq(xi_j,M1i,gamma_I))); %j deviation
    plotPolar(M1k,gamma_II,ny,2,xi_j,-j_dev); %plotting k polar
    
    fid=fopen('prec_pol3.txt','r');
    extract_plots;
    semilogy(points(1,2:points(1,1)+1),points(2,2:points(1,1)+1),'mo')
    semilogy(points(3,2:points(3,1)+1),points(4,2:points(3,1)+1),'go')
    
    legend("Incident and j precursor shocks","Transmitted Shock",...
        "Incident \xi","Reflected shock", "k shock","Article r polar",...
        "Article k polar");
    title(['Polars for t-precursor wave ' name_I '->'...
         name_II ' refraction with \chi=' num2str(chi) ' and \omega='...
        num2str(omega_deg) ' deg'])
    xlabel('\delta (rad)')
    set(gca,'yscale','log')
    ylabel('\xi')
    hold off
else
    'Could not compute j shock wave Mach!'
end
