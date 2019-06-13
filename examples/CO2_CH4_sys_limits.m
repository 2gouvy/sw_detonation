%Computes boundaries between different systems for CO2-CH4 system

%gas parameters
    %ratios of specific heat
    gamma_CO2=1.288; %slow material, first phase
    gamma_CH4=1.303; %fast material, second phase
    %molecular masses:
    mu_CO2=44.01;
    mu_CH4=16.04;

%plots regular refraction and not regular refraction fontier in delta-xi...
    %...plane
nchis=100;ndeltas=100;nys=100;%number of points along chi axis; along...
    %... delta axis; along nys axis
chis=linspace(0,1,nchis);deltas=linspace(30,40,ndeltas);
betas=pi/180*deltas;
regular_irregular_bound=[];
for i=1:nchis
    for j=1:ndeltas
        is_sf_RR=isSfRR(chis(i),betas(j),gamma_CO2,gamma_CH4,mu_CO2,mu_CH4,nys);
            %determines if refraction is regular or irregular
        if (j~=1) && (is_sf_RR~=previous_is_sf_RR)
            regular_irregular_bound=[regular_irregular_bound,[deltas(j);chis(i)]];
        end
        previous_is_sf_RR=is_sf_RR;
    end
end
hold on
plot(regular_irregular_bound(1,:),regular_irregular_bound(2,:))%plot

%ploting RRE<->RRR boundary
nysp=100;% number of points along chi axis
chis=linspace(.35,1,nysp);
deltas=zeros(1,nysp);
syms delta %solving f(xi,delta)=0 boundary for xi going from
    %...1/chis(1) to 1/chis(nysp)
for i=1:nysp
    chi=chis(end-i+1);
    eqn = rrrToRreBoundary(1/chi,...
        delta,gamma_I,gamma_II,M_I,M_II)==0;
    deltas(end-i+1)=vpasolve(eqn,delta,[0,pi/2]);
end
hold on
plot(180/pi*deltas,chis)% Plot
legend("Regular<->non-regular boundary","RRE<->RRR boundary")
xlabel("\delta (°)")
ylabel("\xi")
title("Computed boundaries for CO2->CH4 refraction")
hold off