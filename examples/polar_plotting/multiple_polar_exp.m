%%Plot polar for regular reflection in a homgeneous...
    %%... perfect gas flow
    %%Theory explained in Zhaoyan Han and Xiezhen Yin: ...
        %%... Shock Dynamics (1993) p237-239
    

M1=1.3; %pre-shock Mach
gamma=1.4; %gas ratio of constant heat
ny=200; %number of points parameter
xi_lim=xiLim(M1,gamma); %maximum pressure ratio at M1
xi1=1+(xi_lim-1)*(.1);
plotPolar(M1,gamma,ny); %plotting incident polar
delta1=atan(sqrt(tanDefSq(xi1,M1,gamma))); %deflection...
    %... after incident shock
M2=sqrt(postShockMachSq(xi1,M1,gamma)); %Mach post-incident
plotPolar(M2,gamma,ny,2,xi1,delta1); %reflected polar

%plot description and legend
title("Shock polars")
xlabel('\delta (deg)')
ylabel('\xi')
hold on
hold off