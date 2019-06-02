M1=1.3;
gamma=1.4;
ny=200;
xilim=xi_lim(M1,gamma);
xi1=1+(xilim-1)*(.1)
plot_polar(M1,gamma,ny)
delta1=atan(sqrt(tan_def_sq(xi1,M1,gamma)));
M2=sqrt(new_sq_mach(xi1,M1,gamma));
plot_polar(M2,gamma,ny,xi1,delta1)
title("Shock polars")
xlabel('\delta')
ylabel('\xi')
hold on
hold off