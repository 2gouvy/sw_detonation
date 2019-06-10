M1=1.3;
gamma=1.4;
ny=200;
xi_lim=xiLim(M1,gamma);
xi1=1+(xi_lim-1)*(.1);
plotPolar(M1,gamma,ny)
delta1=atan(sqrt(tanDefSq(xi1,M1,gamma)));
M2=sqrt(postShockMachSq(xi1,M1,gamma));
plotPolar(M2,gamma,ny,2,xi1,delta1)
title("Shock polars")
xlabel('\delta (rad)')
ylabel('\xi')
hold on
hold off