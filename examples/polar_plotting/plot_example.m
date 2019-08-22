%flow parameters
gamma=1.4;
M1=1.3;

%limits of plot
xi_lim=xiLim(M1,gamma); %max value of xi in shock polar
ny=100;%number of xi points
xi_log_plot_step=(xi_lim)^(1/ny);
xi_axis=xi_log_plot_step.^(0:ny);%log-scale xi points

%polar plot
delta_pos=zeros(1,ny+1);
for i=1:ny+1
    delta_pos(i)=atan(sqrt(tanDefSq(xi_axis(i),M1,gamma)));
end %calculate corresponding deviation angles for xi
delta_neg=flip(-delta_pos); %calculate corresponding negative angles
delta_plot=[delta_neg,delta_pos];
xi_plot=[flip(xi_axis),xi_axis];
semilogy(180/pi*delta_plot,xi_plot)% plot full polar

