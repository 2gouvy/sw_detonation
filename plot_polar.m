function plot_polar(M1,gamma,ny,varargin)
%Plots shock polar
%Inputs:
    %ny: number of values on y-axis
    %prev_xi: pressure ratio from previous shock, if any,
        %entered in first varargin
    %prev_dev: deviation in previous shock, if any,
        %entered in second varargin
gamma
M1
xilim=xi_lim(M1,gamma) %max value of xi in shock
        %polar
if nargin==5
    prev_xi=varargin{1}
    prev_dev=varargin{2}
    hold on
    plot([0,0],[1,2],'--')
else
    prev_xi=1;
    prev_dev=0;
end
%[M1,gamma]
xi_log_plot_step=(xilim)^(1/ny);
xi_axis=xi_log_plot_step.^(0:ny);%log-scale xi points
xi_axis(end)=xilim;

%polar plot
delta_pos=zeros(1,ny+1);
for i=1:ny+1
    delta_pos(i)=atan(sqrt(tan_def_sq(xi_axis(i),M1,gamma)));
end %calculate corresponding deviation angles for xi
delta_neg=flip(-delta_pos); %calculate corresponding negative angles
delta_plot=[delta_neg,delta_pos]+prev_dev;
xi_plot=prev_xi*[flip(xi_axis),xi_axis];
hold on
semilogy(delta_plot,xi_plot)% plot full polar
end

