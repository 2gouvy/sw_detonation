function [xis,deltas]=getPolar(M1,gamma,ny,varargin)
%Computes shock polar
%Inputs:
    %ny: number of values on y-axis
    %mode: first varargin
        %There are multiple modes for the function:
            %0, default: calcultes usual, polar
            %1: choosing xi_lim, instead of taking the usual maximum xi...
                %... for y-axis boundary
                %xi_lim=varargin{2}: new xi limit
            %2: calculating polar as second deviation without changing
                %prev_xi: pressure ratio from previous shock, if any,...
                    %...entered in second varargin
                %prev_dev: deviation in previous shock, if any,...
                    %...entered in third varargin
                %xi_lim: optional argument, 4th varargin...
                    %... new xi limit
%Outputs:
    %xis: xi-coordinates of polar points
    %deltas: delta-coordinates of polar points (radian)
gamma;
M1;
xi_lim=xiLim(M1,gamma); %max value of xi in shock polar

%identifying mode
mode=0;
if nargin>=4
    mode=varargin{1};
    if mode==1
        xi_lim=varargin{2};
    elseif mode==2
       prev_xi=varargin{2};
       prev_dev=varargin{3};
       if nargin==7;
           xi_lim=varargin{4};
       end
    end
end
if mode ~= 2
    prev_xi=1;
    prev_dev=0;
end

%[M1,gamma]
xi_log_plot_step=(xi_lim)^(1/ny);
xi_axis=xi_log_plot_step.^(0:ny);%log-scale xi points
xi_axis(end)=xi_lim;

%polar calculation
delta_pos=zeros(1,ny+1);
for i=1:ny+1
    delta_pos(i)=atan(sqrt(tanDefSq(xi_axis(i),M1,gamma)));
    %if tanDefSq(xi_axis(i),M1,gamma)<0
    %    [xi_axis(i),xi_lim,(xi_axis(i)==xi_lim)]
    %end
end %calculate corresponding deviation angles for xi
delta_neg=flip(-delta_pos); %calculate corresponding negative angles
deltas=[delta_neg,delta_pos]+prev_dev;
xis=prev_xi*[flip(xi_axis),xi_axis];
end