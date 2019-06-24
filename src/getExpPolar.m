function [xis_exp,deltas_exp]=getExpPolar(Me,gamma,ny_exp,prev_xi,prev_dev,varargin)
%gets expansion polar
%Inputs:
    %Me: pre-expansion Mach
    %gamma: ratio_of_constant_heats
    %ny_exp: number of computed points along prev_xi-aprev_xis
    %prev_xi: pressure jump before expansion
    %prev_dev: deviation before expansion
    %Optional argument (varargin):
        %min_xi: lower-xi bound
if nargin>5
    min_xi=varargin{1};
else
    min_xi=1;
end
low_bound=log(min_xi/prev_xi)/log(10);
xis_exp=logspace(low_bound,0,ny_exp);
deltas_exp=zeros(1,ny_exp);
for i=1:ny_exp
    M_post_exp=sqrt(2/(gamma-1)*((1+(gamma-1)/2*Me^2)*xis_exp(i)^...
        ((1-gamma)/gamma)-1));
    deltas_exp(i)=prandtlMeyer(M_post_exp,gamma)-...
        prandtlMeyer(Me,gamma);
end
xis_exp=xis_exp*prev_xi;
deltas_exp=deltas_exp+prev_dev;
end