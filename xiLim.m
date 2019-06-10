function [xi_lim]=xiLim(M1,gamma)
%calculates maximum xi for given Mach number
%Inputs:
    %M1: pre-shock Mach
    %gamma: ratio of specific heats
%Outputs:
    %xi_lim: maximum pressure ratio for M1
    xi_lim=(2*gamma)/(gamma+1)*(M1^2-1)+1;
end
