function nu=prandtlMeyer(M,gamma)
%Prandtl-Meyer function for expansion Waves
%Inputs:
  %M: Mach number
  %gamma: ratio of constant heats
    nu=sqrt((gamma+1)/(gamma-1))*atan(sqrt((gamma-1)/(gamma+1)*(M^2-1)))...
        -atan(sqrt(M^2-1));
end
