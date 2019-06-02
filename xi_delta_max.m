function xideltam=xi_delta_max(M1,gamma)
%calculates maximum xi for given Mach number
    xideltam=.5*M1^2-1+sqrt(4/(gamma+1)+...
        (2*(gamma-1))/(gamma+1)*M1^2+.25*M1^4);
end