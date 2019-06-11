function M1sq=xiToSqMach(xi,gamma,phi)
    %computes square of pre-shock Mach when knowing...
        %... pressure ratio
    %Inputs:
        %xi: pressure ratio p2/p1
        %gamma: ratio of specific heats
        %phi: shock - pre-shock flow angle
    %Output:
        %M1sq: square of preshock-Mach number
    M1sq=1/(2*gamma*sin(phi)^2)*((gamma-1)+(1+gamma)*xi);
end