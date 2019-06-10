function M1sq=xiToSqMach(xi,gamma,phi)
    M1sq=1/(2*gamma*sin(phi)^2)*((gamma-1)+(1+gamma)*xi);
end