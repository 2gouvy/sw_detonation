function [xilim]=xi_lim(M1,gamma)
    xilim=(2*gamma)/(gamma+1)*(M1^2-1)+1;
end