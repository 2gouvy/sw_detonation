function [xi_lim]=xiLim(M1,gamma)
    xi_lim=(2*gamma)/(gamma+1)*(M1^2-1)+1;
end
