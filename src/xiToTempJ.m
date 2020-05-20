function temp_jump=xiToTempJ(xi,gamma)
    temp_jump=xi*((gamma-1)*xi+(gamma+1))/((gamma+1)*xi+(gamma-1));
end