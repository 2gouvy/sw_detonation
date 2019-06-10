function M2_sqr=postShockMachSq(xi,M1,gamma)
    M2_sqr=(M1^2*[(gamma+1)*xi+(gamma-1)]-2*(xi^2-1))...
        /(xi*[(gamma-1)*xi+(gamma+1)]);
end
