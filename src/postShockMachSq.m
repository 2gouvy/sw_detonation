function M2_sqr=postShockMachSq(xi,M1,gamma)
    %computes mach after shock, having pressure jump and...
        %...Mach before shock
    %Input:
        %xi: pressure jump
        %M1: pre-shock Mach
        %gamma: ratio of specific heats
    %Output:
        %M2_sqr: square of post-shock Mach
    M2_sqr=(M1^2*[(gamma+1)*xi+(gamma-1)]-2*(xi^2-1))...
        /(xi*[(gamma-1)*xi+(gamma+1)]);
end
