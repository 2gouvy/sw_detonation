function [t2_delta]=tanDefSq(xi,M1,gamma)
%inputs:
    %xi=P2/P1
    %M1=mach before shock
    %gamma= ratio of specific heats
%output:
    %t2_delta: square tangeant of deviation
    t2_delta=(xi-1)^2/(gamma*M1^2-(xi-1))^2*...
        [M1^2-(1+(gamma+1)/(2*gamma)*(xi-1))]/...
        [1+(gamma+1)/(2*gamma)*(xi-1)];
end
