function xi_at_max_def=xiDeltaMax(M1,gamma)
%calcultes xi for maximum delfection angle at given Mach
%Inputs:
    %M1: pre-shock Mach
    %gamma: ratio of specific heats
%Output:
    %xi_at_max: pressure ratio for maximum deflection angle
    xi_at_max_def=.5*M1^2-1+sqrt(4/(gamma+1)+...
        (2*(gamma-1))/(gamma+1)*M1^2+.25*M1^4);
end