function delta_max=deltaMax(M1,gamma)
    %computes maximum possible deflection angle...
        %... at specific pre-shock Mach
    %Input:
        %M1: pre-shock Mach
        %gamma: ratio of specific heats
    %Output:
        %delta_max: maximum deflection angle in radian
    xi=xiDeltaMax(M1,gamma);
    delta_max=sqrt(tanDefSq(xi,M1,gamma));
end
