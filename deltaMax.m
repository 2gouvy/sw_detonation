function delta_max=deltaMax(M1,gamma)
    xi=xiDeltaMax(M1,gamma);
    delta_max=sqrt(tanDefSq(xi,M1,gamma));
end
