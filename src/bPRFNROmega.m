function omega_deg=bPRFNROmega(xi,gamma_I,gamma_II,mu_I,mu_II)
%Gets cooresonding interface inclination at BPR-FNR boundary at...
  %... constant xi
%Inputs:
  %xi: incident pressure jump
  %gamma_I, gamma_II: ratios of specific heats for each phase
  %mu_I, mu_II: molecular weighs of each phase
%Outputs:
  %omega_deg: corresponding inclination angle in degrees
    Msh=sqrt(xiToSqMach(xi,gamma_I,pi/2));
    c=(gamma_II+1)/(gamma_I+1)*sqrt((gamma_I*mu_II)/(gamma_II*mu_I))...
        *(Msh^2-1)/Msh;
    Vi_d_Vt=2/(c+sqrt(c^2+4))*Msh*sqrt((gamma_I*mu_II)/(gamma_II*mu_I));
    omega_deg=180/pi*asin(Vi_d_Vt);
end
