function val=rrrToRreBoundary(xi,omega_rad,gamma_I,gamma_II,mu_I,mu_II)
%Computes RRE<->RRR boundary function, equals zeros at bounday
%Inputs:
  %xi: incident pressure jump
  %omega_rad: interface inclination in radian
  %gamma_I, gamma_II: ratios of specific heats for each phase
  %mu_I, mu_II: molecular weighs of each phase
%Outputs:
  %val: value of boundary function, is zeros ar boundary

Msh=sqrt(xiToSqMach(xi,gamma_I,pi/2)); %incident shock Mach
Mi=Msh/sin(omega_rad); %free-stream Mach
Mt=Mi*sqrt((gamma_I*mu_II)/(gamma_II*mu_I));
val=(1+gamma_I*Mi^2-xi)^2*((gamma_I-1)/(gamma_I+1)+xi)/...
    ((2*gamma_I)/(gamma_I+1)*Mi^2-(gamma_I-1)/(gamma_I+1)-xi)...
    -(1+gamma_II*Mt^2-xi)^2*((gamma_II-1)/(gamma_II+1)+xi)/...
    ((2*gamma_II)/(gamma_II+1)*Mt^2-(gamma_II-1)/(gamma_II+1)-xi);%value...
        %... function
end
