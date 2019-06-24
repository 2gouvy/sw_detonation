function is_sf_RR=isBeforeSfRRRToBPR(xi,omega_rad,gamma_I,gamma_II,mu_I,mu_II,ny)
%determines whether refraction will be regular or not...
    %... does so by checking if a reflected wave polar would have points...
    %... inside transmitted wave polar
    %Inputs:
      %xi: incident pressure jump
      %omega_rad: interface inclination in radian
      %gamma_I, gamma_II: ratios of specific heats for each phase
      %mu_I, mu_II: molecular weighs of each phase
      %ny: number of points on xi axis when computing polars
%computing incident polar
Msh=sqrt(xiToSqMach(xi,gamma_I,pi/2));
Mi=Msh/sin(omega_rad);
[i_xis,i_deltas]=getPolar(Mi,gamma_I,ny);
%computing transmitted polar
Mt=sqrt(gamma_I*mu_II/(gamma_II*mu_I))*Mi;
[t_xis,t_deltas]=getPolar(Mt,gamma_II,ny);
max_t_xi=max(t_xis);
if max_t_xi<=xi
    is_sf_RR=false;%case if pressure jump is too high for transmitted shock
else
    %computing reflected polar
    Mr=sqrt(postShockMachSq(xi,Mi,gamma_I));
    deltai=atan(sqrt(tanDefSq(xi,Mi,gamma_I)));
    [r_xis,r_deltas]=getPolar(Mr,gamma_I,ny,2,xi,deltai,max_t_xi);
    %looking for points from reflected polar inside transmitted
    i=ny+1;
    t_xi=t_xis(i);t_delta=t_deltas(i);
    [xi_coors,delta_coors,ind]=getPolarPoint(r_xis,r_deltas,0,t_xi);
    r_xi=xi_coors(1);r_delta=delta_coors(1);
    while (r_xi~=-1) && (r_delta>t_delta)
        i=i+1;
        t_xi=t_xis(i);
        t_delta=t_deltas(i);
        [xi_coors,delta_coors,ind]=getPolarPoint(r_xis,r_deltas,0,t_xi);
        r_xi=xi_coors(1);
        r_delta=delta_coors(1);
    end
    if r_xi==-1
        is_sf_RR=false;
    else
        is_sf_RR=true;
    end
end
end
