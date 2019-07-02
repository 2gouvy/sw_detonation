function is_sf_strong_RRE=isSfStrongRRE(xi,omega_rad,gamma_I,gamma_II,...
        mu_I,mu_II,ny_exp,varargin)
    %Determines if system is RRE with graphical resolution
    %Inputs:
      %xi: incident pressure jump
      %omega_rad: interface inclination in radian
      %gamma_I, gamma_II: ratios of specific heats for each phase
      %mu_I, mu_II: molecular weighs of each phase
      %ny_exp: number of points on xi axis when computing polars
      %varargin{1}: temperature ratio T_I/T_II
    %Output:
      %is_sf_strong_RRE: bool. true if RRE, false if not
    temp_ratio=1;
    if nargin==8
        temp_ratio=varargin{1};
    end
    Msh=sqrt(xiToSqMach(xi,gamma_I,pi/2)); %computing incident shock Mach
    Mi=Msh/sin(omega_rad); %computing incident free-stream Mach
    Mt=sqrt(temp_ratio*(gamma_I*mu_II)/(gamma_II*mu_I))*Mi; %transmitted free-stream Mach
    t_xi_delta_max=xiDeltaMax(Mt,gamma_II); %xi value at maximum deviation for ...
      %... transmitted wave polar
    t_delta_max=atan(sqrt(tanDefSq(t_xi_delta_max,Mt,gamma_II))); %maximum deviation...
      %... for transmitted polar
    t_xi_delta=atan(sqrt(tanDefSq(xi,Mt,gamma_II))); %deviation at intersection...
      %... between transmitted polar and incident xi
    prev_delta=atan(sqrt(tanDefSq(xi,Mi,gamma_I))); %deviation after incident shock
    t_delta_max<prev_delta;
    Me=sqrt(postShockMachSq(xi,Mi,gamma_I)); %Mach number behind incident shock...
      %... and before supposed expansion wave
    if Me<1 %no expansion if Mach subsonic
        is_sf_strong_RRE=false;
    %executing graphical method
    elseif t_delta_max<prev_delta
        is_sf_strong_RRE=false;
    elseif prev_delta<=t_xi_delta
        is_sf_strong_RRE=true;
    elseif (prev_delta>t_xi_delta) && (xi<t_xi_delta_max)
        is_sf_strong_RRE=false;
    else
        [t_xis,t_deltas]=getPolar(Mt,gamma_II,ny_exp,3,t_xi_delta_max);
        %Me=sqrt(postShockMachSq(xi,Mi,gamma_I));
        [r_xis,r_deltas]=getExpPolar(Me,gamma_I,ny_exp,xi,prev_delta,...
            t_xi_delta_max);
%         figure
%         hold on
%         plot(r_deltas,r_xis)
%         hold on
%         plot(t_deltas,t_xis)
%         hold off
%         legend("transmitted","reflected")
        i=2*(ny_exp+1);
        t_xi=t_xis(i);t_delta=t_deltas(i);
        [xi_coors,delta_coors,ind]=getPolarPoint(r_xis,r_deltas,0,t_xi,true);
        r_xi=xi_coors(1);r_delta=delta_coors(1);
        (r_delta>t_delta);
        while (i>ny_exp+1) && (r_delta>t_delta) && (r_xi~=-1)
            i=i-1;
            t_xi=t_xis(i);
            t_delta=t_deltas(i);
            [xi_coors,delta_coors,ind]=getPolarPoint(r_xis,r_deltas,...
                0,t_xi,true);
            r_xi=xi_coors(1);
            r_delta=delta_coors(1);
        end
        if (i==ny_exp+1)||(r_xi==-1)
            is_sf_strong_RRE=false;
        else
            is_sf_strong_RRE=true;
        end
    end
end
