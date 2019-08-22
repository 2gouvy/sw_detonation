function [before_transition,Msj,Mst,Msi]=isBeforeFPRToTNR(xi_i,omega_rad,gamma_I,...
    gamma_II,mu_I,mu_II,ny,varargin)
    %Determines if FPR to TNR transition was made
    %Inputs:
        %xi_i: incident shock wave pressure jump
        %omega_rad: angle of incidence (rad)
        %gamma_I, gamma_II: ratios of constant heat
        %mus: molecular masses
        %ny: xi-axis resolution when getting polars
        %vargain{1}, solve_Msj_eq: bool determining if Msj 
            %... equation needs to be solved or not, used to reduce...
            %... computation time
        %varargin{2}, Msj: Msj value, used if solve_Msj_eq==false
        %varargin{3}, temp_ratio: ratio of temperatures=T_I/T_II
    solve_Msj_eq=true;
    temp_ratio=1;
    if nargin>=8
        solve_Msj_eq=varargin{1};
        if nargin>=10
            temp_ratio=varargin{3};
        end
    end
    Msi=sqrt(xiToSqMach(xi_i,gamma_I,pi/2)); %incident shock Mach
    Mi=Msi/sin(omega_rad); %incident free stream Mach

    c=(gamma_II+1)/(gamma_I+1)*sqrt(temp_ratio*(gamma_I*mu_II)...
        /(gamma_II*mu_I))*(Msi^2-1)/Msi;
    Mst=1/2*(c+sqrt(c^2+4)); %transmitted shock Mach
    xi_t=((1-gamma_II)+2*gamma_II*Mst^2)/(1+gamma_II); %transmitted pressure...
        %...jump

    %computing j precursor shock Mach
    if solve_Msj_eq
        syms Msj
        eqn= ((1-gamma_I)+2*gamma_I*Msj^2)/(1+gamma_I)...
            ==xi_t*(1-(gamma_II-1)/(gamma_I+1)...
                *sqrt((gamma_I*mu_II)/(gamma_II*mu_I))*...
            (Msj-1/Msj))^(2*gamma_II/(gamma_II-1));
        a=[vpasolve(eqn,Msj,[1,Inf])];
    else
        a=[varargin{2}];
    end
    if isempty(a)
        before_transition=true;
    else
        i_dev=atan(sqrt(tanDefSq(xi_i,Mi,gamma_I)));
        M1r=sqrt(postShockMachSq(xi_i,Mi,gamma_I));
        min_r_dev=i_dev-deltaMax(M1r,gamma_I); %minimum of r deviation
        Msj=a(1);
        xi_j=((1-gamma_I)+2*gamma_I*Msj^2)/(1+gamma_I);
        j_dev=atan(sqrt(tanDefSq(xi_j,Mi,gamma_I)));
        M1k=sqrt(postShockMachSq(xi_j,Mi,gamma_I));
        max_k_dev=-j_dev+deltaMax(M1k,gamma_I); %max of k deviation
        if max_k_dev<min_r_dev
            before_transition=false;
        else
            [xis_k,deltas_k]=getPolar(M1k,gamma_I,ny,2,xi_j,j_dev);
            [xis_r,deltas_r]=getPolar(M1r,gamma_I,ny,2,xi_i,i_dev);
            i=ny+2;
            xi_k=xis_k(i);delta_k=deltas_k(i);
            [xi_coors,delta_coors,inds]=getPolarPoint(xis_r,deltas_r,0,xi_k);
            xi_r=xi_coors(1);delta_r=delta_coors(1);
            while (i<=2*ny+1) && ((xi_r==-1) || (delta_k<delta_r))...
                    && (xi_k>=xi_i)
                i=i+1;
                xi_k=xis_k(i);
                delta_k=deltas_k(i);
                [xi_coors,delta_coors,inds]=...
                    getPolarPoint(xis_r,deltas_r,0,xi_k);
                xi_r=xi_coors(1);delta_r=delta_coors(1);
            end
            if (i==2*ny+1) || (xi_k<xi_i)
                before_transition=false;
            else
                before_transition=true;
            end
        end
    end
end