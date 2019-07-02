%  Using Chemkin fits to calculate H2_O2 mixuture gamma

%  Fits where computied using the following database:
    %  Alexander Burcat and Branko Ruscic
    %  Ideal Gas Thermochemical Database with updates from Active...
        %... Thermochemical Tables
    % 
    %  <http://garfield.chem.elte.hu/Burcat/burcat.html>; 9 January 2013

a_H2=[2.34433112E+00,7.98052075E-03,...
    -1.94781510E-05,2.01572094E-08,...
    -7.37611761E-12,-9.17935173E+02,...
    6.83010238E-01]; %200-6000K fit
a_O2=[3.78245636E+00,-2.99673416E-03,...
    9.84730201E-06,-9.68129509E-09,...
    3.24372837E-12,-1.06394356E+03,...
    3.65767573E+00]; %200-6000K fit
R=8.314; %gas constant
T=298; %temp
x_H2=2/3; x_O2=1/3; %gas proportions
phase_heat_ratio_fit=x_H2*standardHeatRatioPoly(T,a_H2)...
    +x_O2*standardHeatRatioPoly(T,a_O2); %heat ratio for mix

%gamma:
phase_gamma=phase_heat_ratio_fit/...
    (phase_heat_ratio_fit-1)

%Gibbs free energy
P=1e5;P_0=1e5;
a_reac=x_H2*a_H2+x_O2*a_O2;
standard_free_enthalpy_poly=R*T*(standardEnthalpyPoly(T,a_reac)...
    -standardEntropyPoly(T,a_reac));
free_enthalpy_=standard_free_enthalpy_poly+...
    R*T*log(P_0/((1/3)^.5*(2/3)*P));


function val=standardHeatRatioPoly(T,a)
    val=0;
    for k=1:5
        val=val+a(k)*T^(k-1);
    end
end

function val=standardEnthalpyPoly(T,a)
    val=0;
    for k=1:5
        val=val+(a(k)/k)*T^(k-1);
    end
    val=val+a(6)/T;
end

function val=standardEntropyPoly(T,a)
    val=a(1)*log(T);
    for k=2:5
        val=val+a(k)*T^(k-1);
    end
    val=val+a(7);
end