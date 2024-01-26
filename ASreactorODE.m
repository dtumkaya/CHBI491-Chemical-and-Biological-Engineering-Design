% Simulation of AS with Platinum catalyst
% TO DO: Correct the relation between the inputs and outputs

function dydV = ASreactorODE(V, y, R, V_react, U, A_surface, d_h)

% Extract variables
C_H2O = y(1);
C_O2 = y(2);
C_N2 = y(3);
C_NH3 = y(4);
T = y(5);

% R: NH3 oxidation
% 4 NH3 + 3 O2  -> 2 N2 + 6 H2O
a_c = 4/d_h; % 1/m
Sh = 3;
D_ab_NH3 = 1.215 *(10^-9)*T^1.7389; % m2/s
k_m = Sh*D_ab_NH3/d_h; % m/s
t_w = 5*10^-5; % m  

D_eff_NH3 = 100*1.215 *(10^-9)*T^1.7389; % m2/s
A_NH3 = 1.72 * 10^15;
Ea_NH3 = 99500; % J/mol

k_s_NH3 = A_NH3 * exp(-(Ea_NH3/(R*T)));
tiele_NH3 = t_w * sqrt(k_s_NH3/D_eff_NH3);
eff_n_NH3 = tanh(tiele_NH3)/tiele_NH3;

C_NH3_w = (k_m * a_c * C_NH3) / (eff_n_NH3 * k_s_NH3 + k_m * a_c); % (mol/m^3*s) / (1/s) -> mol/m^3
r_NH3 = - k_s_NH3 * C_NH3_w*0.5;

C_tot = 1/(8.21*10^(-5) * (T + 273.15)); % mol/m3

MW_air = 28.96; % g/mol
mass_flow = 1300000/3600; % g/s
F_tot = (mass_flow/MW_air); % mol/s
vflow = F_tot/C_tot; % m^3/s

dC_NH3_dV = r_NH3/vflow;
dC_N2_dV = -0.5*r_NH3/vflow;

F_H2O = C_H2O * vflow; 
F_O2 = C_O2 * vflow; 
F_N2 = C_N2 * vflow; 
F_NH3 = C_NH3 * vflow; 

T_ref = 298.15; % in K

% Specific heat capacities - Shomate equation
% https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Mask=1&Type=JANAFG&Table=on
Cp_NH3 = 19.996 + 49.77*(T/1000) - 15.38*(T/1000)^2 + 1.92 *(T/1000)^3 + 0.19/(T/1000)^2;
Cp_N2 = 19.51 + 19.89*(T/1000) - 8.60*(T/1000)^2 + 1.37 *(T/1000)^3 + 0.53/(T/1000)^2;
Cp_O2 = 31.32 - 20.24*(T/1000) + 57.87*(T/1000)^2 - 36.5 *(T/1000)^3 - 0.007374/(T/1000)^2;
Cp_H2O = 30.09 + 6.83*(T/1000) +  6.79*(T/1000)^2 - 2.53 *(T/1000)^3 + 0.082/(T/1000)^2;

% No! T bağlantılı ekle (ya da 300 derece)
delta_H_NH3ox = -1266000;  % joule/mol

dT_dV = ((U*A_surface*(T_ref-T)) + (-delta_H_NH3ox) * (-r_NH3))  / (F_NH3 * Cp_NH3 + F_O2 * Cp_O2 + F_N2 * Cp_N2 + F_H2O * Cp_H2O);

dydV = [C_H2O; C_O2; dC_N2_dV; dC_NH3_dV; dT_dV];

end
