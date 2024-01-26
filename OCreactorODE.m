% Simulation of oxidation reactor (NO and H2) with Platinum catalyst
% NO, H2, O2, H2O, N2 in
% N2, NO, H2, NO2, H2, O2, H2O out
% NO and NO2 data will be transferred to next.
% Material balance for H2 includes mass transfer.

function dydV = OCreactorODE(V, y, R, V_react, U, A_surface, d_h)

% Extract variables
C_NO = y(1);
C_NO2 = y(2);
C_H2 = y(3);
T = y(4);

% Arrhenius equation for reactions
% k = A * exp(-Ea / (R * T));

% Instead of C_A, we will put values from mass balance terms
a_c = 4/d_h; % 1/m
Sh = 3;
D_AB = 5.486 * 10^-9 * T^1.6816; % m2/s
k_m = Sh*D_AB/d_h; % m/s
t_w = 5*10^-5; % m  

% R1: H2 oxidation
% H2 + 1/2 O2 -> H2O
D_eff_H2 = 3.6608 * 10^-7 * T^0.664; % m2/s 
E_a_H2 = 44600; % j/mol

k_s_H2 = 2.794 * 10^11 * exp(-E_a_H2/(R*T)); % on Platinum, 1/s
tiele_H2 = t_w * sqrt(k_s_H2/D_eff_H2); % unitless
eff_n_H2 = tanh(tiele_H2)/tiele_H2; % unitless

C_H2_w = (k_m * a_c * C_H2) / (eff_n_H2 * k_s_H2 + k_m * a_c); % (mol/m^3*s) / (1/s) -> mol/m^3
r_H2 = -k_s_H2*C_H2_w;

% R2: NO oxidation
% NO + 1/2 O2 -> NO2
delta_H = -58300; % j/mol
delta_S = -76.1; % j/mol*K

% num_Pt = 3.8*10^-3; % mol/kg cat.
C_tot = 1/(8.21*10^(-5) * (T + 273.15)); % mol/m3
C_O2 = C_tot * 0.11; % mol/m3

A7 = 3.65 * 10^8 / T;
Ea_7 = 52000; % j/mol

K_P = exp((-delta_H + delta_S*T) / (R*T)); % unitless
K_C = K_P /(8.21*10^(-5)*T)^(0.5); % 1/(m^3/mol)^(-0.5) = (mol/m^3)^(-0.5)
  
k_7 = A7 * exp(-(Ea_7/(R*T))); % 1/s 
r_NO = -k_7 * (C_NO * C_O2^0.5 - C_NO2/K_C); 

% fprintf("K: %d %d %d \n", K_P, k_7, k_s_H2 );

% Mass balance equations
% volumetric flow rate
MW_air = 28.96; % g/mol
mass_flow = 1300000/3600; % g/s
F_tot = (mass_flow/MW_air); % mol/s
vflow = F_tot/C_tot; % m^3/s

% Rate expression
dC_H2_dV = r_H2 /vflow; % mol/s
dC_NO_dV = r_NO/vflow;
dC_NO2_dV = -r_NO/vflow;

F_NO = C_NO * vflow; % mol/s
F_O2 = C_O2 * vflow;     
F_H2 = C_H2 * vflow;   
F_NO2 = C_NO2 * vflow;

% Specific heat capacity of NO & O2 & H2 (J/(mol·K))
% https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Mask=1&Type=JANAFG&Table=on
Cp_NO = 23.84 + 12.6*(T/1000) - 1.14*(T/1000)^2 - 1.5*(T/1000)^3 + 0.22/(T/1000)^2; 
Cp_O2 = 31.32 - 20.24*(T/1000) + 57.87*(T/1000)^2 - 36.5 *(T/1000)^3 - 0.007374/(T/1000)^2;
Cp_H2 = 33.07 - 11.36 *(T/1000) + 11.4*(T/1000)^2 - 2.77*(T/1000)^3 - 0.16/(T/1000)^2;
Cp_NO2 = 16.11 + 75.9 *(T/1000) - 54.4 *(T/1000)^2 + 14.3*(T/1000)^3 + 0.24/(T/1000)^2;

% Assuming constant enthalpy change for simplicity
T_ref = 298.15; % in K

% No! T bağlantılı ekle (ya da 300 derece)
delta_H_NO = -58300;  % By hess law and corresponding formation enthalpies (NO ox)
delta_H_H2 = -241834; % joule/ mol -> 57.8 kcal/mol for H2 ox
dT_dV = ((U*A_surface*(T_ref-T)) + (-delta_H_H2) * (-r_H2) + (-delta_H_NO) * (-r_NO))  / (F_NO * Cp_NO + F_O2 * Cp_O2 + F_H2 * Cp_H2 + F_NO2 * Cp_NO2);

dydV = [dC_NO_dV; dC_NO2_dV; dC_H2_dV; dT_dV];

end
