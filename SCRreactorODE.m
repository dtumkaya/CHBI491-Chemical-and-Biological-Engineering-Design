% Simulation of SCR with copper zeolite catalyst

function dydV = SCRreactorODE(V, y, R, V_react, U, A_surface, d_h)

% Extract variables
C_NH3 = y(1);
C_N2 = y(2);
C_H2 = y(3);
C_O2 = y(4);
C_H2O = y(5);
C_NO2 = y(6);
C_NO = y(7);
C_N2O = y(8);
T = y(9);

T_ref = 298.15; % in K
C_tot = 1/(8.21*10^(-5) * (T + 273.15)); % mol/m3

% R1: NH3 oxidation
% 4 NH3 + 3 O2 -> 2 N2 + 6 H2O
a_c = 4/d_h; % 1/m
Sh = 3;
D_ab_NH3 = 1.215 *(10^-9)*T^1.7389; % m2/s
k_m = Sh*D_ab_NH3/d_h; % m/s
t_w = 5*10^-5; % m  

D_eff_NH3 = 100*1.215 *(10^-9)*T^1.7389; % m2/s
A_NH3 = 1.79 * 10^12;
Ea_NH3 = 72700;

k_s_NH3 = A_NH3 * exp(-(Ea_NH3/(R*T)));
tiele_NH3 = t_w * sqrt(k_s_NH3/D_eff_NH3);
eff_n_NH3 = tanh(tiele_NH3)/tiele_NH3;

C_NH3_w = (k_m * a_c * C_NH3) / (eff_n_NH3 * k_s_NH3 + k_m * a_c); % (mol/m^3*s) / (1/s) -> mol/m^3
r_NH3 = - k_s_NH3 * C_NH3_w;

% R2: Standard SCR
% 4 NH3 + 4 NO + O2 -> 4 N2 + 6 H2O
a_c = 4/d_h; % 1/m
Sh = 3;
D_ab_NO = 1.2365*(10^-9)*T^1.7006; % m2/s
k_m = Sh*D_ab_NO/d_h; % m/s
t_w = 5*10^-5; % m  

D_eff_NO = 100*1.2365*(10^-9)*T^1.7006; % m2/s
A_NO = 8.32*10^5 * 258 * C_O2^0.5;

Ea_NO = 64600; % j/mol
k_s_NO = A_NO * exp(-(Ea_NO/(R*T)));

tiele_NO = t_w * sqrt(k_s_NO/D_eff_NO);
eff_n_NO = tanh(tiele_NO)/tiele_NO;

C_NO_w = (k_m * a_c * C_NO) / (eff_n_NO * k_s_NO + k_m * a_c); % (mol/m^3*s) / (1/s) -> mol/m^3
r_NO_st = - k_s_NO * C_NO_w * 0.7;

% R3: Fast SCR
% 2NO + 2NO2 + 4NH3 → 4N2 + 6H2O
a_c = 4/d_h; % 1/m
Sh = 3;
D_ab_NO2 = 7.9236 *(10^-10)*T^1.7297; % m2/s
k_m = Sh*D_ab_NO2/d_h; % m/s
t_w = 5*10^-5; % m 

D_eff_NO2 = 100*7.9236 *(10^-10)*T^1.7297; % m2/s
A_fast = 3.3*10^11 * 258 * 0.0082;  % C_NO as constant
Ea_fast = 82000; % j/mol

k_s_fast = A_fast * exp(-(Ea_fast/(R*T)));

tiele_fast = t_w * sqrt(k_s_fast/D_eff_NO2);
eff_n_fast = tanh(tiele_fast)/tiele_fast;

C_NO2_w = (k_m * a_c * C_NO2) / (eff_n_fast * k_s_fast + k_m * a_c); % (mol/m^3*s) / (1/s) -> mol/m^3
r_fast = - k_s_fast * C_NO2_w * 0.7;

% D_ab_NO = 1.2365 *(10^-9)*T^1.7297; % m2/s
% D_eff_NO = 100*1.2365 *(10^-9)*T^1.7297; % m2/s

% R4: N2O Formation
% 2NH3 + 2NO2 → N2 + N2O + 3H2O
a_c = 4/d_h; % 1/m
Sh = 3;
D_ab_NO2 = 7.9236 *(10^-10)*T^1.7297; % m2/s

k_m = Sh*D_ab_NO2/d_h; % m/s
t_w = 5*10^-5; % m 
D_eff_NO2 = 100*7.9236 *(10^-10)*T^1.7297; % m2/s

A_N2O = 1.17*10^3*258;
Ea_N2O = 31800; % j/mol

k_s_N2O = A_N2O * exp(-(Ea_N2O/(R*T)));
tiele_N2O = t_w * sqrt(k_s_N2O/D_eff_NO2);
eff_n_N2O = tanh(tiele_N2O)/tiele_N2O;

C_N2O_w = (k_m * a_c * C_NO2) / (eff_n_N2O * k_s_N2O + k_m * a_c); % (mol/m^3*s) / (1/s) -> mol/m^3
r_N2O = - k_s_N2O * C_N2O_w * 0.7;

%%
% R5: NO oxidation
% delta_S = -76.1; % j/mol*K

A7 = 6.81*10^5; 
Ea_7 = 12390; % j/mol

% K_P = exp((-delta_H + delta_S*T) / (R*T)); % unitless
% K_C = K_P /(8.21*10^(-5)*T)^(0.5); % 1/(m^3/mol)^(-0.5) = (mol/m^3)^(-0.5)

K_eq = 8.61*10^-4*exp(57280/R*T);
k_7 = A7 * exp(-(Ea_7/(R*T))); % 1/s 
r_NO = - k_7 * (C_NO * C_O2^0.5 - C_NO2/K_eq); 

% Mass balance equations
% volumetric flow rate
MW_air = 28.96; % g/mol
mass_flow = 1300000/3600; % g/s
F_tot = (mass_flow/MW_air); % mol/s
vflow = F_tot/C_tot; % m^3/s

F_NH3 = C_NH3 * vflow; % mol/s
F_N2 = C_N2 * vflow;     
F_H2 = C_H2 * vflow;   
F_O2 = C_O2 * vflow;
F_H2O = C_H2O * vflow; % mol/s
F_NO2 = C_NO2 * vflow;     
F_NO = C_NO * vflow;   
F_N2O = C_N2O * vflow;

dC_NH3_dV = (r_NH3 + r_NO_st + 2*r_fast + r_N2O)/vflow; % mol/s
dC_N2_dV = -(0.5 *r_NH3 + r_NO_st+ 2*r_fast + 0.5*r_N2O) /vflow; % mol/s
dC_H2_dV = C_H2;
dC_O2_dV = C_O2;
dC_H2O_dV = C_H2O;
dC_NO2_dV = ((r_fast+ r_N2O) - r_NO)/vflow; % mol/s;
dC_NO_dV =  (r_NO_st + r_fast + r_NO)/vflow; % mol/s;
dC_N2O_dV = -(0.5*r_N2O)/vflow; % mol/s;

% Specific heat capacity of NO & O2 (J/(mol·K))
% https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Mask=1&Type=JANAFG&Table=on
Cp_NH3 = 19.99 + 49.77 *(T/1000) - 15.38*(T/1000)^2 + 1.92*(T/1000)^3 + 0.189/(T/1000)^2;
Cp_N2 = 18.98 + 1.85 *(T/1000) - 9.65*(T/1000)^2 + 16.64*(T/1000)^3 + 0.0001/(T/1000)^2;
Cp_H2 = 33.07 - 11.36 *(T/1000) + 11.4*(T/1000)^2 - 2.77*(T/1000)^3 - 0.16/(T/1000)^2;
Cp_O2 = 31.32 - 20.24*(T/1000) + 57.87*(T/1000)^2 - 36.5 *(T/1000)^3 - 0.007374/(T/1000)^2;
Cp_H2O = 30.09 + 6.83*(T/1000) + 6.79*(T/1000)^2 - 2.53 *(T/1000)^3 + 0.082/(T/1000)^2;
Cp_NO2 = 16.11 + 75.9 *(T/1000) - 54.4 *(T/1000)^2 + 14.3*(T/1000)^3 + 0.24/(T/1000)^2;
Cp_NO = 23.84 + 12.6*(T/1000) - 1.14*(T/1000)^2 - 1.5*(T/1000)^3 + 0.22/(T/1000)^2; 
Cp_N2O = 27.68 + 51.15*(T/1000) - 30.64*(T/1000)^2 - 6.85*(T/1000)^3 - 0.16/(T/1000)^2; 

% Assuming constant enthalpy change for simplicity
delta_NH3_ox = -312000;
delta_H_NO = -57000; 
delta_stan = -407000;
delta_N2O= -328000;
delta_fast = -378000;
dT_dV = ((U*A_surface*(T_ref-T)) + (-delta_N2O) * (-r_N2O) + (-delta_fast)* (-r_fast) + (-delta_stan) * (-r_NO_st) + (-delta_H_NO) * (-r_NO) + (-delta_NH3_ox) * (-r_NH3))  / (F_NH3 * Cp_NH3 + F_N2 * Cp_N2 + F_N2O * Cp_N2O + F_H2O*Cp_H2O + F_NO * Cp_NO + F_O2 * Cp_O2 + F_H2 * Cp_H2 + F_NO2 * Cp_NO2);

dydV = [dC_NH3_dV; dC_N2_dV; dC_H2_dV; dC_O2_dV; dC_H2O_dV; dC_NO2_dV; dC_NO_dV; dC_N2O_dV; dT_dV];

end
