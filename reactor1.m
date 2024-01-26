function reactor1()
V_react = 3 * 10^-3;          % Volume of the reactor (m^3) -> 5 L
U = 10;       % Overall heat transfer coefficient (J/(s * m^2 * K))
cpsi = 400; % From paper: Single and multisite detailed kinetic models for the adsorption and desorption of NO2 over Cu based NH3-SCR catalyst
A_surface = (6.54*10^-4)/cpsi; % Heat transfer surface area (m^2)

d_h = 0.001; % hydrolic diameter - m (1 mm)
R = 8.314; % J/K*mol
T0 = 300 + 273.15;  % Initial temperature (K)

C_tot = 1/(8.21*10^(-5) * T0); % mol/m3
C_NO = C_tot * 5 * 10^-4;  % Initial concentration of NO (mol/m^3 - 500 ppm)
C_H2 = C_tot * 10^-3; % Initial concentration of H2 (mol/m^3 - 1000 ppm)
C_NO2 = 0;      % Initial concentration of NO2 (mol/m^3)
T = T0;         % Initial temperature
    
%% FEEDBACK FROM CLASS
% feed flow rate 1300 kg/h
% feed inlet = 300 derece
% feed molar fraction h2 = 1000 ppm (mole fraction * 1e6)
% feed molar fraction no = 500 ppm
% engine out %11 o2, 20% h2o, balance with N2 (rest)
% conversion higher that 95%

Vspan = [0 3 * 10^-3]; % Simulation over volume (m^-3) -> 5 L

% Solve ODEs
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
[V, y] = ode45(@OCreactorODE, Vspan, [C_NO C_NO2 C_H2 T], options, ...
 R, V_react, U, A_surface, d_h);

% Extract results
reactor1.C_NO = y(:, 1);
reactor1.C_NO2 = y(:, 2);
reactor1.C_H2 = y(:, 3);
reactor1.T = y(:, 4);

% Plot results
figure;
subplot(4, 1, 1);
plot(V, reactor1.C_NO, 'b', 'LineWidth', 2);
xlabel('Volume (m^3)');
ylabel('C_{NO} (mol/m^3)');
title('Concentration of NO vs Volume');

subplot(4, 1, 2);
plot(V, reactor1.C_NO2, 'g', 'LineWidth', 2);
xlabel('Volume (m^3)');
ylabel('C_{NO2} (mol/m^3)');
title('Concentration of NO2 vs Volume');

subplot(4, 1, 3);
plot(V, reactor1.C_H2, 'g', 'LineWidth', 2);
xlabel('Volume (m^3)');
ylabel('C_{H2} (mol/m^3)');
title('Concentration of H2 vs Volume');

subplot(4, 1, 4);
plot(V, reactor1.T, 'm', 'LineWidth', 2);
xlabel('Volume (m^3)');
ylabel('Temperature (K)');
title('Temperature vs Volume');

end