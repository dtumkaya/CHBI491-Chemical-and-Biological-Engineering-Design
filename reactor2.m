function reactor2()
V_react = 20 * 10^-3;          % Volume of the reactor (m^3) -> 20 L
U = 10;       % Overall heat transfer coefficient (J/(s * m^2 * K))
cpsi = 400; % From paper: Single and multisite detailed kinetic models for the adsorption and desorption of NO2 over Cu based NH3-SCR catalyst
A_surface = (6.54*10^-4)/cpsi; % Heat transfer surface area (m^2)

d_h = 0.001; % m (1 mm)
R = 8.314; % J/K*mol
T0 = 370 + 273.15;       % Initial temperature (K)

C_tot = 1/(8.21*10^(-5) * T0); % mol/m3
C_NH3 = 0.0122; % from ammonia decomposition
C_N2 = C_tot * 0.69; % given
C_H2 = C_tot * 10^-3; % given
C_O2 = C_tot * 0.11; % given
C_H2O = C_tot * 0.2; % given
C_NO = 0.0082; % from OCR 
C_NO2 = 0.0024; % from OCR
C_N2O = 0;
T = T0;

%%
    
Vspan = [0 20 * 10^-3]; % Simulation over volume (m^-3) -> 5 L

% Solve ODEs
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
[V, y] = ode45(@SCRreactorODE, Vspan, [C_NH3 C_N2 C_H2 C_O2 C_H2O C_NO2 C_NO C_N2O T], options, ...
 R, V_react, U, A_surface, d_h);

% Extract results
C_NH3 = y(:, 1);
C_N2 = y(:, 2);
C_H2 = y(:, 3);
C_O2 = y(:, 4);
C_H2O = y(:, 5);
C_NO2 = y(:, 6);
C_NO = y(:, 7);
C_N2O = y(:, 8);
T = y(:, 9);

% Plot results
figure;
subplot(6, 1, 1);
plot(V, C_NH3, 'b', 'LineWidth', 2);
xlabel('Volume (m^3)');
ylabel('C_{NH3} (mol/m^3)');
title('Concentration of NH3 vs Volume');

subplot(6, 1, 2);
plot(V, C_NO2, 'g', 'LineWidth', 2);
xlabel('Volume (m^3)');
ylabel('C_{NO2} (mol/m^3)');
title('Concentration of NO2 vs Volume');

subplot(6, 1, 3);
plot(V, C_NO, 'g', 'LineWidth', 2);
xlabel('Volume (m^3)');
ylabel('C_{NO} (mol/m^3)');
title('Concentration of NO vs Volume');

subplot(6, 1, 4);
plot(V, C_N2O, 'g', 'LineWidth', 2);
xlabel('Volume (m^3)');
ylabel('C_{N2O} (mol/m^3)');
title('Concentration of N2O vs Volume');

subplot(6, 1, 5);
plot(V, C_N2, 'g', 'LineWidth', 2);
xlabel('Volume (m^3)');
ylabel('C_{N2} (mol/m^3)');
title('Concentration of N2 vs Volume');

subplot(6, 1, 6);
plot(V, T, 'm', 'LineWidth', 2);
xlabel('Volume (m^3)');
ylabel('Temperature (K)');
title('Temperature vs Volume');

end