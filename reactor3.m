function reactor3()
V_react = 3 * 10^-3;          % Volume of the reactor (m^3) -> 5 L
U = 10;       % Overall heat transfer coefficient (J/(s * m^2 * K))
cpsi = 400; % From paper: Single and multisite detailed kinetic models for the adsorption and desorption of NO2 over Cu based NH3-SCR catalyst
A_surface = (6.54*10^-4)/cpsi; % Heat transfer surface area (m^2)

d_h = 0.001; % m (1 mm)
R = 8.314; % J/K*mol
T0 = 376 + 273.15;       % Initial temperature (K)

C_tot = 1/(8.21*10^(-5) * T0); % mol/m3
C_H2O = C_tot * 0.2;
C_O2 = C_tot * 0.11;
C_N2 = C_tot * 0.69;
C_NH3 = 5.49*10^-9;
T = T0;         % Initial temperature

%% Function calls to solve ODE
    
Vspan = [0 3 * 10^-3]; % Simulation over volume (m^-3) -> 5 L

% Solve ODEs
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
[V, y] = ode45(@ASreactorODE, Vspan, [C_H2O C_O2 C_N2 C_NH3 T], options, ...
 R, V_react, U, A_surface, d_h);

% Extract results
reactor3.C_H2O = y(:, 1);
reactor3.C_O2 = y(:, 2);
reactor3.C_N2 = y(:, 3);
reactor3.C_NH3 = y(:, 4);
reactor3.T = y(:, 5);

% Plot results
figure;
subplot(2, 1, 1);
plot(V, reactor3.C_NH3, 'g', 'LineWidth', 2);
xlabel('Volume (m^3)');
ylabel('C_{NH3} (mol/m^3)');
title('Concentration of NH3 vs Volume');

subplot(2, 1, 2);
plot(V, reactor3.T, 'm', 'LineWidth', 2);
xlabel('Volume (m^3)');
ylabel('Temperature (K)');
title('Temperature vs Volume');

end