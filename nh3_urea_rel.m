% Reactor sizing for urea decomposition unit
% With 0.02 s (residence time) -> L = 65 cm

% Universal gas constant in J/(mol*K)
R = 8.314;

% Temperature in Kelvin 
T = 298.15;

k1 = 4.9 * 10^3 * exp(-5505 / (R * T)); % 1/s

% We need to supply CA_0 / 0.3 AdBlue since it has %30 urea.
CA_0 =  1.15 * 1.0625*10^-2; % mol/m^3 -> urea -> 530 ppm
tau_range = linspace(0, 0.1, 100); % seconds

Q = 1.6 * 1.67 * 10^-5; % m^3/s (1.6 L/min - vol. flow rate)
A = pi*0.005^2; % m^2 (r - 5 cm ?)
L_values = tau_range * Q / A;

fprintf("\n %d ", 0.02 * Q);
fprintf("\n %d ", Q);
fprintf("\n %d ", A);
fprintf("\n %d ", 1.5 * A);

[CA_values, CB_values] = arrayfun(@(tau) concentrationA(tau, CA_0, k1), tau_range);

% Plotting the concentration results
figure;
plot(tau_range, CA_values);
xlabel('Residence Time (s)');
ylabel('Concentration of A (mol/m3)');
title('Concentration of urea in PFR');

figure;
plot(tau_range, CB_values);
xlabel('Residence Time (s)');
ylabel('Concentration of B (mol/m3)');
title('Concentration of NH3 in PFR');

% Plotting the reactor length results
figure;
plot(tau_range, L_values * 10^2);
xlabel('Residence Time (s)');
ylabel('Reactor Length (m)');
title('Reactor Length vs. Residence Time');

% Function to calculate the concentration of urea
function [CA, CB] = concentrationA(tau, C0_urea, k1)
    CA = C0_urea * exp(-k1 * tau);
    CB = C0_urea * (1 - exp(-k1 * tau));
end
