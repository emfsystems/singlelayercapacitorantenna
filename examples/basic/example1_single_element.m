%% Example 1: Single Element Antenna at 2.4 GHz
% This script demonstrates basic usage of the SingleLayerCapacitorAntenna_new class
% for a single capacitor-based antenna with microstrip feed line.

clear; clc; close all;

%% --- Constants & Parameters ---
mu0    = 4*pi*1e-7;
c      = physconst('LightSpeed');
f0     = 2.4e9;           % Operating frequency (Hz)
lambda = c / f0;

% Antenna Element Parameters
L      = 0.02;            % Physical length of the capacitor element (m)
T      = 0.001;           % Thickness (m)
C      = 470e-12;         % Capacitance (F)
R      = 0.0028;               % Series resistance (Ohm)
Ind    = 1e-9;            % Lumped inductance (H)
V0     = 1;               % Feed voltage (V)
Z0     = 14;              % Source impedance (Ohm)
D      = 10;             % Far-field observation distance (m)

% Array dimensions (arbitrary for single element)
Al     = 0.03;
Aw     = 0.03;

%% --- Create the Antenna Object ---
antenna = SingleLayerCapacitorAntenna_new(1, L, T, Al, Aw, ...
    f0, V0, C, mu0, D, Ind, R, Z0, ...
    'ConfigurationType','single');

%% --- Visualize Geometry ---
antenna.show();

%% --- Impedance Calculation ---
[Zin, ~] = antenna.lineBasedImpedance(f0, 1);
disp(['Input Impedance at 2.4 GHz: ', num2str(real(Zin), '%.2f'), ' + j', num2str(imag(Zin), '%.2f'), ' Ohms']);

%% --- VSWR & Return Loss ---
VSWR = antenna.calculateVSWR(f0);
RL   = antenna.calculateReturnLoss(f0);
fprintf('VSWR: %.2f | Return Loss: %.2f dB\n', VSWR, RL);

%% --- Current Distribution ---
I = antenna.currentDistribution(f0);
disp(['Current Magnitude: ', num2str(abs(I)), ' A']);
disp(['Current Phase: ', num2str(rad2deg(angle(I))), ' degrees']);

%% --- Radiation Pattern Visualization ---
antenna.visualizeRadiationPattern();
antenna.patternAzimuth(0);
antenna.patternElevation(0);
antenna.visualize3DPattern();

%% --- Efficiency & Directivity ---
D = antenna.calculateDirectivity();
G = antenna.calculateGain();
eff = antenna.calculateRadiationEfficiency(f0);
fprintf('Directivity: %.2f dBi | Gain: %.2f dBi | Efficiency: %.2f%%\n', D, G, eff*100);
antenna.estimateBandwidth([], 2);

[t, h] = antenna.impulseResponse();
figure;
plot(t*1e9, real(h)); % Convert t to nanoseconds
xlabel('Time (ns)');
ylabel('Impulse Response (Real)');
title('Time-Domain Reflection');
grid on;

%% --- S-Parameter Export ---
antenna.exportSParameters('SLC_SingleElement.s1p');

%% --- Simulate Reception ---
antenna.OperatingMode = 'receive';
incField = struct('Az', 0, 'El', 0, 'Amplitude', 1);  % Broadside incidence
Vrx = antenna.simulateReception(incField, f0);
disp(['Received voltage (loaded): ', num2str(abs(Vrx)), ' V']);

