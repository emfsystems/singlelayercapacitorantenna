%% Example 2: 8-Element Rectangular Array at 2.4 GHz
% This script demonstrates a small antenna array composed of 8 elements
% arranged in a 2x4 rectangular grid, excited via a common feed.

clear; clc; close all;

%% --- Constants & Parameters ---
mu0    = 4*pi*1e-7;
c      = physconst('LightSpeed');
f0     = 2.4e9;              % Operating frequency
lambda = c / f0;

% Antenna Element Parameters
L      = 0.02;            % Physical length of the capacitor element (m)
T      = 0.001;           % Thickness (m)
C      = 470e-12;         % Capacitance (F)
R      = 0.0028;               % Series resistance (Ohm)
Ind    = 1e-9;            % Lumped inductance (H)
V0     = 1;               % Feed voltage (V)
Z0     = 5;              % Source impedance (Ohm)
D      = 10;             % Far-field observation distance (m)

% Array Layout
nRows  = 2;
nCols  = 4;
N      = nRows * nCols;
Al     = 0.06;               % Total array length (m)
Aw     = 0.04;               % Total array width (m)

%% --- Create Rectangular Array Antenna ---
antenna = SingleLayerCapacitorAntenna_new(N, L, T, Al, Aw, ...
    f0, V0, C, mu0, D, Ind, R, Z0, ...
    'ConfigurationType','rectangular', ...
    'NumRows', nRows, ...
    'NumCols', nCols);

%% --- Visualize Array Layout ---
antenna.show();

%% --- Impedance & Coupling ---
Zmat = zeros(N);
for i = 1:N
    for j = 1:N
        if i == j
            Zmat(i,j) = antenna.lineBasedImpedance(f0, i);
        else
            Zmat(i,j) = antenna.mutualImpedance(f0, i, j);
        end
    end
end
disp('Mutual Impedance Matrix (Ohms):');
disp(Zmat);

%% --- Solve Currents and Display ---
I = antenna.currentDistribution(f0);
disp('Element Currents (Magnitude & Phase):');
disp([abs(I), rad2deg(angle(I))]);

%% --- Visualizations ---
antenna.visualizeCurrentDistributionInteractive();
antenna.visualizeRadiationPattern();
antenna.patternAzimuth(0);
antenna.patternElevation(0);
antenna.visualize3DPattern();

%% --- Performance Metrics ---
D = antenna.calculateDirectivity();
G = antenna.calculateGain();
eff = antenna.calculateRadiationEfficiency(f0);
fprintf('Array Directivity: %.2f dBi | Gain: %.2f dBi | Efficiency: %.2f%%\n', D, G, eff*100);

%% --- Bandwidth Estimation ---
antenna.estimateBandwidth([], 2);  % Default freq range, VSWR ≤ 2

%% --- S-Parameter Export ---
antenna.exportSParameters('SLC_RectangularArray.s1p');
%% --- Simulate Reception ---
antenna.OperatingMode = 'receive';
incField = struct('Az', 0, 'El', 0, 'Amplitude', 1);  % Broadside incidence
Vrx = antenna.simulateReception(incField, f0);
disp(['Received voltage (loaded): ', num2str(abs(Vrx)), ' V']);


