%% Example 3: Compound Antenna with 4 Subarrays (4x4 Total Elements)
% Demonstrates compound array configuration with 4 subarrays, each 2x2,
% using multi-port feeding and reception modeling.

clear; clc; close all;

%% --- Constants & Parameters ---
mu0    = 4*pi*1e-7;
c      = physconst('LightSpeed');
f0     = 2.4e9;
lambda = c / f0;

% Element Parameters
L      = 0.02;            % Physical length of the capacitor element (m)
T      = 0.001;           % Thickness (m)
C      = 470e-12;         % Capacitance (F)
R      = 0.0028;               % Series resistance (Ohm)
Ind    = 1e-9;            % Lumped inductance (H)
V0     = 1;               % Feed voltage (V)
Z0     = 5;              % Source impedance (Ohm)
D      = 10;  

% Subarray Settings
nSub  = 4;
nR    = 2;        % 2x2 subarrays
nC    = 2;
nEl   = nR * nC;
subL  = 0.04;     % Subarray length
subW  = 0.04;
spacing = 0.06;   % Spacing between subarray centers

%% --- Define Subarray Layout ---
CompoundDefs = repmat(struct, nSub, 1);
offsets = [-1, -1; 1, -1; -1, 1; 1, 1] * spacing/2;

for i = 1:nSub
    CompoundDefs(i).NumRows = nR;
    CompoundDefs(i).NumCols = nC;
    CompoundDefs(i).Alength = subL;
    CompoundDefs(i).Awidth  = subW;
    CompoundDefs(i).ArrayOffset = [offsets(i,1), offsets(i,2), 0];
    CompoundDefs(i).Voltage = V0;  % Optional: vary for beamforming
    CompoundDefs(i).Label = ['Sub' num2str(i)];
end

N = nSub * nEl;

%% --- Create Compound Antenna ---
antenna = SingleLayerCapacitorAntenna_new(N, L, T, subL, subW, ...
    f0, V0, C, mu0, D, Ind, R, Z0, ...
    'ConfigurationType','compound', ...
    'CompoundArrayDefinitions', CompoundDefs, ...
    'MultiPort', true, ...
    'CenterCompound', true);

%% --- Visualize Geometry ---
antenna.show();

%% --- Analyze Currents ---
I = antenna.currentDistribution(f0);
fprintf('Element Currents (Abs / Phase deg):\n');
disp([abs(I), rad2deg(angle(I))]);

%% --- Pattern Visualization ---
antenna.visualize3DPattern();
antenna.patternAzimuth(0);
antenna.patternElevation(0);

%% --- Performance Evaluation ---
D = antenna.calculateDirectivity();
G = antenna.calculateGain();
eff = antenna.calculateRadiationEfficiency(f0);
fprintf('Directivity: %.2f dBi | Gain: %.2f dBi | Efficiency: %.2f%%\n', D, G, eff*100);

%% --- Imaging-Style Reception Simulation ---
antenna.OperatingMode = 'sensor';
azScan = -60:10:60;
elScan = -20:10:20;
data = antenna.performImagingScan(azScan, elScan, f0, 1);

% Find max response direction
Vmax = 0; bestIdx = 0;
for k = 1:numel(data)
    V = sum(abs(data(k).FeedVolt));  % total voltage across feeds
    if V > Vmax
        Vmax = V;
        bestIdx = k;
    end
end
fprintf('Max received voltage = %.3f V at Az = %d°, El = %d°\n', ...
    Vmax, data(bestIdx).Az, data(bestIdx).El);

%% --- Save Pattern for Post-Processing ---
antenna.exportPatternData('SLC_Compound_Pattern.mat', f0);

