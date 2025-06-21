% FMCW_PhaseCoded_Visualization.m
% Simulates and visualizes:
% 1. Aircraft-mounted FMCW signal (transmitted)
% 2. Ground-based phase-coded radar signal (transmitted)
% 3. Received signal at the aircraft RWR (from coded pulse reflection)

clear; clc;

%% Constants
c = physconst('LightSpeed');
fs = 100e6;               % Common sampling frequency (Hz)
T = 100e-6;               % Duration for single pulse
t = (0:1/fs:T-1/fs).';    % Time vector

%% Geometry
% Aircraft position (transmitter and receiver)
aircraftPos = [0; 0; 1000];  % 1 km altitude
% Ground radar position (center of 5x5x5m cube)
groundRadarPos = [1000; 0; 0];  % 1 km away, level ground

% Approximate range for FMCW signal (to ground cube center)
rngTrue = norm(aircraftPos - groundRadarPos);

%% 1. Aircraft FMCW signal generation
sweepTime = T;  % Match total duration to sweep time
fmcwWav = phased.FMCWWaveform( ...
    'SweepTime', sweepTime, ...
    'SweepBandwidth', 200e6, ...
    'SampleRate', fs, ...
    'SweepDirection','Up');

fmcwSig = 0.2 * fmcwWav();  % Weakened aircraft signal amplitude

%% 2. Ground-based Phase-Coded radar signal
phaseCodeWav = phased.PhaseCodedWaveform( ...
    'Code','Barker', ...
    'ChipWidth',1e-6, ...
    'SampleRate',fs);

barkerPulse = phaseCodeWav();
numRepeats = ceil(length(t) / length(barkerPulse));
barkerSig = repmat(barkerPulse, numRepeats, 1);
barkerSig = barkerSig(1:length(t));


%% 3. Reflection of Coded Signal at aircraft
roundTripDelay = 2 * norm(aircraftPos - groundRadarPos) / c;
delaySamples = round(roundTripDelay * fs);

fc_barker = 1e9;  % 1 GHz carrier frequency for path loss
fspl = phased.FreeSpace( ...
    'SampleRate', fs, ...
    'OperatingFrequency', fc_barker, ...
    'TwoWayPropagation', true);

attenuatedBarker = fspl(barkerSig, groundRadarPos, aircraftPos, [0;0;0], [0;0;0]);

delayedBarker = [zeros(delaySamples,1); attenuatedBarker];
delayedBarker = delayedBarker(1:length(t));

% Final received signal
rxSig = fmcwSig + delayedBarker;

%% Visualization
figure;
subplot(4,1,1);
plot(t*1e6, real(fmcwSig), 'LineWidth', 1.5);
title('Aircraft FMCW Transmitted Signal');
xlabel('Time (microseconds)'); ylabel('Amplitude (Normalized)');
grid on;
ylim([min(real(fmcwSig)) max(real(fmcwSig))]*1.1);

subplot(4,1,2);
plot(t*1e6, real(barkerSig), 'LineWidth', 1.5);
title('Ground Phase-Coded Radar Transmitted Signal');
xlabel('Time (microseconds)'); ylabel('Amplitude (Normalized)');
grid on;
ylim([-1.5 1.5]);

subplot(4,1,3);
plot(t*1e6, real(rxSig), 'LineWidth', 1.5);
title('Received Signal at Aircraft (FMCW + Reflected Coded Pulse)');
xlabel('Time (microseconds)'); ylabel('Amplitude (Normalized)');
grid on;
ylim([min(real(rxSig)) max(real(rxSig))]*1.1);

subplot(4,1,4);
spectrogram(rxSig, 256, 200, 256, fs, 'yaxis');
title('Spectrogram of Received Signal');
colormap('jet');
