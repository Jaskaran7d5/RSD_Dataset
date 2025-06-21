% FMCW_PhaseCoded_Visualization_with_Two_Ground_Emitters.m
% Simulates and visualizes:
% 1. Aircraft-mounted FMCW signal (transmitted)
% 2. Ground-based phase-coded radar signal 1 (transmitted)
% 3. Ground-based phase-coded radar signal 2 (transmitted)
% 4. Received signal at the aircraft RWR (from both coded pulse reflections)

clear; clc;

%% Constants
c = physconst('LightSpeed');
fs = 100e6;               % Common sampling frequency (Hz)
T = 100e-6;               % Duration for single pulse
t = (0:1/fs:T-1/fs).';    % Time vector

%% Geometry
% Aircraft position (transmitter and receiver)
aircraftPos = [0; 0; 1000];  % 1 km altitude

% Ground radar positions
groundRadarPos1 = [1000; 0; 0];     % 1 km away, level ground
groundRadarPos2 = [1200; 500; 0];   % 1.2 km away, 0.5 km lateral offset

%% 1. Aircraft FMCW signal generation
sweepTime = T;  % Match total duration to sweep time
fmcwWav = phased.FMCWWaveform( ...
    'SweepTime', sweepTime, ...
    'SweepBandwidth', 200e6, ...
    'SampleRate', fs, ...
    'SweepDirection','Up');

fmcwSig = 0.2 * fmcwWav();  % Weakened aircraft signal amplitude

%% 2. Ground-based Phase-Coded radar signal 1
phaseCodeWav1 = phased.PhaseCodedWaveform( ...
    'Code','Barker', ...
    'ChipWidth',1e-6, ...
    'SampleRate',fs);

barkerPulse1 = phaseCodeWav1();
numRepeats1 = ceil(length(t) / length(barkerPulse1));
barkerSig1 = repmat(barkerPulse1, numRepeats1, 1);
barkerSig1 = barkerSig1(1:length(t));

%% 3. Ground-based Phase-Coded radar signal 2 (new emitter)
phaseCodeWav2 = phased.PhaseCodedWaveform( ...
    'Code','Barker', ...
    'ChipWidth',1e-6, ...
    'SampleRate',fs);

barkerPulse2 = phaseCodeWav2();
numRepeats2 = ceil(length(t) / length(barkerPulse2));
barkerSig2 = repmat(barkerPulse2, numRepeats2, 1);
barkerSig2 = barkerSig2(1:length(t));

%% 4. Reflection of Coded Signals at aircraft

% Free-space path loss system object
fc_barker = 1e9;  % 1 GHz carrier frequency for path loss
fspl = phased.FreeSpace( ...
    'SampleRate', fs, ...
    'OperatingFrequency', fc_barker, ...
    'TwoWayPropagation', true);

% For emitter 1
roundTripDelay1 = 2 * norm(aircraftPos - groundRadarPos1) / c;
delaySamples1 = round(roundTripDelay1 * fs);

attenuatedBarker1 = fspl(barkerSig1, groundRadarPos1, aircraftPos, [0;0;0], [0;0;0]);
delayedBarker1 = [zeros(delaySamples1,1); attenuatedBarker1];
delayedBarker1 = delayedBarker1(1:length(t));

% For emitter 2
roundTripDelay2 = 2 * norm(aircraftPos - groundRadarPos2) / c;
delaySamples2 = round(roundTripDelay2 * fs);

attenuatedBarker2 = fspl(barkerSig2, groundRadarPos2, aircraftPos, [0;0;0], [0;0;0]);
delayedBarker2 = [zeros(delaySamples2,1); attenuatedBarker2];
delayedBarker2 = delayedBarker2(1:length(t));

% Final received signal (sum of all)
rxSig = fmcwSig + delayedBarker1 + delayedBarker2;

%% Visualization
figure;

subplot(5,1,1);
plot(t*1e6, real(fmcwSig), 'LineWidth', 1.5);
title('Aircraft FMCW Transmitted Signal');
xlabel('Time (microseconds)'); ylabel('Amplitude (Normalized)');
grid on;
ylim([min(real(fmcwSig)) max(real(fmcwSig))]*1.1);

subplot(5,1,2);
plot(t*1e6, real(barkerSig1), 'LineWidth', 1.5);
title('Ground Phase-Coded Radar Transmitted Signal 1');
xlabel('Time (microseconds)'); ylabel('Amplitude (Normalized)');
grid on;
ylim([-1.5 1.5]);

subplot(5,1,3);
plot(t*1e6, real(barkerSig2), 'LineWidth', 1.5);
title('Ground Phase-Coded Radar Transmitted Signal 2');
xlabel('Time (microseconds)'); ylabel('Amplitude (Normalized)');
grid on;
ylim([-1.5 1.5]);

subplot(5,1,4);
plot(t*1e6, real(rxSig), 'LineWidth', 1.5);
title('Received Signal at Aircraft (FMCW + Both Reflected Coded Pulses)');
xlabel('Time (microseconds)'); ylabel('Amplitude (Normalized)');
grid on;
ylim([min(real(rxSig)) max(real(rxSig))]*1.1);

subplot(5,1,5);
spectrogram(rxSig, 256, 200, 256, fs, 'yaxis');
title('Spectrogram of Received Signal');
colormap('jet');

