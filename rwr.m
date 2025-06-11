% Define the transmitted waveform parameters
fs = 4e9;                     % Sampling frequency for the systems (Hz)
fc = 1.8e9;                   % Operating frequency of the surveillance radar (Hz)
T = 3e-6;                     % Chirp duration (s)
PRF = 1/(15e-6);              % Pulse repetition frequency (Hz)
BW = 30e6;                    % Chirp bandwidth (Hz)
c= physconst('LightSpeed');   % Speed of light in air (m/s)

% Assume the surveillance radar is at the origin and is stationary
radarPos= [0;0;0];          % Radar position (m)
radarVel= [0;0;0];          % Radar speed (m/s)

% Assume aircraft is moving with constant velocity
rwrPos= [-3000;1000;1000];   % Aircraft position (m)
rwrVel= [200; 0; 0];        % Aircraft speed (m/s)

% Configure objects to model ground radar and aircraft's relative motion
rwrPose = phased.Platform(rwrPos, rwrVel);
radarPose = phased.Platform(radarPos, radarVel);

% Configure the LFM waveform using the waveform parameters defined above
wavGen= phased.LinearFMWaveform('SampleRate',fs,'PulseWidth',T,'SweepBandwidth',BW,'PRF',PRF);

% Configure the Uniform Rectangular Array
antennaTx = phased.URA('ElementSpacing',repmat((c/fc)/2, 1, 2), 'Size', [8,8]);

% Configure objects for transmitting and propagating the radar signal
tx = phased.Transmitter('Gain', 5, 'PeakPower',100);
radiator = phased.Radiator( 'Sensor', antennaTx, 'OperatingFrequency', fc);
envIn = phased.FreeSpace('TwoWayPropagation',false,'SampleRate', fs,'OperatingFrequency',fc);

% Transmit a train of pulses
numPulses = 4;
txPulseTrain = helperRWR('simulateTransmission',numPulses, wavGen, rwrPos, ...
    radarPos, rwrVel, radarVel, rwrPose, radarPose, tx, radiator, envIn,fs,fc,PRF);

% Observe the signal arriving at the RWR
pspectrum(txPulseTrain,fs,'spectrogram','FrequencyLimits',[1.7e9 1.9e9], 'Leakage',0.65)
title('Transmitted Pulse Train Spectrogram')
caxis([-110 -90])