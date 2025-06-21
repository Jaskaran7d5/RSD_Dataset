%==========================================================================
% Dual-Emitter RWR Simulation with Basic Waveform Implementation
% Features:
%   - Two ground emitters at [0,0,0] and [2000,1000,0]
%   - 8 aircraft flight paths
%   - 7 radar waveform types implemented without Phased Array Toolbox
%   - Saves dataset to .mat file
%==========================================================================

clearvars; close all; clc;

%% 1) Constants & common time axis
c = 3e8; % Speed of light [m/s]
fs_common = 50e6; % Common sampling rate [Hz]
t_max = 1e-3; % Capture window [s]
t_vec = (0:1/fs_common:t_max-1/fs_common).'; % time samples

%% 2) Two Ground-based emitter positions (x,y,z in meters)
emitters = [0, 2000;  % x [m]
            0, 1000;  % y [m]
            0, 0];    % z [m] (both ground level)

%% 3) Radar waveform definitions (7 types - basic MATLAB implementation)
radars = { ...
    % Rectangular Pulse
    struct('type','Rectangular','fs',20e6, 'PRF',1000, 'pulseWidth',30e-6), ...
    % LFM
    struct('type','LFM', 'fs',50e6, 'PRF',2000, 'pulseWidth',100e-6,'sweepBW',20e6), ...
    % PhaseCoded
    struct('type','PhaseCoded', 'fs',15e6, 'PRF',1500, 'chipWidth',2e-6, 'code',[1 1 1 -1 -1 1 -1]), ...
    % CW
    struct('type','CW', 'fs',10e6, 'fc',24.125e9,'duration',t_max), ...
    % FMCW
    struct('type','FMCW', 'fs',100e6, 'sweepTime',200e-6,'sweepBW',150e6,'fc',77e9), ...
    % PulseDoppler
    struct('type','PulseDoppler','fs',100e6, 'PRF',5000, 'pulseWidth',10e-6,'numPulses',64,'fc',3e9), ...
    % FSK
    struct('type','FSK', 'fs',200e6, 'baseFreq',10e9,'freqSteps',[0 50e6 100e6],...
        'pulseWidth',20e-6,'PRI',300e-6) ...
};
numRadars = numel(radars);

%% 4) Aircraft flight-paths (8 linear, level paths)
flights = { ...
    struct('initPos',[-3000; 1000; 1000],'vel',[200; 0; 0]), ...
    struct('initPos',[-4000; 500; 1500],'vel',[150; 50; 0]), ...
    struct('initPos',[-2500; 1500; 1200],'vel',[180; -30; 0]), ...
    struct('initPos',[2000; -2000; 1000],'vel',[-220; 60; 0]), ...
    struct('initPos',[5000; 0; 1500],'vel',[-180; 0; 0]), ...
    struct('initPos',[-1000; -3000; 1200],'vel',[0; 210; 0]), ...
    struct('initPos',[0; 4000; 1000],'vel',[100; -150; 0]), ...
    struct('initPos',[-4500; 2500; 1500],'vel',[200; -80; 0]) ...
};
numFlights = numel(flights);

%% 5) Radar combinations (all pairs of radar types)
validCombos = {};
for i = 1:numRadars
    for j = 1:numRadars
        validCombos{end+1} = [i, j]; % All combinations of two radars
    end
end

%% 6) Pre-allocate dataset struct
N = numel(validCombos) * numFlights;
dataset = struct('rxData_ground',cell(N,1), 'rxData_aircraft',[], 'fs',[], ...
    'radarIdx',[], 'radarTypes',[], 'txSigs',[], 'txParams',[], ...
    'aircraftTxSig',[], 'aircraftTxParams',[], ...
    'flightIdx',[], 'flightPath',[], 'emitterPositions',[]);
idx = 1;

%% 7) Loop over each radar combination & each flight
for c = 1:numel(validCombos)
    combo = validCombos{c}; % [radarIdx1, radarIdx2]
    
    for f = 1:numFlights
        F = flights{f};
        rxTotal_ground = zeros(size(t_vec)); % received at aircraft
        rxTotal_aircraft = zeros(size(t_vec)); % received at ground
        txSigs = cell(2,1);
        txParams = cell(2,1);
        types = cell(1,2);
        
        % Process both emitters
        for emitterIdx = 1:2
            R = radars{combo(emitterIdx)};
            txPos = emitters(:,emitterIdx);
            types{emitterIdx} = R.type;
            
            % Generate waveform using basic MATLAB functions
            switch R.type
                case 'Rectangular'
                    [txSig, t_radar] = generateRectangularPulse(R, t_max);
                case 'LFM'
                    [txSig, t_radar] = generateLFMPulse(R, t_max);
                case 'PhaseCoded'
                    [txSig, t_radar] = generatePhaseCodedPulse(R, t_max);
                case 'CW'
                    [txSig, t_radar] = generateCWWaveform(R, t_max);
                case 'FMCW'
                    [txSig, t_radar] = generateFMCWWaveform(R, t_max);
                case 'PulseDoppler'
                    [txSig, t_radar] = generatePulseDopplerWaveform(R, t_max);
                case 'FSK'
                    [txSig, t_radar] = generateFSKPulse(R, t_max);
                otherwise
                    error('Unknown radar type: %s', R.type);
            end
            
            % Calculate propagation delay (midpoint approximation)
            midPos = F.initPos + F.vel*(t_max/2);
            rng = norm(midPos - txPos);
            delay_sec = rng / c;  % One-way delay
            
            % Apply delay and resample to common rate
            delayedSig = applyDelayAndResample(txSig, t_radar, delay_sec, R.fs, fs_common, t_vec);
            
            % Accumulate signal
            rxTotal_ground = rxTotal_ground + delayedSig;
            
            % Store TX details
            txSigs{emitterIdx} = txSig;
            R.txPos = txPos;
            txParams{emitterIdx} = R;
        end
        
        % Aircraft transmitter (simple CW)
        aircraftTx = struct('fc',10e9, 'fs',fs_common, 'type','CW');
        aircraftTxSig = exp(1j*2*pi*aircraftTx.fc*t_vec);
        
        % Propagate aircraft signal to ground (midpoint approximation)
        midPos = F.initPos + F.vel*(t_max/2);
        rng_aircraft = norm(midPos - [0;0;0]); % To origin
        delay_aircraft = rng_aircraft / c;
        d_samp_aircraft = round(delay_aircraft * fs_common);
        
        if d_samp_aircraft >= length(t_vec)
            rxTotal_aircraft = zeros(size(t_vec));
        else
            rxTotal_aircraft = [zeros(d_samp_aircraft,1); aircraftTxSig(1:end-d_samp_aircraft)];
        end
        
        % Store scenario
        dataset(idx).rxData_ground = rxTotal_ground;
        dataset(idx).rxData_aircraft = rxTotal_aircraft;
        dataset(idx).fs = fs_common;
        dataset(idx).radarIdx = combo;
        dataset(idx).radarTypes = types;
        dataset(idx).txSigs = txSigs;
        dataset(idx).txParams = txParams;
        dataset(idx).aircraftTxSig = aircraftTxSig;
        dataset(idx).aircraftTxParams = aircraftTx;
        dataset(idx).flightIdx = f;
        dataset(idx).flightPath = F;
        dataset(idx).emitterPositions = emitters;
        idx = idx + 1;
    end
end

%% 8) Save dataset
save('RSD_dual_emitter_basic_waveforms.mat', 'dataset', '-v7.3');
fprintf('Generated %d scenarios with dual emitters; saved to RSD_dual_emitter_basic_waveforms.mat\n', N);

%% Waveform Generation Functions (Basic MATLAB Implementation)
function [txSig, t] = generateRectangularPulse(R, t_max)
    fs = R.fs;
    t = 0:1/fs:t_max-1/fs;
    pulse_width_samples = round(R.pulseWidth * fs);
    pri_samples = round(1/R.PRF * fs);
    
    % Generate pulse train
    txSig = zeros(size(t));
    num_pulses = floor(t_max * R.PRF);
    for i = 0:num_pulses-1
        start_idx = i * pri_samples + 1;
        end_idx = min(start_idx + pulse_width_samples - 1, length(t));
        txSig(start_idx:end_idx) = 1;
    end
    txSig = txSig(:);
end

function [txSig, t] = generateLFMPulse(R, t_max)
    fs = R.fs;
    t = 0:1/fs:t_max-1/fs;
    pulse_width = R.pulseWidth;
    pri = 1/R.PRF;
    
    % Chirp parameters
    f0 = 0; % Start frequency
    f1 = R.sweepBW; % End frequency
    k = (f1 - f0)/pulse_width; % Chirp rate
    
    % Generate pulse train
    txSig = zeros(size(t));
    num_pulses = floor(t_max / pri);
    for i = 0:num_pulses-1
        t_start = i * pri;
        t_pulse = t(t >= t_start & t < t_start + pulse_width) - t_start;
        chirp = exp(1j*pi*k*t_pulse.^2); % LFM signal
        idx = find(t >= t_start & t < t_start + pulse_width);
        txSig(idx(1:length(chirp))) = chirp;
    end
    txSig = txSig(:);
end

function [txSig, t] = generatePhaseCodedPulse(R, t_max)
    fs = R.fs;
    t = 0:1/fs:t_max-1/fs;
    chip_dur = R.chipWidth;
    code = R.code;
    pri = 1/R.PRF;
    
    % Generate pulse train
    txSig = zeros(size(t));
    num_pulses = floor(t_max / pri);
    for i = 0:num_pulses-1
        t_start = i * pri;
        for chip = 1:length(code)
            chip_start = t_start + (chip-1)*chip_dur;
            chip_end = chip_start + chip_dur;
            idx = t >= chip_start & t < chip_end;
            txSig(idx) = code(chip);
        end
    end
    txSig = txSig(:);
end

function [txSig, t] = generateCWWaveform(R, t_max)
    fs = R.fs;
    t = 0:1/fs:t_max-1/fs;
    txSig = exp(1j*2*pi*R.fc*t);
    txSig = txSig(:);
end

function [txSig, t] = generateFMCWWaveform(R, t_max)
    fs = R.fs;
    t = 0:1/fs:t_max-1/fs;
    sweep_time = R.sweepTime;
    num_sweeps = ceil(t_max / sweep_time);
    
    % Generate FMCW signal
    txSig = [];
    for i = 1:num_sweeps
        t_sweep = 0:1/fs:min(sweep_time, t_max - (i-1)*sweep_time) - 1/fs;
        chirp = exp(1j*2*pi*(R.fc*t_sweep + (R.sweepBW/(2*sweep_time))*t_sweep.^2));
        txSig = [txSig, chirp];
    end
    txSig = txSig(1:length(t)).';
end

function [txSig, t] = generatePulseDopplerWaveform(R, t_max)
    fs = R.fs;
    t = 0:1/fs:t_max-1/fs;
    pulse_width = R.pulseWidth;
    pri = 1/R.PRF;
    num_pulses = R.numPulses;
    
    % Generate pulse Doppler waveform
    txSig = zeros(size(t));
    for i = 0:min(num_pulses-1, floor(t_max/pri)-1)
        t_start = i * pri;
        t_pulse = t(t >= t_start & t < t_start + pulse_width) - t_start;
        pulse = exp(1j*2*pi*R.fc*t_pulse);
        idx = find(t >= t_start & t < t_start + pulse_width);
        txSig(idx(1:length(pulse))) = pulse;
    end
    txSig = txSig(:);
end

function [txSig, t] = generateFSKPulse(R, t_max)
    fs = R.fs;
    t = 0:1/fs:t_max-1/fs;
    pri = R.PRI;
    pulse_width = R.pulseWidth;
    freqs = R.baseFreq + R.freqSteps;
    
    % Generate FSK waveform
    txSig = zeros(size(t));
    num_pulses = min(length(freqs), floor(t_max/pri));
    for i = 1:num_pulses
        t_start = (i-1)*pri;
        t_pulse = t(t >= t_start & t < t_start + pulse_width) - t_start;
        pulse = exp(1j*2*pi*freqs(i)*t_pulse);
        idx = find(t >= t_start & t < t_start + pulse_width);
        txSig(idx(1:length(pulse))) = pulse;
    end
    txSig = txSig(:);
end

function delayedSig = applyDelayAndResample(txSig, t_radar, delay_sec, fs_radar, fs_common, t_common)
    % Apply time delay
    delay_samples_radar = round(delay_sec * fs_radar);
    if delay_samples_radar >= length(txSig)
        delayed = zeros(size(txSig));
    else
        delayed = [zeros(delay_samples_radar,1); txSig(1:end-delay_samples_radar)];
    end
    
    % Resample to common rate
    if fs_radar ~= fs_common
        [P,Q] = rat(fs_common/fs_radar);
        resampled = resample(delayed, P, Q);
    else
        resampled = delayed;
    end
    
    % Align with common time vector
    if length(resampled) < length(t_common)
        delayedSig = [resampled; zeros(length(t_common)-length(resampled),1)];
    else
        delayedSig = resampled(1:length(t_common));
    end
end
