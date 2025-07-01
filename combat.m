clc; clear; close all;

% === Constants ===
c = 3e8;                % Speed of light (m/s)
Fs = 20e6;              % Sampling rate (Hz)
T = 0.1;                % Duration (s)
dt = 1/Fs;              % Time step
t = 0:dt:T;             % Time vector
fc = 1e9;               % Carrier frequency (Hz)

% === Target ===
target_init = [0; 1000; 10000];
target_vel = [200; 0; 0];
target_pos = target_init + target_vel .* t;

% === Pursuers (Combat Spread Formation) ===
pursuers = {
    % Leader (left of center, 500m behind, 500m left)
    struct('offset', [-500; -500; 0], 'speed_factor', 0.8), 
    % Wingman (right of center, 500m behind, 500m right)
    struct('offset', [-500; +500; 0], 'speed_factor', 0.8)
};

% === Pulse Parameters ===
pulse_width = 10e-6;
bandwidth = 2e6;
t_pulse = 0:1/Fs:pulse_width;
pulse = chirp(t_pulse, 0, pulse_width, bandwidth);

pri = 100e-6; % Pulse Repetition Interval
pulse_starts = 0:pri:T;
tx_signal = zeros(size(t));

% === Generate Pulse Train ===
for k = 1:length(pulse_starts)
    idx_start = round(pulse_starts(k)*Fs) + 1;
    idx_end = min(idx_start + length(pulse) - 1, length(t));
    tx_signal(idx_start:idx_end) = pulse(1:(idx_end - idx_start + 1));
end

% === Preallocate Scenario ===
num_pursuers = length(pursuers);
Emitters(num_pursuers) = struct( ...
    'Position', [], ...
    'Velocity', [], ...
    'Type', '', ...
    'ScanType', '', ...
    'Size', [], ...
    'Frequency', [], ...
    'PRF', [], ...
    'Power', [], ...
    'PDWs', [], ...
    'PulseDetails', [] );


% === Parallel Simulation for Each Pursuer ===
parfor pursuerIdx = 1:num_pursuers
    fprintf('Simulating pursuer %d...\n', pursuerIdx);
    
    purs = pursuers{pursuerIdx};
    purs_init = target_init + purs.offset;
    purs_vel = purs.speed_factor * target_vel;
    purs_pos = purs_init + purs_vel .* t;

    % Distance, delay, doppler
    d = sqrt(sum((target_pos - purs_pos).^2, 1));
    delays = d / c;
    delay_samples = round(delays * Fs);
    rel_vel = target_vel - purs_vel;

    rx_signal = zeros(size(t));
    for i = 1:length(t)
        delay_idx = i + delay_samples(i);
        if delay_idx <= length(t)
            r_vec = (target_pos(:,i) - purs_pos(:,i)) / d(i);
            r_vel = sum(rel_vel .* r_vec);
            fd = (2 * r_vel / c) * fc;
            path_loss = (c / (4*pi*fc*d(i)))^2;

            rx_signal(delay_idx) = rx_signal(delay_idx) + ...
                tx_signal(i) * path_loss * exp(1j * 2 * pi * fd * t(i));
        end
    end

    rx_signal = awgn(rx_signal, 20);  % Add noise (SNR = 20 dB)

    [pdws, pulse_details] = extractPDWs(rx_signal, Fs, t);

    % Store Emitter
    emitter = struct();
    emitter.Position = purs_init;
    emitter.Velocity = purs_vel;
    emitter.Type = 'LFM';
    emitter.ScanType = 'Circular';
    emitter.Size = [1 1];
    emitter.Frequency = fc;
    emitter.PRF = 1/pri;
    emitter.Power = 1000;
    emitter.PDWs = pdws;
    emitter.PulseDetails = pulse_details;

    Emitters(pursuerIdx) = emitter;
end

% === Store Full Scenario ===
scenario = struct();
scenario.Emitters = Emitters;

save('cpu_multi_pursuer_dataset.mat', 'scenario', '-v7.3');
disp('CPU-based simulation complete. Dataset saved.');

% === PDW Extraction ===
function [pdwTable, pulseDetails] = extractPDWs(signal, Fs, t)
    envelope = abs(hilbert(real(signal)));
    threshold = 0.5 * max(envelope);
    above_threshold = envelope > threshold;

    pulse_starts = find(diff(above_threshold) > 0) + 1;
    pulse_ends = find(diff(above_threshold) < 0);

    num_pulses = min(length(pulse_starts), length(pulse_ends));
    pulse_starts = pulse_starts(1:num_pulses);
    pulse_ends = pulse_ends(1:num_pulses);

    TOA = zeros(num_pulses, 1);
    PW = zeros(num_pulses, 1);
    PA = zeros(num_pulses, 1);
    AMP = zeros(num_pulses, 1);
    RF = zeros(num_pulses, 1);
    pulseDetails = cell(num_pulses, 1);

    for i = 1:num_pulses
        idx = pulse_starts(i):pulse_ends(i);
        pulse_data = signal(idx);
        TOA(i) = t(pulse_starts(i));
        PW(i) = t(pulse_ends(i)) - TOA(i);
        PA(i) = max(abs(pulse_data));
        AMP(i) = 20*log10(PA(i));
        RF(i) = mean(diff(unwrap(angle(pulse_data)))) * Fs / (2*pi);

        pulseDetails{i} = struct('Signal', pulse_data, ...
                                 'Time', t(idx), ...
                                 'Index', idx);
    end

    pdwTable = table(TOA, PW, PA, AMP, RF);
end
