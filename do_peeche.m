%% Radar PDW Simulation: Two Aircraft Pursuing One Target from Different Angles
clc; clear; close all;

% Constants
c = 3e8;            % Speed of light (m/s)
Fs = 50e6;           % Sampling rate (Hz)
T = 1;               % Simulation duration (seconds)
dt = 1/Fs;           % Time step
t = 0:dt:T;          % Time vector
fc = 1e9;            % Carrier frequency (Hz)

% Initialize scenario structure
scenario = struct();
scenario.Emitters = struct('Position', {}, 'Velocity', {}, 'Type', {}, ...
                           'ScanType', {}, 'Size', {}, 'Frequency', {}, ...
                           'PRF', {}, 'Power', {}, 'PDWs', {}, 'PulseDetails', {});

% Define target aircraft
target_init = [0; 1000; 10000];  % Initial position (m)
target_vel = [200; 0; 0];        % Velocity vector (m/s)

% Define two pursuing aircraft with different approach angles
pursuers = {
    % Pursuer 1: 45° approach
    struct(...
        'angle', 45, ...          % Approach angle (degrees)
        'distance', 5000, ...      % Initial distance (m)
        'speed_factor', 0.8, ...   % Relative to target speed
        'offset', [-5000*cosd(45); -5000*sind(45); 0] ...
    ), 
    % Pursuer 2: -30° approach
    struct(...
        'angle', -30, ...          % Approach angle (degrees)
        'distance', 6000, ...      % Initial distance (m)
        'speed_factor', 0.9, ...   % Relative to target speed
        'offset', [-6000*cosd(-30); -6000*sind(-30); 0] ...
    )
};

% Calculate target position over time
target_pos = target_init + target_vel .* t;

% Simulate each pursuer
for pursuerIdx = 1:length(pursuers)
    fprintf('Simulating pursuer %d at %d°...\n', pursuerIdx, pursuers{pursuerIdx}.angle);
    
    % Calculate pursuer initial position and velocity
    pursuer_init = target_init + pursuers{pursuerIdx}.offset;
    pursuer_vel = pursuers{pursuerIdx}.speed_factor * target_vel;
    
    % Calculate pursuer position over time
    pursuer_pos = pursuer_init + pursuer_vel .* t;
    
    % Calculate distance and delay
    distances = sqrt(sum((target_pos - pursuer_pos).^2, 1));
    delays = distances / c;
    delay_samples = round(delays * Fs);
    
    % Generate LFM pulse
    pulse_width = 10e-6;
    bandwidth = 2e6;
    t_pulse = 0:1/Fs:pulse_width;
    pulse = chirp(t_pulse, 0, pulse_width, bandwidth);
    
    % Create transmitted signal (pulse train)
    tx_signal = zeros(size(t));
    pri = 100e-6;  % Pulse repetition interval
    pulse_starts = 0:pri:T;
    
    for start_time = pulse_starts
        start_idx = round(start_time * Fs) + 1;
        end_idx = min(start_idx + length(pulse) - 1, length(t));
        tx_signal(start_idx:end_idx) = pulse(1:(end_idx-start_idx+1));
    end
    
    % Simulate received signal with Doppler
    rx_signal = zeros(size(t));
    for i = 1:length(t)
        if delay_samples(i) > 0 && (i + delay_samples(i)) <= length(t)
            % Doppler calculation
            rel_vel = target_vel - pursuer_vel;
            radial_vec = (target_pos(:,i) - pursuer_pos(:,i)) / distances(i);
            radial_vel = sum(rel_vel .* radial_vec);
            doppler_shift = (2 * radial_vel / c) * fc;
            path_loss = (c/(4*pi*fc*distances(i)))^2;
            rx_signal(i + delay_samples(i)) = rx_signal(i + delay_samples(i)) + ...
                tx_signal(i) * path_loss * exp(1j*2*pi*doppler_shift*t(i));
        end
    end
    
    % Add noise
    rx_signal = awgn(rx_signal, 20);  % 20dB SNR
    
    % Extract PDWs
    [pdws, pulse_details] = extractPDWs(rx_signal, Fs, t);
    
    % Create emitter struct for storage
    emitter = struct();
    emitter.Position = pursuer_init;  % Initial position
    emitter.Velocity = pursuer_vel;   % Velocity vector
    emitter.Type = 'LFM';             % Signal type
    emitter.ScanType = 'Circular';    % Scan pattern
    emitter.Size = [1 1];             % Emitter size
    emitter.Frequency = fc;            % Carrier frequency
    emitter.PRF = 1/pri;              % Pulse Repetition Frequency
    emitter.Power = 1000;             % Transmission power (W)
    emitter.PDWs = pdws;              % Pulse Descriptor Words
    emitter.PulseDetails = pulse_details; % Raw pulse data
    
    % Append to scenario emitters (using dynamic array expansion)
    scenario.Emitters(end+1) = emitter;
end

% Save dataset
save('multi_pursuer_pdw_dataset.mat', 'scenario', '-v7.3');
disp('Simulation complete. Dataset saved.');

%% --- PDW Extraction Function --- 
function [pdwTable, pulseDetails] = extractPDWs(signal, Fs, t)
    % Detect pulses using envelope detection
    envelope = abs(hilbert(real(signal)));
    threshold = 0.5 * max(envelope);
    above_threshold = envelope > threshold;
    
    % Find pulse start/end indices
    pulse_starts = find(diff(above_threshold) > 0) + 1;
    pulse_ends = find(diff(above_threshold) < 0);
    
    % Ensure equal number of starts/ends
    num_pulses = min(length(pulse_starts), length(pulse_ends));
    pulse_starts = pulse_starts(1:num_pulses);
    pulse_ends = pulse_ends(1:num_pulses);
    
    % Initialize PDW arrays
    TOA = zeros(num_pulses, 1);
    PW = zeros(num_pulses, 1);
    PA = zeros(num_pulses, 1);
    AMP = zeros(num_pulses, 1);
    RF = zeros(num_pulses, 1);
    pulseDetails = cell(num_pulses, 1);
    
    for i = 1:num_pulses
        % Extract pulse parameters
        pulse_samples = pulse_starts(i):pulse_ends(i);
        pulse_data = signal(pulse_samples);
        
        % Time of Arrival (TOA)
        TOA(i) = t(pulse_starts(i));
        
        % Pulse Width (PW)
        PW(i) = t(pulse_ends(i)) - t(pulse_starts(i));
        
        % Pulse Amplitude (PA) and overall Amplitude (AMP)
        PA(i) = max(abs(pulse_data));
        AMP(i) = 20*log10(PA(i)); % Convert to dB
        
        % Carrier Frequency (RF) - using phase analysis
        phase_diff = diff(unwrap(angle(pulse_data)));
        RF(i) = mean(phase_diff) * Fs / (2*pi);
        
        % Store pulse details
        pulseDetails{i} = struct(...
            'Signal', pulse_data, ...
            'Time', t(pulse_samples), ...
            'Index', pulse_samples);
    end
    
    % Create PDW table
    pdwTable = table(TOA, PW, PA, AMP, RF);
end
