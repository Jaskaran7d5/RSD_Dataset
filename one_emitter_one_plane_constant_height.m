%==========================================================================

% rwr_deinterleave_full_dataset_v2.m

% Now: One ground emitter at the origin, and the aircraft also emits
%      Eight linear, level flight paths

% Outputs:
% dataset(n).rxData_ground, .rxData_aircraft, .fs, .radarIdx, .radarTypes,
% .txSigs, .txParams, .aircraftTxSig, .aircraftTxParams,
% .flightIdx, .flightPath

%==========================================================================

clearvars; close all; clc;

%% 1) Constants & common time axis
c = physconst('LightSpeed'); % Speed of light [m/s]
fs_common = 50e6; % Common sampling rate [Hz]
t_max = 1e-3; % Capture window [s]
t_vec = (0:1/fs_common:t_max-1/fs_common).'; % time samples

%% 2) 1 Ground‐based emitter position (x,y,z in meters)
emitters = [0; 0; 0];   % Only one emitter at the origin

%% 3) Radar waveform definitions (7 types)
radars = { ...
    struct('type','Rectangular','fs',20e6, 'PRF',1000, 'pulseWidth',30e-6), ...
    struct('type','LFM', 'fs',50e6, 'PRF',2000, 'pulseWidth',100e-6,'sweepBW',20e6), ...
    struct('type','PhaseCoded', 'fs',15e6, 'PRF',1500, 'chipWidth',2e-6, 'code',[1 1 1 -1 -1 1 -1]), ...
    struct('type','CW', 'fs',10e6, 'fc',24.125e9,'duration',t_max), ...
    struct('type','FMCW', 'fs',100e6, 'sweepTime',200e-6,'sweepBW',150e6,'fc',77e9), ...
    struct('type','PulseDoppler','fs',100e6, 'PRF',5000, 'pulseWidth',10e-6,'numPulses',64,'fc',3e9), ...
    struct('type','FSK', 'fs',200e6, 'baseFreq',10e9,'freqSteps',[0 50e6 100e6],...
        'pulseWidth',20e-6,'PRI',300e-6) ...
};
numRadars = numel(radars);

%% 4) Aircraft flight‐paths (8 linear, level paths)
flights = { ...
    struct('initPos',[-3000; 1000; 1000],'vel',[200; 0; 0]), ...      % original
    struct('initPos',[-4000; 500; 1500],'vel',[150; 50; 0]), ...     % original
    struct('initPos',[-2500; 1500; 1200],'vel',[180; -30; 0]), ...   % original
    struct('initPos',[2000; -2000; 1000],'vel',[-220; 60; 0]), ...   % new: west-north
    struct('initPos',[5000; 0; 1500],'vel',[-180; 0; 0]), ...        % new: straight west
    struct('initPos',[-1000; -3000; 1200],'vel',[0; 210; 0]), ...    % new: straight north
    struct('initPos',[0; 4000; 1000],'vel',[100; -150; 0]), ...      % new: southeast
    struct('initPos',[-4500; 2500; 1500],'vel',[200; -80; 0]) ...    % new: east-south
};
numFlights = numel(flights);

%% 5) Only single-radar scenarios (since only one ground emitter)
validCombos = {};
for i = 1:numRadars
    validCombos{end+1} = i; % Each radar type is its own scenario
end

%% 6) Pre‐allocate dataset struct
N = numel(validCombos) * numFlights;
dataset = struct('rxData_ground',cell(N,1), 'rxData_aircraft',[], 'fs',[], 'radarIdx',[], ...
    'radarTypes',[], 'txSigs',[], 'txParams',[], ...
    'aircraftTxSig',[], 'aircraftTxParams',[], ...
    'flightIdx',[], 'flightPath',[]);
idx = 1;

%% 7) Loop over each radar & each flight
for c = 1:numel(validCombos)
    combo = validCombos{c};
    k = 1; % Only one radar per scenario
    txPos_k = emitters; % Only one emitter position

    for f = 1:numFlights
        F = flights{f};
        rwrPos_t = F.initPos + F.vel * t_vec.'; % [3×N] dynamic RWR path
        rxTotal_ground = zeros(size(t_vec)); % received at aircraft from ground
        rxTotal_aircraft = zeros(size(t_vec)); % received at ground from aircraft
        txSigs = cell(k,1);
        txParams = cell(k,1);
        types = cell(1,k);

        % --- Ground emitter waveform ---
        for m = 1:k
            R = radars{ combo(m) };
            txPos = txPos_k;
            types{m} = R.type;

            % --- Generate TX waveform at R.fs samples/s ---
            switch R.type
                case 'Rectangular'
                    wf = phased.RectangularWaveform('SampleRate',R.fs, ...
                        'PRF',R.PRF,'PulseWidth',R.pulseWidth,'OutputFormat','Samples',...
                        'NumSamples', length(t_vec));
                    txSig = wf();
                case 'LFM'
                    wf = phased.LinearFMWaveform('SampleRate',R.fs, ...
                        'PRF',R.PRF,'PulseWidth',R.pulseWidth, ...
                        'SweepBandwidth',R.sweepBW,'OutputFormat','Samples',...
                        'NumSamples', length(t_vec));
                    txSig = wf();
                case 'PhaseCoded'
                    wf = phased.PhaseCodedWaveform('SampleRate',R.fs, ...
                        'PRF',R.PRF,'ChipWidth',R.chipWidth,'Code','Barker', ...
                        'OutputFormat','Samples', ...
                        'NumSamples', length(t_vec));
                    txSig = wf();
                case 'CW'
                    txSig = exp(1j*2*pi*R.fc*t_vec);
                case 'FMCW'
                    wf = phased.FMCWWaveform('SampleRate',R.fs, ...
                        'SweepTime',R.sweepTime,'SweepBandwidth',R.sweepBW);
                    txSig = wf();
                case 'PulseDoppler'
                    wf = phased.RectangularWaveform('SampleRate',R.fs, ...
                        'PRF',R.PRF,'PulseWidth',R.pulseWidth, ...
                        'NumPulses',R.numPulses,'OutputFormat','Pulses');
                    temp = wf();
                    txSig = interp1(linspace(0,1,length(temp)),temp,linspace(0,1,length(t_vec))).';
                case 'FSK'
                    priS = round(R.PRI * R.fs);
                    pwS = round(R.pulseWidth * R.fs);
                    raw = zeros(priS*numel(R.freqSteps),1);
                    tt = (0:pwS-1)'/R.fs;
                    for n=1:numel(R.freqSteps)
                        f0 = R.baseFreq + R.freqSteps(n);
                        idx0 = (n-1)*priS + (1:pwS);
                        raw(idx0) = exp(1j*2*pi*f0*tt);
                    end
                    txSig = interp1(linspace(0,1,length(raw)),raw,linspace(0,1,length(t_vec))).';
                otherwise
                    error('Unknown type %s',R.type);
            end

            % --- Propagate: static delay based on midpoint of flight ---
            midPos = F.initPos + F.vel*(t_max/2);
            rng = norm(midPos - txPos);

            % Compute integer sample delay
            d_s = 2 * rng / c;
            d_samp = round(d_s * R.fs);

            % Delay signal (or skip if it's out of frame)
            if d_samp >= length(t_vec)
                delayed = zeros(size(t_vec)); % signal arrives after capture ends
            else
                sig = [zeros(d_samp,1); txSig];
                sig = sig(1:min(length(sig), length(t_vec))); % truncate
                % Resample to common sample rate
                resampled = resample(sig, fs_common, R.fs);
                % Pad or truncate to match t_vec
                if length(resampled) < length(t_vec)
                    resampled(end+1:length(t_vec)) = 0;
                else
                    resampled = resampled(1:length(t_vec));
                end
                delayed = resampled;
            end

            % Accumulate signal at aircraft
            rxTotal_ground = rxTotal_ground + delayed;
            txSigs{m} = txSig;
            R.txPos = txPos; % Attach emitter position to this radar's struct
            txParams{m} = R;
        end

        % --- Aircraft as emitter: simple CW example, can be any waveform ---
        aircraftTx.fc = 10e9; % Example: 10 GHz CW
        aircraftTx.fs = fs_common;
        aircraftTx.type = 'CW';
        aircraftTxSig = exp(1j*2*pi*aircraftTx.fc*t_vec);

        % Propagate from moving aircraft to ground receiver at [0;0;0]
        % Use midpoint for static delay (as above)
        midPos_aircraft = F.initPos + F.vel*(t_max/2);
        rng_aircraft = norm(midPos_aircraft - emitters);

        d_s_aircraft = 2 * rng_aircraft / c;
        d_samp_aircraft = round(d_s_aircraft * fs_common);

        if d_samp_aircraft >= length(t_vec)
            delayed_aircraft = zeros(size(t_vec));
        else
            sig_air = [zeros(d_samp_aircraft,1); aircraftTxSig];
            sig_air = sig_air(1:min(length(sig_air), length(t_vec)));
            delayed_aircraft = sig_air;
        end

        rxTotal_aircraft = rxTotal_aircraft + delayed_aircraft;

        % --- Store scenario in dataset ---
        dataset(idx).rxData_ground = rxTotal_ground;      % What aircraft receives from ground
        dataset(idx).rxData_aircraft = rxTotal_aircraft;  % What ground receives from aircraft
        dataset(idx).fs = fs_common;
        dataset(idx).radarIdx = combo;
        dataset(idx).radarTypes = types;
        dataset(idx).txSigs = txSigs;
        dataset(idx).txParams = txParams;
        dataset(idx).aircraftTxSig = aircraftTxSig;
        dataset(idx).aircraftTxParams = aircraftTx;
        dataset(idx).flightIdx = f;
        dataset(idx).flightPath = F;
        idx = idx + 1;
    end
end

%% 8) Save
save('RSD_full_dataset_v2_withAircraftTx.mat','dataset');
fprintf('Generated %d scenarios with aircraft emitter; saved to RSD_full_dataset_v2_withAircraftTx.mat\n',idx-1);
