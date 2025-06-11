%==========================================================================
% rwr_deinterleave_full_dataset_v2.m
%
% Extended dataset generator:
%  • All combinations of 1–4 simultaneous radars
%  • 5 fixed emitter sites; sites #1–2 restricted to “older tech”
%  • 3 aircraft flight lines
%  • No randomness—fully reproducible
%
% Outputs:
%   dataset(n).rxData, .fs, .radarIdx, .radarTypes,
%   .txSigs, .txParams, .flightIdx, .flightPath
%==========================================================================

clearvars; close all; clc;

%% 1) Constants & common time axis
c        = physconst('LightSpeed');  % Speed of light [m/s]
fs_common= 50e6;                     % Common sampling rate [Hz]
t_max    = 1e-3;                     % Capture window [s]
t_vec    = (0:1/fs_common:t_max-1/fs_common).';  % time samples

%% 2) 5 Ground‐based emitter positions (x,y,z in meters)
emitters = [...
     0,     0,     0;    % Site #1 – old‐tech only
  5000,  1000,     0;    % Site #2 – old‐tech only
 -2000,  3000,     0;    % Site #3 – any tech
 10000, -2000,     0;    % Site #4 – any tech
   -500,  8000,    0]';  % Site #5 – any tech

%% 3) Radar waveform definitions (7 types)
radars = { ...
  struct('type','Rectangular','fs',20e6,   'PRF',1000, 'pulseWidth',30e-6), ...
  struct('type','LFM',         'fs',50e6,   'PRF',2000, 'pulseWidth',100e-6,'sweepBW',20e6), ...
  struct('type','PhaseCoded',  'fs',15e6,   'PRF',1500, 'chipWidth',2e-6, 'code',[1 1 1 -1 -1 1 -1]), ...
  struct('type','CW',          'fs',10e6,                  'fc',24.125e9,'duration',t_max), ...
  struct('type','FMCW',       'fs',100e6,                 'sweepTime',200e-6,'sweepBW',150e6,'fc',77e9), ...
  struct('type','PulseDoppler','fs',100e6,  'PRF',5000, 'pulseWidth',10e-6,'numPulses',64,'fc',3e9), ...
  struct('type','FSK',         'fs',200e6,                 'baseFreq',10e9,'freqSteps',[0 50e6 100e6],...
                                  'pulseWidth',20e-6,'PRI',300e-6) ...
};
numRadars = numel(radars);

%% 4) Define which sites are “older‐tech only” and what counts as old
oldSites       = [1,2];  % emitter indices restricted to old tech
oldTechTypes   = {'Rectangular','CW','PhaseCoded'};

%% 5) Aircraft flight‐paths (initPos & vel)
flights = { ...
  struct('initPos',[-3000; 1000; 1000],'vel',[200;   0;   0]), ...
  struct('initPos',[-4000;  500; 1500],'vel',[150;  50;    0]), ...
  struct('initPos',[-2500; 1500; 1200],'vel',[180; -30;    0])  ...
};
numFlights = numel(flights);

%% 6) Build all combos of 1–4 radars, then filter by old‐tech rule
allCombos = {};
for k = 1:4
    Ck = nchoosek(1:numRadars, k);       % [numCombos x k]
    for row = 1:size(Ck,1)
        allCombos{end+1} = Ck(row,:);    % store each combo as row vector
    end
end

validCombos = {};
for i = 1:numel(allCombos)
    combo = allCombos{i};
    ok = true;
    for p = 1:numel(combo)
        siteIdx = p;   % emitter position index (1–5)
        rType   = radars{combo(p)}.type;
        if ismember(siteIdx, oldSites) && ~ismember(rType, oldTechTypes)
            ok = false; break;
        end
    end
    if ok
        validCombos{end+1} = combo; %#ok<SAGROW>
    end
end

%% 7) Pre‐allocate dataset struct
N = numel(validCombos) * numFlights;
dataset = struct('rxData',cell(N,1), 'fs',[], 'radarIdx',[], ...
                 'radarTypes',[], 'txSigs',[], 'txParams',[], ...
                 'flightIdx',[], 'flightPath',[]);
idx = 1;

%% 8) Loop over each valid combo & each flight
for c = 1:numel(validCombos)
    combo   = validCombos{c};
    k       = numel(combo);
    % Pre‐assign the k emitter positions (sites 1..k)
    txPos_k = emitters(:,1:k);
    
    for f = 1:numFlights
        F        = flights{f};
        rwrPos_t = F.initPos + F.vel * t_vec.';  % [3×N] dynamic RWR path
        rxTotal  = zeros(size(t_vec));          % composite received signal
        txSigs   = cell(k,1);
        txParams = cell(k,1);
        types    = cell(1,k);
        
        for m = 1:k
            R     = radars{ combo(m) };
            txPos = txPos_k(:,m);
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
                pwS  = round(R.pulseWidth * R.fs);
                raw  = zeros(priS*numel(R.freqSteps),1);
                tt   = (0:pwS-1)'/R.fs;
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
            midPos  = F.initPos + F.vel*(t_max/2);
            rng     = norm(midPos - txPos);
            % Compute integer sample delay
            d_s = 2 * rng / c;
            d_samp = round(d_s * R.fs);
            
            % Delay signal (or skip if it's out of frame)
            if d_samp >= length(t_vec)
                delayed = zeros(size(t_vec));  % signal arrives after capture ends
            else
                sig = [zeros(d_samp,1); txSig];
                sig = sig(1:min(length(sig), length(t_vec)));  % truncate
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
            
            % Accumulate signal
            rxTotal = rxTotal + delayed;
            
            txSigs{m}   = txSig;
            R.txPos = txPos;  % Attach emitter position to this radar's struct
            txParams{m} = R;
        end
        
        % --- Store scenario in dataset ---
        dataset(idx).rxData     = rxTotal;
        dataset(idx).fs         = fs_common;
        dataset(idx).radarIdx   = combo;
        dataset(idx).radarTypes = types;
        dataset(idx).txSigs     = txSigs;
        dataset(idx).txParams   = txParams;
        dataset(idx).flightIdx  = f;
        dataset(idx).flightPath = F;
        idx = idx + 1;
    end
end

%% 9) Save
save('RSD_full_dataset_v2.mat','dataset');
fprintf('Generated %d scenarios; saved to RSD_full_dataset_v2.mat\n',idx-1);