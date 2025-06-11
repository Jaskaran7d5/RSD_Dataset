## Radar Signal Deinterleaving Dataset Generation

This repository contains MATLAB code for generating a comprehensive, labeled dataset intended for Radar Signal Deinterleaving (RSD) tasks. The dataset simulates the reception of radar pulse signals at a Radar Warning Receiver (RWR) onboard an aircraft flying along multiple predefined trajectories. The received signals are generated as combinations of multiple ground-based radar emitters transmitting different pulse waveforms.

---

## Scenario Overview

The dataset simulates a scenario with:
- **5 ground-based radar emitter sites**, with known fixed coordinates in 3D space.
- **7 radar waveform types** representing modern and legacy radar systems.
- **3 different aircraft trajectories**, each defined by a starting position and constant velocity vector.

During each simulation, a random combination of 1 to 4 radars transmits simultaneously. The aircraft's RWR receives the superimposed signals, simulating realistic electromagnetic interference between multiple emitters.

---

## Radar Types Considered

The dataset includes 7 radar waveform types:
1. **Rectangular Pulse**
2. **Linear Frequency Modulated (LFM)**
3. **Phase-Coded Pulse (7-chip Barker)**
4. **Continuous Wave (CW)**
5. **Frequency Modulated Continuous Wave (FMCW)**
6. **Pulse-Doppler**
7. **Frequency Shift Keying (FSK)**

Radar sites 1 and 2 are constrained to use legacy waveforms (Rectangular, CW, Phase-Coded) to reflect realistic technology deployment.

---

## Aircraft Trajectories

Three flight paths are defined, each with a unique starting location and constant velocity. The RWR continuously samples received signals over a fixed time window (1 ms) at a uniform sampling rate (50 MHz).

---

## Dataset Contents

Each sample in the generated dataset includes:
- The received baseband signal (`rxData`) as a 1 ms complex time series.
- Metadata on radar types active in the scenario.
- Transmit waveform parameters and identities.
- Aircraft trajectory information.

This dataset can be used for training and evaluating supervised learning models for tasks such as:
- Radar waveform classification
- Emitter localization
- Pulse deinterleaving and clustering

---

## Code Flexibility

The dataset generation code is modular and configurable. Users can modify:
- **Number of emitters**
- **Radar waveform parameters** (e.g., PRF, pulse width, bandwidth)
- **Aircraft flight paths**
- **Sampling rate and duration**
- **Waveform combinations (e.g., max number of emitters per scenario)**

The code enforces key physical constraints such as integer-valued sample delays, valid waveform configurations (e.g., SampleRate/PRF must be an integer), and proper signal alignment during propagation simulation.

---

## Requirements

- MATLAB R2021a or later
- Phased Array System Toolbox
- Signal Processing Toolbox (for resampling and filtering)

---

## Usage

To generate the dataset, run the primary script (e.g., `RSD_DG.m`). The output is a MATLAB struct array (`dataset`) saved as a `.mat` file, containing the signals and all associated metadata.

For analysis, training, or export to other formats, the `dataset` can be processed further using MATLAB or Python (after conversion).

---

## License

This code is intended for academic and research use. Please cite appropriately when used in publications or shared publicly.

