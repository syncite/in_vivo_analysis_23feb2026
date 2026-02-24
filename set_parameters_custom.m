function par = set_parameters_custom()
% SET_PARAMETERS_CUSTOM  wave_clus parameters for Plexon in vivo recordings.
%
% Optimised for:
%   - 40 kHz sampling rate
%   - Halorhodopsin/ChR2 optogenetics experiments
%   - Detecting low-amplitude L2/3 pyramidal cells (~0.1 Hz firing rates)
%   - Data pre-filtered at 250-5000 Hz with CAR (extract_and_filter.m)
%   - Long recordings (~45 min) on a 64 GB RAM machine
%
% Usage:
%   Place this file in your wave_clus directory (or on the MATLAB path)
%   and rename to set_parameters.m to override the defaults, OR pass it
%   explicitly:
%       par = set_parameters_custom();
%       Get_spikes({filename}, 'par', par);
%       Do_clustering(spkfile, 'par', par);
%
% To use with the wave_clus GUI:
%   Copy this file to your wave_clus folder and rename it set_parameters.m
%   (back up the original first). The GUI calls set_parameters() internally.

%% SEGMENT LENGTH
% Default is 5 min, which fragments sparse units across segment boundaries.
% At 0.1 Hz a 5-min segment has ~30 spikes — borderline for clustering.
% With 64 GB RAM, the full recording fits in one SPC pass. The distance
% matrix for ~20k spikes at 8 Hz average detection rate is ~3.5 GB.
% If you get memory errors, try 20 (minutes).
par.segments_length = Inf;

%% SAMPLING RATE
% Only used if the data file doesn't contain sr.
% Plexon default is 40 kHz.
par.sr = 40000;

%% DETECTION PARAMETERS
% stdmin is the critical parameter for catching low-amplitude units.
% Default (5) misses L2/3 pyramidals whose spikes sit at 3.5-4.5x noise.
% 4 is a good balance — more noise crossings reach the clustering stage
% but end up in cluster 0. Go to 3.5 if you're still missing units and
% can tolerate longer sorting sessions.
par.stdmin = 4;
par.stdmax = 50;

% Detection polarity. 'both' catches units regardless of which phase
% crosses threshold first. Costs ~2x detection time but important if
% your electrode orientation varies across channels.
par.detection = 'both';

% Refractory period. 1.5 ms is fine for cortical neurons (absolute
% refractory ~1 ms). Don't go lower or you'll get double-detections
% from ringing.
par.ref_ms = 1.5;

% Spike window. At 40 kHz, w_pre=20 = 0.5 ms before peak, w_post=44 =
% 1.1 ms after. This captures the full repolarization of pyramidal cells
% which matters for half-width and trough-to-peak classification.
% 60 samples total at 40 kHz = 1.5 ms window. Increase w_post to 60 if
% you see clipped repolarizations on wide-waveform units.
par.w_pre = 20;
par.w_post = 60;                     % increased from 44 for full pyramidal
                                     % repolarisation at 40 kHz
par.alignment_window = 10;

%% FILTER PARAMETERS
% The data is already bandpassed at 250-5000 Hz by extract_and_filter.m.
% wave_clus re-filters internally. We set 300-5000 Hz here:
%   - 300 Hz lowcut trims residual LFP bleed from the external 250 Hz edge
%   - 5000 Hz highcut matches the external filter so we don't lose
%     high-frequency content from the initial spike downstroke
% The cascade of two bandpass filters at similar cutoffs causes negligible
% waveform distortion for classification purposes.
par.detect_fmin = 300;
par.detect_fmax = 5000;              % raised from 3000 default
par.detect_order = 4;

par.sort_fmin = 300;
par.sort_fmax = 5000;                % raised from 3000 default
par.sort_order = 2;

%% SPC CLUSTERING PARAMETERS
par.mintemp = 0.00;
par.maxtemp = 0.251;
par.tempstep = 0.01;
par.SWCycles = 100;
par.KNearNeighb = 11;

% Minimum cluster size. Default (60) is too aggressive for sparse neurons.
% At 0.1 Hz over 45 min = ~270 spikes total. With Inf segment length,
% all 270 are available for clustering. 20 is a reasonable floor — any
% fewer and you can't meaningfully estimate a waveform template.
par.min_clus = 20;

par.randomseed = 0;                  % clock-based random seed
par.temp_plot = 'log';
par.c_ov = 0.7;                      % overlap coefficient for inclusion
par.elbow_min = 0.4;                 % regime border detection threshold

%% INTERPOLATION
% Cubic spline interpolation at 5x improves alignment precision.
% At 40 kHz base rate this gives 200 kHz effective resolution for
% peak alignment — 5 microsecond precision.
par.int_factor = 5;
par.interpolation = 'y';

%% FEATURE EXTRACTION
par.min_inputs = 10;
par.max_inputs = 0.75;               % proportion of max if < 1
par.scales = 4;                      % wavelet decomposition scales
par.features = 'wav';                % wavelets (better than PCA for
                                     % non-stationary waveform shapes)

%% FORCE MEMBERSHIP (template matching)
% After initial SPC clustering, unassigned spikes are force-assigned to
% the nearest cluster. This recovers spikes that fell outside the SPC
% temperature regime — important for sparse units where every spike counts.
par.template_sdnum = 3;              % max radius in std devs
par.template_k = 10;
par.template_k_min = 10;
par.template_type = 'center';        % distance to cluster center
par.force_feature = 'spk';           % match on full spike shape
par.force_auto = true;               % auto-force in batch mode

%% TEMPLATE MATCHING
% For recordings with >40k detected events, template matching kicks in
% to avoid SPC memory explosion. With Inf segments and 40 kHz, a 45-min
% recording at moderate MUA rates can hit this. Leave enabled.
par.match = 'y';
par.max_spk = 40000;
par.permut = 'y';                    % random subset for template building

%% TIME RANGE
par.tmax = 'all';                    % load entire recording
par.tmin = 0;

%% PLOTTING
par.cont_segment = true;
par.max_spikes_plot = 1000;
par.print2file = true;
par.cont_plot_samples = 100000;
par.to_plot_std = 1;
par.all_classes_ax = 'mean';
par.plot_feature_stats = false;

%% HISTOGRAM
par.nbins = 100;
par.bin_step = 1;

end
