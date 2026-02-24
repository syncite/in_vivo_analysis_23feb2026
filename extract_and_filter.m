function extract_and_filter(plx_file, varargin)
% EXTRACT_AND_FILTER  Extract .plx → per-channel .mat, then bandpass + CAR.
%
%   extract_and_filter()                         — opens file picker
%   extract_and_filter('path/to/file.plx')
%   extract_and_filter('path/to/file.plx', 'param', value, ...)
%
% This combines Plexon extraction and filtering into one step. Output files
% are ready to open directly in the wave_clus GUI.
%
% Parameters:
%   'sdk_path'      — path to Plexon Offline SDK [default: see below]
%   'freq_low'      — highpass cutoff in Hz      [default: 250]
%   'freq_high'     — lowpass cutoff in Hz        [default: 5000]
%   'filter_order'  — Butterworth order           [default: 3]
%   'channels'      — which channels to process (empty = all) [default: []]
%   'keep_raw'      — keep unfiltered channel files [default: false]
%   'max_trials'    — keep only first N strobed events (empty = all) [default: []]
%   'post_trial_padding_s' — seconds after Nth event to include [default: 5]
%
% Output structure:
%   <plxname>.mat                                    — master file (events, metadata)
%   <plxname>_channels/channel_<N>_filtered_CAR.mat  — filtered+CAR (single)
%                                                       ↑ open these in wave_clus
%
% After running this, open wave_clus GUI and load the _filtered_CAR.mat
% files for manual spike sorting. Then run analyse_peri_event.m on the
% wave_clus output.

    %% Parse parameters
    p = inputParser;
    addOptional(p, 'plx_file', '', @ischar);
    addParameter(p, 'sdk_path', ...
        'C:\Users\Madhu\Desktop\PlexonToMat\OmniPlex and MAP Offline SDK Bundle', @ischar);
    addParameter(p, 'freq_low', 250, @isnumeric);
    addParameter(p, 'freq_high', 5000, @isnumeric);
    addParameter(p, 'filter_order', 3, @isnumeric);
    addParameter(p, 'channels', [], @isnumeric);
    addParameter(p, 'keep_raw', false, @islogical);
    addParameter(p, 'max_trials', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    addParameter(p, 'post_trial_padding_s', 5, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    parse(p, plx_file, varargin{:});
    opts = p.Results;

    %% ========== PART 1: PLEXON EXTRACTION ==========

    addpath(genpath(opts.sdk_path));

    % Select .plx file
    if isempty(opts.plx_file)
        [PlexonDataFile, PlexonDataPath] = uigetfile( ...
            {'*.plx','Plexon files'; '*.*','All Files'}, ...
            'Select the Plexon data file');
        if isequal(PlexonDataFile, 0), return; end
        plx_full = fullfile(PlexonDataPath, PlexonDataFile);
    else
        plx_full = opts.plx_file;
    end

    [plxDir, plxName, ~] = fileparts(plx_full);

    % Read metadata
    [~, ~, FsPlexon, ~, ~, ~, ~, ~, ~, ~, ~, Duration, DateTime] = ...
        plx_information(plx_full);

    fprintf('\n========== EXTRACTION ==========\n');
    fprintf('File: %s\n', plx_full);
    fprintf('Sampling rate: %d Hz | Duration: %.1f s\n', FsPlexon, Duration);

    % Find active channels
    [~, samplecounts] = plx_adchan_samplecounts(plx_full);
    activeChannels = find(samplecounts > 0);
    nChannels = numel(activeChannels);
    fprintf('Active channels: %d\n', nChannels);

    % Create output directory (<plxname>_channels)
    outDir = fullfile(plxDir, [plxName '_channels']);
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    % Determine which channels to process
    if isempty(opts.channels)
        channels = 1:nChannels;
    else
        channels = opts.channels;
    end

    % Extract strobed events
    [EventsNo, Eventstime, EventTag] = plx_event_ts(plx_full, 257);
    fprintf('Strobed events: %d\n', EventsNo);

    cutoff_time_s = [];
    if ~isempty(opts.max_trials)
        if isempty(Eventstime)
            warning('No strobed events found; ignoring max_trials and processing full recording.');
        else
            n_req = max(1, round(opts.max_trials));
            n_use = min(n_req, numel(Eventstime));
            nth_event_time = Eventstime(n_use);
            cutoff_time_s = nth_event_time + opts.post_trial_padding_s;

            keep_evt = Eventstime <= cutoff_time_s;
            Eventstime = Eventstime(keep_evt);
            EventTag = EventTag(keep_evt);
            EventsNo = numel(Eventstime);

            fprintf(['Trial-limited extraction enabled: first %d events, cutoff at %.3f s ' ...
                '(Nth=%.3f s + %.3f s); keeping %d strobed events.\n'], ...
                n_use, cutoff_time_s, nth_event_time, opts.post_trial_padding_s, EventsNo);
        end
    end

    % Extract per-channel voltage
    firstADsamples = zeros(nChannels, 1);
    nSamples = 0;

    for ii = 1:nChannels
        chIdx = activeChannels(ii) - 1;  % Plexon SDK: 0-based
        [~, ~, ts, ~, ad] = plx_ad_v(plx_full, chIdx);
        firstADsamples(ii) = ts;

        ad = ad(:)';
        if ~isempty(cutoff_time_s)
            keep_n = floor((cutoff_time_s - ts) * FsPlexon) + 1;
            keep_n = min(max(keep_n, 0), numel(ad));
            ad = ad(1:keep_n);
        end

        nSamples = numel(ad);
        data = single(ad);  %#ok<NASGU>
        save(fullfile(outDir, sprintf('channel_%d.mat', ii)), 'data');
        fprintf('  Channel %d (Plexon ch %d): %d samples, firstAD = %.6f s\n', ...
            ii, activeChannels(ii), nSamples, ts);
    end

    % Print event tag summary
    unique_tags = unique(EventTag);
    fprintf('\nEvent tag summary:\n');
    tag_labels = struct('t32385','5V 1000ms', 't32281','0V catch', ...
        't32297','0.45V 1000ms', 't32329','1.1V 1000ms', ...
        't32321','5V 100ms', 't32417','5V 10ms');
    for ti = 1:numel(unique_tags)
        tag = unique_tags(ti);
        n = sum(EventTag == tag);
        field = sprintf('t%d', tag);
        if isfield(tag_labels, field)
            label = tag_labels.(field);
        else
            label = 'unknown';
        end
        fprintf('  Tag %d (%s): %d trials\n', tag, label, n);
    end

    % Save master file
    params = struct( ...
        'pipeline_step', 'extract_and_filter', ...
        'source_file', plx_full, ...
        'freq_low', opts.freq_low, ...
        'freq_high', opts.freq_high, ...
        'filter_order', opts.filter_order, ...
        'max_trials', opts.max_trials, ...
        'post_trial_padding_s', opts.post_trial_padding_s, ...
        'date_processed', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

    masterFile = fullfile(plxDir, [plxName '.mat']);
    save(masterFile, ...
        'firstADsamples', 'activeChannels', 'FsPlexon', ...
        'EventsNo', 'Eventstime', 'EventTag', ...
        'DateTime', 'Duration', 'nChannels', 'params', ...
        '-v7.3');
    fprintf('Master file: %s\n', masterFile);

    %% ========== PART 2: BANDPASS FILTER + CAR ==========

    fprintf('\n========== FILTERING ==========\n');
    Fs = FsPlexon;

    fprintf('Bandpass: %d–%d Hz (Butterworth order %d, effective %d after filtfilt)\n', ...
        opts.freq_low, opts.freq_high, opts.filter_order, opts.filter_order * 2);

    % Design filter
    Wn = [opts.freq_low, opts.freq_high] / (Fs / 2);
    [b, a] = butter(opts.filter_order, Wn, 'bandpass');

    % Pass 1: filter each channel, then align lengths for CAR
    fprintf('Pass 1: bandpass filtering...\n');
    nCh = numel(channels);
    filtered = cell(nCh, 1);
    ch_lengths = zeros(1, nCh);

    for ii = 1:nCh
        ch = channels(ii);
        tmp = load(fullfile(outDir, sprintf('channel_%d.mat', ch)), 'data');
        if isempty(tmp.data)
            error(['Channel %d has no samples after trial truncation. Increase max_trials ' ...
                'or disable max_trials.'], ch);
        end

        filt_data = filtfilt(b, a, double(tmp.data));
        filtered{ii} = filt_data;
        ch_lengths(ii) = numel(filt_data);

        fprintf('  Channel %d filtered\n', ch);
    end

    nSamples = min(ch_lengths);
    if any(ch_lengths ~= nSamples)
        warning('Channel lengths differ after truncation. Using first %d samples for CAR.', nSamples);
    end

    running_sum = zeros(1, nSamples);
    for ii = 1:nCh
        filtered{ii} = filtered{ii}(1:nSamples);
        running_sum = running_sum + filtered{ii};
    end

    % Common average reference
    car_signal = running_sum / nCh;
    clear running_sum;

    % Pass 2: subtract CAR and save
    fprintf('Pass 2: subtracting CAR and saving...\n');
    for ii = 1:nCh
        ch = channels(ii);
        data = single(filtered{ii} - car_signal);  %#ok<NASGU>
        filtered{ii} = [];  % free memory

        outFile = fullfile(outDir, sprintf('channel_%d_filtered_CAR.mat', ch));
        save(outFile, 'data');
        fprintf('  Saved: channel_%d_filtered_CAR.mat\n', ch);
    end

    %% Clean up raw channel files
    if ~opts.keep_raw
        fprintf('\nDeleting raw channel files (not needed after filtering)...\n');
        for ii = 1:nCh
            ch = channels(ii);
            rawFile = fullfile(outDir, sprintf('channel_%d.mat', ch));
            if exist(rawFile, 'file')
                delete(rawFile);
                fprintf('  Deleted: channel_%d.mat\n', ch);
            end
        end
    end

    %% Done
    fprintf('\n========== DONE ==========\n');
    fprintf('Channel files: %s\n', outDir);
    fprintf('Master file:   %s\n', masterFile);
    fprintf('\nNext step:\n');
    fprintf('  1. Open wave_clus GUI\n');
    fprintf('  2. Load the channel_*_filtered_CAR.mat files from:\n');
    fprintf('     %s\n', outDir);
    fprintf('  3. Set parameters: sr=%d, stdFactor=4, segments_length=Inf\n', Fs);
    fprintf('  4. Sort each channel manually\n');
    fprintf('  5. Run analyse_peri_event on the wave_clus output\n');
end
