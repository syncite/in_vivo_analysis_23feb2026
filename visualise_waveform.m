function visualise_waveform(master_mat_file, varargin)
% VISUALISE_WAVEFORM  Plot filtered trace + unit raster around random LED trials.
%
%   visualise_waveform()                         % file picker for master .mat
%   visualise_waveform('path/to/master.mat')
%   visualise_waveform('path/to/master.mat', 'channels', [3 8 12])
%
% For each requested channel, this function:
%   1) loads channel_<N>_filtered_CAR.mat
%   2) loads spike times from times_channel_<N>_*.mat (wave_clus output)
%   3) picks 3 random LED trials
%   4) plots -10 ms to +10 ms around each onset (3 columns)
%      - row 1: filtered waveform
%      - row 2: raster (different colors per unit)

    p = inputParser;
    addOptional(p, 'master_mat_file', '', @ischar);
    addParameter(p, 'channels', [], @isnumeric);
    addParameter(p, 'n_trials', 3, @isnumeric);
    addParameter(p, 'window_ms', 10, @isnumeric);
    addParameter(p, 'timing_correction', 0.02295, @isnumeric);
    addParameter(p, 'led_tags', [32385 32297 32329 32321 32417], @isnumeric);
    addParameter(p, 'seed', [], @(x) isempty(x) || isnumeric(x));
    parse(p, master_mat_file, varargin{:});
    opts = p.Results;

    if isempty(opts.master_mat_file)
        [matFile, matPath] = uigetfile('*.mat', ...
            'Select the master .mat file (from extract_and_filter)');
        if isequal(matFile, 0), return; end
        opts.master_mat_file = fullfile(matPath, matFile);
    end

    master = load(opts.master_mat_file, ...
        'Eventstime', 'EventTag', 'nChannels', 'FsPlexon', 'firstADsamples');

    all_event_times = master.Eventstime(:)' - opts.timing_correction;
    all_event_tags = master.EventTag(:)';

    led_mask = ismember(all_event_tags, opts.led_tags);
    if ~any(led_mask)
        warning('No requested LED tags found. Falling back to all non-catch tags.');
        led_mask = all_event_tags ~= 32281;
    end
    led_events = all_event_times(led_mask);
    if isempty(led_events)
        error('No LED events found in master file.');
    end

    if isempty(opts.seed)
        rng('shuffle');
    else
        rng(opts.seed);
    end

    n_trials = min(max(1, round(opts.n_trials)), numel(led_events));
    trial_ids = randperm(numel(led_events), n_trials);
    trial_times = led_events(trial_ids);

    [masterPath, masterName, ~] = fileparts(opts.master_mat_file);
    channelDir = fullfile(masterPath, [masterName '_channels']);
    if ~exist(channelDir, 'dir')
        channelDir = fullfile(masterPath, masterName);
    end

    if isempty(opts.channels)
        channels = 1:master.nChannels;
    else
        channels = opts.channels(:)';
    end

    Fs = master.FsPlexon;
    n_before = round((opts.window_ms / 1000) * Fs);
    n_after = round((opts.window_ms / 1000) * Fs);
    sample_offsets = -n_before:n_after;
    t_ms = sample_offsets / Fs * 1000;

    for ch = channels
        filt_candidates = {
            fullfile(channelDir, sprintf('channel_%d_filtered_CAR.mat', ch))
            fullfile(channelDir, sprintf('channel_%d.mat', ch))
        };
        filt_file = pick_first_existing(filt_candidates);
        if isempty(filt_file)
            warning('Channel %d: filtered file not found. Skipping.', ch);
            continue;
        end

        wc_candidates = {
            fullfile(channelDir, sprintf('times_channel_%d_filtered_CAR.mat', ch))
            fullfile(channelDir, sprintf('times_channel_%d_filtered_CARfiltered.mat', ch))
            fullfile(masterPath, sprintf('times_channel_%d.mat', ch))
            fullfile(channelDir, sprintf('times_channel_%d.mat', ch))
        };
        wc_file = pick_first_existing(wc_candidates);

        td = load(filt_file, 'data');
        trace = double(td.data(:)');

        first_sample_s = 0;
        if isfield(master, 'firstADsamples') && numel(master.firstADsamples) >= ch
            first_sample_s = double(master.firstADsamples(ch));
        end

        cluster_ids = [];
        spike_by_unit = {};
        if ~isempty(wc_file)
            wc = load(wc_file, 'cluster_class');
            cc = wc.cluster_class;
            cluster_ids = unique(cc(:, 1));
            cluster_ids = cluster_ids(cluster_ids > 0);
            spike_by_unit = cell(numel(cluster_ids), 1);
            for ui = 1:numel(cluster_ids)
                mask = cc(:, 1) == cluster_ids(ui);
                spike_by_unit{ui} = cc(mask, 2) / 1000; % wave_clus stores ms
            end
        end

        seg = nan(n_trials, numel(sample_offsets));
        valid_trial = false(1, n_trials);
        for ti = 1:n_trials
            event_sample = round((trial_times(ti) - first_sample_s) * Fs) + 1;
            idx = event_sample + sample_offsets;
            if idx(1) >= 1 && idx(end) <= numel(trace)
                seg(ti, :) = trace(idx);
                valid_trial(ti) = true;
            end
        end

        fig = figure('Color', 'w', 'Position', [80 120 1500 700]);
        tiled = tiledlayout(2, n_trials, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tiled, sprintf('Channel %d | %d random LED trials | window = +/- %d ms', ...
            ch, n_trials, round(opts.window_ms)));

        if any(valid_trial)
            y_min = min(seg(valid_trial, :), [], 'all');
            y_max = max(seg(valid_trial, :), [], 'all');
            y_pad = max(1, 0.05 * (y_max - y_min + eps));
            wave_ylim = [y_min - y_pad, y_max + y_pad];
        else
            wave_ylim = [];
        end

        n_units = numel(cluster_ids);
        unit_colors = lines(max(n_units, 1));

        for ti = 1:n_trials
            ax_wave = nexttile(tiled, ti);
            hold(ax_wave, 'on');
            if valid_trial(ti)
                plot(ax_wave, t_ms, seg(ti, :), 'k', 'LineWidth', 1);
            else
                text(ax_wave, 0, 0, 'Window exceeds data bounds', ...
                    'HorizontalAlignment', 'center', 'Color', [0.4 0.4 0.4]);
            end
            xline(ax_wave, 0, '--', 'Color', [0.8 0 0], 'LineWidth', 1);
            xlim(ax_wave, [-opts.window_ms opts.window_ms]);
            if ~isempty(wave_ylim), ylim(ax_wave, wave_ylim); end
            grid(ax_wave, 'on');
            if ti == 1
                ylabel(ax_wave, 'Filtered signal');
            end
            title(ax_wave, sprintf('LED trial #%d', trial_ids(ti)));

            ax_raster = nexttile(tiled, n_trials + ti);
            hold(ax_raster, 'on');
            xline(ax_raster, 0, '--', 'Color', [0.8 0 0], 'LineWidth', 1);
            xlim(ax_raster, [-opts.window_ms opts.window_ms]);
            grid(ax_raster, 'on');
            xlabel(ax_raster, 'Time from LED onset (ms)');

            if n_units == 0
                text(ax_raster, 0, 0.5, 'No sorted units found for this channel', ...
                    'HorizontalAlignment', 'center', 'Color', [0.4 0.4 0.4]);
                ylim(ax_raster, [0 1]);
                set(ax_raster, 'YTick', []);
            else
                for ui = 1:n_units
                    rel_ms = (spike_by_unit{ui} - trial_times(ti)) * 1000;
                    rel_ms = rel_ms(rel_ms >= -opts.window_ms & rel_ms <= opts.window_ms);
                    if ~isempty(rel_ms)
                        scatter(ax_raster, rel_ms, ui * ones(size(rel_ms)), ...
                            22, unit_colors(ui, :), 'filled', 'MarkerFaceAlpha', 0.8);
                    end
                end
                ylim(ax_raster, [0.5 n_units + 0.5]);
                yticks(ax_raster, 1:n_units);
                yticklabels(ax_raster, arrayfun(@(x) sprintf('Cl%d', x), ...
                    cluster_ids, 'UniformOutput', false));
                if ti == 1
                    ylabel(ax_raster, 'Unit');
                end
            end
        end

        set(fig, 'Name', sprintf('visualise_waveform_ch%d', ch));
    end
end

function f = pick_first_existing(candidates)
    f = '';
    for i = 1:numel(candidates)
        if exist(candidates{i}, 'file')
            f = candidates{i};
            return;
        end
    end
end
