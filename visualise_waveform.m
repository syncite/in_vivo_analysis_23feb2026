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
%   3) picks random LED trials where both units have >=1 spike in window
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
        'Eventstime', 'EventTag', 'nChannels', 'FsPlexon', ...
        'firstADsamples', 'firstADsample');

    all_event_times = master.Eventstime(:)' - opts.timing_correction;
    all_event_tags = master.EventTag(:)';

    led_mask = ismember(all_event_tags, opts.led_tags);
    if ~any(led_mask)
        warning('No requested LED tags found. Falling back to all non-catch tags.');
        led_mask = all_event_tags ~= 32281;
    end
    led_events = all_event_times(led_mask);
    led_event_idx = find(led_mask);
    if isempty(led_events)
        error('No LED events found in master file.');
    end

    if isempty(opts.seed)
        rng('shuffle');
    else
        rng(opts.seed);
    end

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

        first_sample_s = get_first_sample_seconds(master, ch);

        cluster_ids = [];
        spike_by_unit = {};
        if ~isempty(wc_file)
            wc = load(wc_file, 'cluster_class');
            cc = wc.cluster_class;
            cluster_ids = unique(cc(:, 1));
            cluster_ids = cluster_ids(cluster_ids > 0);
            cc_times_s = cc(:, 2)' / 1000;
            [cc_times_s, ref_mode] = normalize_spike_time_reference( ...
                cc_times_s, first_sample_s, all_event_times);
            fprintf('  Channel %d: spike time reference = %s\n', ch, ref_mode);

            spike_by_unit = cell(numel(cluster_ids), 1);
            for ui = 1:numel(cluster_ids)
                mask = cc(:, 1) == cluster_ids(ui);
                spike_by_unit{ui} = cc_times_s(mask); % normalized to event base (s)
            end
        end

        n_units = numel(cluster_ids);
        if n_units < 2
            warning(['Channel %d: found %d unit(s). Need at least 2 units to ' ...
                'select trials where both units spike. Skipping.'], ch, n_units);
            continue;
        end

        % Use the first two sorted clusters as the required "both units".
        required_units = 1:2;
        eligible = false(1, numel(led_events));
        for ei = 1:numel(led_events)
            t0 = led_events(ei);

            event_sample = round((t0 - first_sample_s) * Fs) + 1;
            idx = event_sample + sample_offsets;
            if idx(1) < 1 || idx(end) > numel(trace)
                continue;
            end

            both_units_spike = true;
            for ui = required_units
                rel_ms = (spike_by_unit{ui} - t0) * 1000;
                in_win = rel_ms >= -opts.window_ms & rel_ms <= opts.window_ms;
                if ~any(in_win)
                    both_units_spike = false;
                    break;
                end
            end
            eligible(ei) = both_units_spike;
        end

        eligible_idx = find(eligible);
        if isempty(eligible_idx)
            warning(['Channel %d: no LED trials where both required units have ' ...
                'spikes within +/- %.1f ms. Skipping.'], ch, opts.window_ms);
            continue;
        end

        n_trials = min(max(1, round(opts.n_trials)), numel(eligible_idx));
        if numel(eligible_idx) < round(opts.n_trials)
            warning('Channel %d: only %d qualifying trials found (requested %d).', ...
                ch, numel(eligible_idx), round(opts.n_trials));
        end
        chosen = eligible_idx(randperm(numel(eligible_idx), n_trials));
        trial_times = led_events(chosen);
        trial_ids = led_event_idx(chosen);

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
            if valid_trial(ti)
                y_min = min(seg(ti, :));
                y_max = max(seg(ti, :));
                y_rng = y_max - y_min;
                if y_rng < eps
                    y_pad = max(1e-6, abs(y_max) * 0.05);
                else
                    y_pad = 0.02 * y_rng;
                end
                ylim(ax_wave, [y_min - y_pad, y_max + y_pad]);
            end
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

function first_sample_s = get_first_sample_seconds(master, ch)
    first_sample_s = 0;

    if isfield(master, 'firstADsamples')
        vals = double(master.firstADsamples(:));
    elseif isfield(master, 'firstADsample')
        vals = double(master.firstADsample(:));
    else
        return;
    end

    if ch >= 1 && ch <= numel(vals)
        first_sample_s = vals(ch);
    end
end

function [spike_s, mode] = normalize_spike_time_reference(spike_s_raw, first_sample_s, event_s)
    spike_rel = double(spike_s_raw(:)');
    spike_abs = spike_rel + first_sample_s;

    if isempty(spike_rel) || isempty(event_s) || abs(first_sample_s) < eps
        spike_s = spike_rel;
        mode = 'as-is';
        return;
    end

    lo = min(event_s) - 2;
    hi = max(event_s) + 2;
    frac_rel = mean(spike_rel >= lo & spike_rel <= hi);
    frac_abs = mean(spike_abs >= lo & spike_abs <= hi);

    if frac_abs > frac_rel + 0.05
        spike_s = spike_abs;
        mode = 'offset+firstAD';
        return;
    end
    if frac_rel > frac_abs + 0.05
        spike_s = spike_rel;
        mode = 'as-is';
        return;
    end

    med_ev = median(event_s);
    score_rel = abs(median(spike_rel) - med_ev);
    score_abs = abs(median(spike_abs) - med_ev);

    if score_abs < score_rel
        spike_s = spike_abs;
        mode = 'offset+firstAD';
    else
        spike_s = spike_rel;
        mode = 'as-is';
    end
end
