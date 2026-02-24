function visualise_waveform(master_mat_file, varargin)
% VISUALISE_WAVEFORM  Plot filtered trace with raster strips around random LED trials.
%
%   visualise_waveform()                         % file picker for master .mat
%   visualise_waveform('path/to/master.mat')
%   visualise_waveform('path/to/master.mat', 'channels', [3 8 12])
%
% For each requested channel, this function:
%   1) loads channel_<N>_filtered_CAR.mat
%   2) loads spike times from times_channel_<N>_*.mat (wave_clus output)
%   3) makes one figure with trials where unit 2 spikes in first 10 ms
%   4) makes one figure with trials where unit 2 does not spike in first 10 ms
%   5) makes one overlay figure of unit 2 spike-aligned traces:
%      early spikes (0-10 ms from LED) vs all other unit 2 spikes.
%      Trial panels contain waveform (top) + raster strip (bottom).

    % Editable plot windows (ms): [start end]
    top_window_ms = [-100 100];
    bottom_window_ms = [-3 10];
    unit2_early_window_ms = [0 10];
    spike_overlay_window_ms = [-1 2];

    p = inputParser;
    addOptional(p, 'master_mat_file', '', @ischar);
    addParameter(p, 'channels', [], @isnumeric);
    addParameter(p, 'n_trials', 3, @isnumeric);
    addParameter(p, 'window_ms', [], @(x) isempty(x) || isnumeric(x));
    addParameter(p, 'led_tags', [32385 32297 32329 32321 32417], @isnumeric);
    addParameter(p, 'seed', [], @(x) isempty(x) || isnumeric(x));
    parse(p, master_mat_file, varargin{:});
    opts = p.Results;

    % Backward compatibility: if provided, use symmetric window for top row.
    if ~isempty(opts.window_ms)
        top_window_ms = [-abs(opts.window_ms) abs(opts.window_ms)];
    end
    top_window_ms = sort(top_window_ms(:)');
    bottom_window_ms = sort(bottom_window_ms(:)');

    if isempty(opts.master_mat_file)
        [matFile, matPath] = uigetfile('*.mat', ...
            'Select the master .mat file (from extract_and_filter)');
        if isequal(matFile, 0), return; end
        opts.master_mat_file = fullfile(matPath, matFile);
    end

    master = load(opts.master_mat_file, ...
        'Eventstime', 'EventTag', 'nChannels', 'FsPlexon', ...
        'firstADsamples', 'firstADsample');

    all_event_times = master.Eventstime(:)';
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
    sample_offsets_top = round((top_window_ms(1) / 1000) * Fs) : ...
                         round((top_window_ms(2) / 1000) * Fs);
    sample_offsets_bottom = round((bottom_window_ms(1) / 1000) * Fs) : ...
                            round((bottom_window_ms(2) / 1000) * Fs);
    t_top_ms = sample_offsets_top / Fs * 1000;
    t_bottom_ms = sample_offsets_bottom / Fs * 1000;

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
                'evaluate unit 2 trial groups. Skipping.'], ch, n_units);
            continue;
        end

        unit2_idx = 2;
        unit2_spikes = spike_by_unit{unit2_idx};
        has_unit2_early = false(1, numel(led_events));
        event_in_bounds = false(1, numel(led_events));

        for ei = 1:numel(led_events)
            t0 = led_events(ei);

            event_sample = round((t0 - first_sample_s) * Fs) + 1;
            idx_top = event_sample + sample_offsets_top;
            idx_bottom = event_sample + sample_offsets_bottom;
            event_in_bounds(ei) = idx_top(1) >= 1 && idx_top(end) <= numel(trace) && ...
                                  idx_bottom(1) >= 1 && idx_bottom(end) <= numel(trace);

            rel2_ms = (unit2_spikes - t0) * 1000;
            has_unit2_early(ei) = any(rel2_ms >= unit2_early_window_ms(1) & ...
                                      rel2_ms <= unit2_early_window_ms(2));
        end

        group_hit = find(event_in_bounds & has_unit2_early);
        group_miss = find(event_in_bounds & ~has_unit2_early);

        if isempty(group_hit) && isempty(group_miss)
            warning(['Channel %d: no LED trials in-bounds for requested windows ' ...
                '(top %.1f to %.1f ms, bottom %.1f to %.1f ms). Skipping.'], ...
                ch, top_window_ms(1), top_window_ms(2), ...
                bottom_window_ms(1), bottom_window_ms(2));
            continue;
        end

        unit_colors = lines(max(n_units, 1));

        n_req = max(1, round(opts.n_trials));
        if isempty(group_hit)
            warning('Channel %d: no trials where unit 2 spikes in first %.1f ms.', ...
                ch, unit2_early_window_ms(2));
        else
            if numel(group_hit) < n_req
                warning('Channel %d: only %d unit2-early trials found (requested %d).', ...
                    ch, numel(group_hit), n_req);
            end
            chosen_hit = group_hit(randperm(numel(group_hit), min(n_req, numel(group_hit))));
            plot_trial_group_figure(ch, 'Unit 2 spike in first 10 ms', ...
                led_events(chosen_hit), led_event_idx(chosen_hit), ...
                trace, first_sample_s, Fs, ...
                sample_offsets_top, sample_offsets_bottom, ...
                t_top_ms, t_bottom_ms, cluster_ids, spike_by_unit, unit_colors, ...
                top_window_ms, bottom_window_ms);
        end

        if isempty(group_miss)
            warning('Channel %d: no trials where unit 2 is silent in first %.1f ms.', ...
                ch, unit2_early_window_ms(2));
        else
            if numel(group_miss) < n_req
                warning('Channel %d: only %d unit2-no-early trials found (requested %d).', ...
                    ch, numel(group_miss), n_req);
            end
            chosen_miss = group_miss(randperm(numel(group_miss), min(n_req, numel(group_miss))));
            plot_trial_group_figure(ch, 'Unit 2 no spike in first 10 ms', ...
                led_events(chosen_miss), led_event_idx(chosen_miss), ...
                trace, first_sample_s, Fs, ...
                sample_offsets_top, sample_offsets_bottom, ...
                t_top_ms, t_bottom_ms, cluster_ids, spike_by_unit, unit_colors, ...
                top_window_ms, bottom_window_ms);
        end

        plot_unit2_overlay_figure(ch, unit2_spikes, led_events, ...
            unit2_early_window_ms, spike_overlay_window_ms, ...
            trace, first_sample_s, Fs);
    end
end

function plot_wave_raster_tile(host_ax, t_ms, seg, is_valid, trial_time_s, ...
    cluster_ids, spike_by_unit, unit_colors, win_ms, panel_title, ...
    show_left_labels, show_x_label)

    host_pos = get(host_ax, 'Position');
    delete(host_ax);

    raster_frac = 0.30;
    gap_frac = 0.05;
    wave_frac = 1 - raster_frac - gap_frac;
    wave_pos = [host_pos(1), host_pos(2) + host_pos(4)*(raster_frac + gap_frac), ...
                host_pos(3), host_pos(4)*wave_frac];
    raster_pos = [host_pos(1), host_pos(2), host_pos(3), host_pos(4)*raster_frac];

    ax_wave = axes('Position', wave_pos);
    hold(ax_wave, 'on');
    if is_valid
        plot(ax_wave, t_ms, seg, 'k', 'LineWidth', 1);
        y_min = min(seg);
        y_max = max(seg);
        y_rng = y_max - y_min;
        if y_rng < eps
            y_pad = max(1e-6, abs(y_max) * 0.05);
        else
            y_pad = 0.02 * y_rng;
        end
        ylim(ax_wave, [y_min - y_pad, y_max + y_pad]);
    else
        text(ax_wave, mean(win_ms), 0, 'Window exceeds data bounds', ...
            'HorizontalAlignment', 'center', 'Color', [0.4 0.4 0.4]);
        ylim(ax_wave, [-1 1]);
    end
    xline(ax_wave, 0, '--', 'Color', [0.8 0 0], 'LineWidth', 1);
    xlim(ax_wave, win_ms);
    grid(ax_wave, 'on');
    box(ax_wave, 'on');
    set(ax_wave, 'XTickLabel', []);
    if show_left_labels
        ylabel(ax_wave, 'Filtered signal');
    else
        set(ax_wave, 'YTickLabel', []);
    end
    title(ax_wave, panel_title, 'FontSize', 9);

    ax_raster = axes('Position', raster_pos);
    hold(ax_raster, 'on');
    xline(ax_raster, 0, '--', 'Color', [0.8 0 0], 'LineWidth', 1);
    xlim(ax_raster, win_ms);
    grid(ax_raster, 'on');
    box(ax_raster, 'on');

    n_units = numel(cluster_ids);
    if n_units == 0
        text(ax_raster, mean(win_ms), 0.5, 'No sorted units found for this channel', ...
            'HorizontalAlignment', 'center', 'Color', [0.4 0.4 0.4]);
        ylim(ax_raster, [0 1]);
        set(ax_raster, 'YTick', []);
    else
        for ui = 1:n_units
            rel_ms = (spike_by_unit{ui} - trial_time_s) * 1000;
            rel_ms = rel_ms(rel_ms >= win_ms(1) & rel_ms <= win_ms(2));
            if ~isempty(rel_ms)
                scatter(ax_raster, rel_ms, ui * ones(size(rel_ms)), ...
                    22, unit_colors(ui, :), 'filled', 'MarkerFaceAlpha', 0.8);
            end
        end
        ylim(ax_raster, [0.5 n_units + 0.5]);
        yticks(ax_raster, 1:n_units);
        yticklabels(ax_raster, arrayfun(@(x) sprintf('Cl%d', x), ...
            cluster_ids, 'UniformOutput', false));
    end

    if show_left_labels
        ylabel(ax_raster, 'Unit');
    else
        set(ax_raster, 'YTickLabel', []);
    end
    if show_x_label
        xlabel(ax_raster, 'Time from LED onset (ms)');
    else
        set(ax_raster, 'XTickLabel', []);
    end
end

function plot_trial_group_figure(ch, group_label, trial_times, trial_ids, ...
    trace, first_sample_s, Fs, sample_offsets_top, sample_offsets_bottom, ...
    t_top_ms, t_bottom_ms, cluster_ids, spike_by_unit, unit_colors, ...
    top_window_ms, bottom_window_ms)

    n_trials = numel(trial_times);
    [seg_top, valid_top] = extract_event_segments(trace, first_sample_s, Fs, ...
        trial_times, sample_offsets_top);
    [seg_bottom, valid_bottom] = extract_event_segments(trace, first_sample_s, Fs, ...
        trial_times, sample_offsets_bottom);

    fig = figure('Color', 'w', 'Position', [80 120 1500 820]);
    tiled = tiledlayout(2, n_trials, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tiled, sprintf(['Channel %d | %s | %d random trials | ' ...
        'top %.1f to %.1f ms | bottom %.1f to %.1f ms'], ...
        ch, group_label, n_trials, ...
        top_window_ms(1), top_window_ms(2), ...
        bottom_window_ms(1), bottom_window_ms(2)));

    for ti = 1:n_trials
        host_top = nexttile(tiled, ti);
        plot_wave_raster_tile(host_top, t_top_ms, seg_top(ti, :), valid_top(ti), ...
            trial_times(ti), cluster_ids, spike_by_unit, unit_colors, ...
            top_window_ms, sprintf('LED trial #%d', trial_ids(ti)), ...
            ti == 1, true);

        host_bottom = nexttile(tiled, n_trials + ti);
        plot_wave_raster_tile(host_bottom, t_bottom_ms, seg_bottom(ti, :), ...
            valid_bottom(ti), trial_times(ti), cluster_ids, spike_by_unit, ...
            unit_colors, bottom_window_ms, ...
            sprintf('LED trial #%d (zoom)', trial_ids(ti)), ...
            ti == 1, true);
    end

    set(fig, 'Name', sprintf('visualise_waveform_ch%d_%s', ...
        ch, regexprep(lower(group_label), '[^a-z0-9]+', '_')));
end

function plot_unit2_overlay_figure(ch, unit2_spikes, led_events, ...
    early_window_ms, overlay_window_ms, trace, first_sample_s, Fs)

    early_mask = spike_mask_in_windows(unit2_spikes, led_events, early_window_ms / 1000);
    spikes_early = unit2_spikes(early_mask);
    spikes_other = unit2_spikes(~early_mask);

    sample_offsets_overlay = round((overlay_window_ms(1) / 1000) * Fs) : ...
                             round((overlay_window_ms(2) / 1000) * Fs);
    t_overlay_ms = sample_offsets_overlay / Fs * 1000;

    snippets_early = extract_spike_snippets(trace, spikes_early, ...
        first_sample_s, Fs, sample_offsets_overlay);
    snippets_other = extract_spike_snippets(trace, spikes_other, ...
        first_sample_s, Fs, sample_offsets_overlay);

    fig = figure('Color', 'w', 'Position', [120 170 920 520]);
    ax = axes('Parent', fig);
    hold(ax, 'on');

    col_early = [0.86 0.20 0.20];
    col_other = [0.15 0.45 0.85];
    line_col_early = blend_with_white(col_early, 0.30);
    line_col_other = blend_with_white(col_other, 0.30);

    for si = 1:size(snippets_other, 1)
        plot(ax, t_overlay_ms, snippets_other(si, :), 'Color', line_col_other, ...
            'LineWidth', 0.8);
    end
    for si = 1:size(snippets_early, 1)
        plot(ax, t_overlay_ms, snippets_early(si, :), 'Color', line_col_early, ...
            'LineWidth', 0.8);
    end

    if ~isempty(snippets_other)
        plot(ax, t_overlay_ms, mean(snippets_other, 1), 'Color', col_other, 'LineWidth', 2);
    end
    if ~isempty(snippets_early)
        plot(ax, t_overlay_ms, mean(snippets_early, 1), 'Color', col_early, 'LineWidth', 2);
    end

    h_other = plot(ax, nan, nan, 'Color', col_other, 'LineWidth', 2);
    h_early = plot(ax, nan, nan, 'Color', col_early, 'LineWidth', 2);

    xline(ax, 0, '--', 'Color', [0.1 0.1 0.1], 'LineWidth', 1);
    xlim(ax, overlay_window_ms);
    xlabel(ax, 'Time from Unit 2 spike (ms)');
    ylabel(ax, 'Filtered signal');
    grid(ax, 'on');
    box(ax, 'on');
    title(ax, sprintf(['Channel %d | Unit 2 spike overlays | early 0-10 ms from LED: n=%d | ' ...
        'other: n=%d'], ch, size(snippets_early, 1), size(snippets_other, 1)));
    legend(ax, [h_other h_early], ...
        {sprintf('Other unit 2 spikes (n=%d)', size(snippets_other, 1)), ...
        sprintf('Early unit 2 spikes (n=%d)', size(snippets_early, 1))}, ...
        'Location', 'best');
    set(fig, 'Name', sprintf('visualise_waveform_ch%d_unit2_overlay', ch));
end

function [segments, valid] = extract_event_segments(trace, first_sample_s, Fs, ...
    event_times_s, sample_offsets)

    n_evt = numel(event_times_s);
    segments = nan(n_evt, numel(sample_offsets));
    valid = false(1, n_evt);

    for ei = 1:n_evt
        event_sample = round((event_times_s(ei) - first_sample_s) * Fs) + 1;
        idx = event_sample + sample_offsets;
        if idx(1) >= 1 && idx(end) <= numel(trace)
            segments(ei, :) = trace(idx);
            valid(ei) = true;
        end
    end
end

function snippets = extract_spike_snippets(trace, spike_times_s, first_sample_s, Fs, sample_offsets)
    snippets = zeros(0, numel(sample_offsets));
    for si = 1:numel(spike_times_s)
        spike_sample = round((spike_times_s(si) - first_sample_s) * Fs) + 1;
        idx = spike_sample + sample_offsets;
        if idx(1) >= 1 && idx(end) <= numel(trace)
            snippets(end+1, :) = trace(idx); %#ok<AGROW>
        end
    end
end

function mask = spike_mask_in_windows(spike_s, event_times_s, win_s)
    mask = false(size(spike_s));
    if isempty(spike_s) || isempty(event_times_s)
        return;
    end

    for ei = 1:numel(event_times_s)
        t0 = event_times_s(ei);
        mask = mask | (spike_s >= t0 + win_s(1) & spike_s <= t0 + win_s(2));
    end
end

function col = blend_with_white(base_col, alpha_val)
    col = alpha_val * base_col + (1 - alpha_val) * [1 1 1];
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
