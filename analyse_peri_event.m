function analyse_peri_event(master_mat_file, varargin)
% ANALYSE_PERI_EVENT  Rasters, PSTHs, z-scored heatmaps, waveform classification.
%
%   analyse_peri_event()                         — opens file picker
%   analyse_peri_event('path/to/master.mat')
%   analyse_peri_event('path/to/master.mat', 'param', value, ...)
%
% Run this AFTER you've manually sorted in wave_clus. Produces:
%   1. Population summary per event tag (rasters, PSTHs, heatmap)
%   2. Waveform classification scatter (half-width vs trough-to-peak time)
%   3. Per-unit figure (waveform, scatter position, raster, z-score + 5V rate)
%   4. Per-unit LED waveform comparison (pre-LED vs early-LED)
%
% Parameters:
%   'channels'           — channel indices [default: all]
%   'window'             — [pre post] in seconds [default: [-1 2]]
%   'bin_size'           — PSTH bin width in seconds [default: 0.02]
%   'timing_correction'  — strobe latency in seconds [default: 0.02295]
%   'baseline_window'    — [start end] for z-score baseline [default: [-1 0]]
%   'heatmap_bin'        — bin width for heatmap [default: 0.1]
%   'min_spikes'         — minimum spikes to include a unit [default: 50]
%   'event_tags'         — specific tags to analyse (empty = all) [default: []]
%   'primary_tag'        — tag used for per-unit rasters/heatmaps [default: 32385]
%
% Event tag decoder (Anesthetised_590nm_v1 protocol):
%   32385 → 5V, 1000ms    (Stage 1 — full power)
%   32281 → 0V, 1000ms    (Catch — negative control)
%   32297 → 0.450V, 1000ms
%   32329 → 1.100V, 1000ms
%   32321 → 5V, 100ms
%   32417 → 5V, 10ms

    %% Parse parameters
    p = inputParser;
    addOptional(p, 'master_mat_file', '', @ischar);
    addParameter(p, 'channels', [], @isnumeric);
    addParameter(p, 'window', [-1 2], @isnumeric);
    addParameter(p, 'bin_size', 0.02, @isnumeric);
    addParameter(p, 'timing_correction', 0.02295, @isnumeric);
    addParameter(p, 'baseline_window', [-1 0], @isnumeric);
    addParameter(p, 'heatmap_bin', 0.1, @isnumeric);
    addParameter(p, 'min_spikes', 50, @isnumeric);
    addParameter(p, 'event_tags', [], @isnumeric);
    addParameter(p, 'primary_tag', 32385, @isnumeric);
    parse(p, master_mat_file, varargin{:});
    opts = p.Results;

    %% Load master file
    if isempty(opts.master_mat_file)
        [matFile, matPath] = uigetfile('*.mat', ...
            'Select the master .mat file (from extract_and_filter)');
        if isequal(matFile, 0), return; end
        opts.master_mat_file = fullfile(matPath, matFile);
    end

    master = load(opts.master_mat_file, 'Eventstime', 'EventTag', ...
        'nChannels', 'FsPlexon', 'firstADsamples', 'firstADsample');

    Fs = master.FsPlexon;
    all_event_times = master.Eventstime(:)' - opts.timing_correction;
    all_event_tags  = master.EventTag(:)';

    [masterPath, masterName, ~] = fileparts(opts.master_mat_file);
    channelDir = fullfile(masterPath, [masterName '_channels']);
    if ~exist(channelDir, 'dir')
        channelDir = fullfile(masterPath, masterName);
    end

    if isempty(opts.channels)
        channels = 1:master.nChannels;
    else
        channels = opts.channels;
    end

    outDir = fullfile(masterPath, 'peri_event_results');
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    if isempty(opts.event_tags)
        unique_tags = unique(all_event_tags);
    else
        unique_tags = opts.event_tags;
    end

    fprintf('Master file: %s\n', opts.master_mat_file);
    fprintf('Channel dir: %s\n', channelDir);
    fprintf('Timing correction: %.4f s\n', opts.timing_correction);
    fprintf('Window: [%.1f, %.1f] s | PSTH bin: %.0f ms | Heatmap bin: %.0f ms\n', ...
        opts.window(1), opts.window(2), opts.bin_size*1000, opts.heatmap_bin*1000);

    %% Tag labels
    tag_info = containers.Map('KeyType', 'int32', 'ValueType', 'char');
    tag_info(int32(32385)) = '5V 1000ms';
    tag_info(int32(32281)) = '0V catch';
    tag_info(int32(32297)) = '0.45V 1000ms';
    tag_info(int32(32329)) = '1.1V 1000ms';
    tag_info(int32(32321)) = '5V 100ms';
    tag_info(int32(32417)) = '5V 10ms';

    % Trial groups used in per-unit panels and waveform comparisons
    fiveV_tags = [32385 32321 32417];
    fiveV_events = all_event_times(ismember(all_event_tags, fiveV_tags));
    led_events = all_event_times(all_event_tags ~= 32281);  % exclude catch trials

    %% Load wave_clus output: spike times + waveforms
    units = struct('label', {}, 'spike_times_s', {}, 'channel', {}, ...
        'cluster', {}, 'waveforms', {}, 'avg_waveform', {}, ...
        'half_width_ms', {}, 'trough_to_peak_ms', {}, 'n_spikes', {});

    for ch = channels
        candidates = {
            fullfile(channelDir, sprintf('times_channel_%d_filtered_CAR.mat', ch))
            fullfile(channelDir, sprintf('times_channel_%d_filtered_CARfiltered.mat', ch))
            fullfile(masterPath, sprintf('times_channel_%d.mat', ch))
            fullfile(channelDir, sprintf('times_channel_%d.mat', ch))
        };

        wcFile = '';
        for ci = 1:numel(candidates)
            if exist(candidates{ci}, 'file')
                wcFile = candidates{ci};
                break;
            end
        end

        if isempty(wcFile)
            fprintf('  Channel %d: no wave_clus output found — skipping\n', ch);
            continue;
        end

        [~, wcName, ~] = fileparts(wcFile);
        fprintf('  Channel %d: loading %s\n', ch, wcName);

        wc = load(wcFile, 'cluster_class', 'spikes');
        cc = wc.cluster_class;
        has_waveforms = isfield(wc, 'spikes');
        cluster_ids = unique(cc(:,1));
        first_sample_s = get_first_sample_seconds(master, ch);
        cc_times_s = cc(:, 2)' / 1000;
        [cc_times_s, ref_mode] = normalize_spike_time_reference( ...
            cc_times_s, first_sample_s, all_event_times);
        fprintf('    Spike time reference: %s\n', ref_mode);

        for cid = cluster_ids(:)'
            if cid == 0, continue; end

            mask = cc(:,1) == cid;
            n_spikes = sum(mask);

            if n_spikes < opts.min_spikes
                fprintf('    Cl%d: %d spikes (< %d) — skipping\n', ...
                    cid, n_spikes, opts.min_spikes);
                continue;
            end

            spike_s_raw = cc_times_s(mask);
            [spike_s, sort_idx] = sort(spike_s_raw);

            units(end+1).label = sprintf('Ch%d Cl%d (n=%d)', ch, cid, n_spikes);
            units(end).spike_times_s = spike_s;
            units(end).channel = ch;
            units(end).cluster = cid;
            units(end).n_spikes = n_spikes;

            % Extract waveforms and compute metrics
            if has_waveforms
                unit_wf = wc.spikes(mask, :);
                unit_wf = unit_wf(sort_idx, :);
                avg_wf = mean(unit_wf, 1);
                units(end).waveforms = unit_wf;
                units(end).avg_waveform = avg_wf;

                [hw, ttp] = compute_waveform_metrics(avg_wf, Fs);
                units(end).half_width_ms = hw;
                units(end).trough_to_peak_ms = ttp;
                fprintf('    Cl%d: %d spikes | HW=%.2f ms | TTP=%.2f ms\n', ...
                    cid, n_spikes, hw, ttp);
            else
                units(end).waveforms = [];
                units(end).avg_waveform = [];
                units(end).half_width_ms = NaN;
                units(end).trough_to_peak_ms = NaN;
                fprintf('    Cl%d: %d spikes (no waveforms)\n', cid, n_spikes);
            end
        end
    end

    n_units = numel(units);
    fprintf('\nTotal units: %d\n', n_units);

    if n_units == 0
        warning('No units found. Did wave_clus produce times_*.mat files in %s?', channelDir);
        return;
    end

    %% Bin edges
    psth_edges = opts.window(1):opts.bin_size:opts.window(2);
    psth_centers = psth_edges(1:end-1) + opts.bin_size/2;

    heat_edges = opts.window(1):opts.heatmap_bin:opts.window(2);
    heat_centers = heat_edges(1:end-1) + opts.heatmap_bin/2;

    bl_mask = heat_centers >= opts.baseline_window(1) & ...
              heat_centers < opts.baseline_window(2);

    %% ================================================================
    %  FIGURE 1: Waveform classification scatter (half-width vs TTP)
    %  ================================================================
    hw_vals = [units.half_width_ms];
    ttp_vals = [units.trough_to_peak_ms];
    pyr_hw_cutoff = 0.3;   % ms
    pyr_ttp_cutoff = 0.8;  % ms
    valid_shape = ~isnan(hw_vals) & ~isnan(ttp_vals);
    is_pyramidal = valid_shape & hw_vals > pyr_hw_cutoff & ttp_vals > pyr_ttp_cutoff;
    if any(valid_shape)
        hw_axis_max = max(hw_vals(valid_shape));
        ttp_axis_max = max(ttp_vals(valid_shape));
    else
        hw_axis_max = 1;
        ttp_axis_max = 1;
    end

    if any(valid_shape)
        fig_class = figure('Position', [100 100 900 650], 'Color', 'w', 'Visible', 'off');
        hold on;

        % Classification boundaries (from your MATLAB analysis script)
        xline(pyr_hw_cutoff, '--k', 'LineWidth', 1.2, 'Alpha', 0.5);
        yline(pyr_ttp_cutoff, '--k', 'LineWidth', 1.2, 'Alpha', 0.5);

        % Plot each unit with its colour
        for ui = 1:n_units
            scatter(hw_vals(ui), ttp_vals(ui), 120, get_unit_color(ui), ...
                'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
            text(hw_vals(ui) + 0.01, ttp_vals(ui) + 0.03, units(ui).label, ...
                'FontSize', 7);
        end

        xlabel('Half-Width (ms)', 'FontSize', 12);
        ylabel('Trough-to-Peak Time (ms)', 'FontSize', 12);
        title('Waveform Classification', 'FontSize', 14, 'FontWeight', 'bold');
        xlim([0 hw_axis_max*1.3 + 0.1]);
        ylim([0 ttp_axis_max*1.3 + 0.2]);
        grid on; box on;

        % Quadrant labels
        x_lim = xlim;
        y_lim = ylim;
        text(0.72*x_lim(2), 0.82*y_lim(2), 'Putative Pyramidal', 'FontSize', 11, ...
            'HorizontalAlignment', 'center', 'Color', [0.4 0.4 0.4]);
        text(0.25*x_lim(2), 0.25*y_lim(2), 'Putative FS/Interneuron', 'FontSize', 11, ...
            'HorizontalAlignment', 'center', 'Color', [0.4 0.4 0.4]);

        % Insets: overlaid mean waveforms for each class
        pyr_idx = find(is_pyramidal);
        fs_idx = find(valid_shape & ~is_pyramidal);

        ax_inset_pyr = axes('Position', [0.60 0.62 0.33 0.26], 'Parent', fig_class); %#ok<LAXES>
        hold(ax_inset_pyr, 'on');
        n_pyr_wf = 0;
        for pi = pyr_idx
            if isempty(units(pi).avg_waveform), continue; end
            wf = units(pi).avg_waveform;
            t_ms = (0:numel(wf)-1) / Fs * 1000;
            plot(ax_inset_pyr, t_ms, wf, 'k', 'LineWidth', 1.0);
            n_pyr_wf = n_pyr_wf + 1;
        end
        if n_pyr_wf == 0
            text(ax_inset_pyr, 0.5, 0.5, 'No pyramidal waveforms', ...
                'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                'FontSize', 8, 'Color', [0.3 0.3 0.3]);
        end
        title(ax_inset_pyr, sprintf('Pyramidal means (n=%d)', n_pyr_wf), 'FontSize', 9);
        xlabel(ax_inset_pyr, 'ms', 'FontSize', 8);
        ylabel(ax_inset_pyr, 'mV', 'FontSize', 8);
        set(ax_inset_pyr, 'FontSize', 7);
        box(ax_inset_pyr, 'on');

        ax_inset_fs = axes('Position', [0.60 0.28 0.33 0.26], 'Parent', fig_class); %#ok<LAXES>
        hold(ax_inset_fs, 'on');
        n_fs_wf = 0;
        for fi = fs_idx
            if isempty(units(fi).avg_waveform), continue; end
            wf = units(fi).avg_waveform;
            t_ms = (0:numel(wf)-1) / Fs * 1000;
            plot(ax_inset_fs, t_ms, wf, 'Color', [0.85 0.1 0.1], 'LineWidth', 1.0);
            n_fs_wf = n_fs_wf + 1;
        end
        if n_fs_wf == 0
            text(ax_inset_fs, 0.5, 0.5, 'No FS waveforms', ...
                'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                'FontSize', 8, 'Color', [0.3 0.3 0.3]);
        end
        title(ax_inset_fs, sprintf('FS means (n=%d)', n_fs_wf), 'FontSize', 9);
        xlabel(ax_inset_fs, 'ms', 'FontSize', 8);
        ylabel(ax_inset_fs, 'mV', 'FontSize', 8);
        set(ax_inset_fs, 'FontSize', 7);
        box(ax_inset_fs, 'on');

        figFile = fullfile(outDir, 'waveform_classification.png');
        saveas(fig_class, figFile);
        fprintf('Saved: %s\n', figFile);
        close(fig_class);
    end

    %% ================================================================
    %  FIGURE 2: Population summary per event tag
    %  ================================================================
    for ti = 1:numel(unique_tags)
        tag = unique_tags(ti);
        event_mask = all_event_tags == tag;
        event_times = all_event_times(event_mask);
        n_trials = numel(event_times);

        if tag_info.isKey(int32(tag))
            tag_label = tag_info(int32(tag));
        else
            tag_label = sprintf('tag %d', tag);
        end

        fprintf('\n=== %s (tag %d): %d trials ===\n', tag_label, tag, n_trials);

        fig = figure('Position', [50 50 1400 300 + n_units*60], ...
            'Color', 'w', 'Visible', 'off');

        % --- RASTER (bigger dots) ---
        ax_raster = subplot(3, 1, 1);
        hold on;
        ytick_pos = [];
        ytick_labels_r = {};
        y_offset = 0;

        for ui = 1:n_units
            spk = units(ui).spike_times_s;
            for trial_i = 1:n_trials
                t0 = event_times(trial_i);
                rel = spk(spk >= t0 + opts.window(1) & spk < t0 + opts.window(2)) - t0;
                if ~isempty(rel)
                    yy = y_offset + trial_i;
                    scatter(rel, repmat(yy, numel(rel), 1), 6, ...
                        get_unit_color(ui), 'filled', 'MarkerFaceAlpha', 0.7);
                end
            end
            ytick_pos(end+1) = y_offset + n_trials/2;  %#ok
            ytick_labels_r{end+1} = units(ui).label;    %#ok
            y_offset = y_offset + n_trials + 3;
        end

        xlim(opts.window); ylim([-1 y_offset]);
        set(ax_raster, 'YTick', ytick_pos, 'YTickLabel', ytick_labels_r, ...
            'FontSize', 7, 'YDir', 'reverse');
        xlabel('Time from event onset (s)');
        title(sprintf('Spike Rasters — %s (%d trials)', tag_label, n_trials), 'FontSize', 11);

        % --- PSTH ---
        ax_psth = subplot(3, 1, 2);
        hold on;

        for ui = 1:n_units
            spk = units(ui).spike_times_s;
            trial_rates = zeros(n_trials, numel(psth_centers));
            for trial_i = 1:n_trials
                t0 = event_times(trial_i);
                rel = spk(spk >= t0 + opts.window(1) & spk < t0 + opts.window(2)) - t0;
                trial_rates(trial_i, :) = histcounts(rel, psth_edges) / opts.bin_size;
            end
            mean_rate = mean(trial_rates, 1);
            sem_rate = std(trial_rates, 0, 1) / sqrt(n_trials);
            col = get_unit_color(ui);

            fill([psth_centers fliplr(psth_centers)], ...
                [mean_rate+sem_rate fliplr(mean_rate-sem_rate)], ...
                col, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
            plot(psth_centers, mean_rate, 'Color', col, 'LineWidth', 1.2, ...
                'DisplayName', units(ui).label);
        end

        xlim(opts.window);
        xlabel('Time from event onset (s)');
        ylabel('Firing rate (Hz)');
        title('PSTH (mean +/- SEM)', 'FontSize', 10);
        legend('show', 'Location', 'eastoutside', 'FontSize', 6);
        grid on;

        % --- Z-SCORED HEATMAP (fixed scale -3 to +3) ---
        ax_heat = subplot(3, 1, 3);
        z_matrix = zeros(n_units, numel(heat_centers));
        heat_labels = {};

        for ui = 1:n_units
            z_matrix(ui, :) = compute_zscore_row(units(ui).spike_times_s, ...
                event_times, heat_edges, opts.heatmap_bin, bl_mask);
            heat_labels{ui} = units(ui).label;
        end

        % Clamp to [-3, +3]
        z_matrix = max(min(z_matrix, 3), -3);

        imagesc(heat_centers, 1:n_units, z_matrix);
        set(gca, 'CLim', [-3 3]);
        colormap(ax_heat, redblue(256));
        cb = colorbar;
        cb.Label.String = 'Z-score';
        set(gca, 'YTick', 1:n_units, 'YTickLabel', heat_labels, 'FontSize', 7);
        xlabel('Time from event onset (s)');
        title(sprintf('Z-scored firing rate (baseline: %.1f to %.1f s)', ...
            opts.baseline_window(1), opts.baseline_window(2)), 'FontSize', 10);

        sgtitle(sprintf('%s — %d units, %d trials', tag_label, n_units, n_trials), ...
            'FontSize', 13, 'FontWeight', 'bold');

        figFile = fullfile(outDir, sprintf('peri_event_%s_tag%d.png', ...
            strrep(tag_label, ' ', '_'), tag));
        saveas(fig, figFile);
        fprintf('  Saved: %s\n', figFile);
        close(fig);
    end

    %% ================================================================
    %  FIGURE 3: Per-unit figures (waveform, scatter, raster, z-score + 5V rate)
    %  ================================================================
    % Use primary_tag events for per-unit rasters/heatmaps
    primary_mask = all_event_tags == opts.primary_tag;
    primary_events = all_event_times(primary_mask);
    n_primary = numel(primary_events);

    if tag_info.isKey(int32(opts.primary_tag))
        primary_label = tag_info(int32(opts.primary_tag));
    else
        primary_label = sprintf('tag %d', opts.primary_tag);
    end

    fprintf('\n=== Per-unit figures (%s, %d trials) ===\n', primary_label, n_primary);

    unitDir = fullfile(outDir, 'per_unit');
    if ~exist(unitDir, 'dir'), mkdir(unitDir); end

    for ui = 1:n_units
        fig_u = figure('Position', [50 50 1200 900], 'Color', 'w', 'Visible', 'off');

        % --- Panel A: average waveform with std shading ---
        ax_wf = subplot(2, 2, 1);
        hold on;

        if ~isempty(units(ui).avg_waveform)
            wf = units(ui).avg_waveform;
            wf_std = std(units(ui).waveforms, 0, 1);
            n_samples = numel(wf);
            time_ms = (0:n_samples-1) / Fs * 1000;

            fill([time_ms fliplr(time_ms)], ...
                [wf+wf_std fliplr(wf-wf_std)], ...
                [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
            plot(time_ms, wf, 'k', 'LineWidth', 2);

            % Mark trough and peak
            [~, trough_idx] = min(wf);
            [~, peak_idx_rel] = max(wf(trough_idx:end));
            peak_idx = trough_idx + peak_idx_rel - 1;
            plot(time_ms(trough_idx), wf(trough_idx), 'bv', ...
                'MarkerFaceColor', 'b', 'MarkerSize', 10);
            plot(time_ms(peak_idx), wf(peak_idx), 'r^', ...
                'MarkerFaceColor', 'r', 'MarkerSize', 10);

            xlabel('Time (ms)');
            ylabel('Amplitude (mV)');
        end
        title(sprintf('%s\nHW=%.2f ms | TTP=%.2f ms', ...
            units(ui).label, units(ui).half_width_ms, units(ui).trough_to_peak_ms), ...
            'FontSize', 10);
        grid on; box on;

        % --- Panel B: position on classification scatter ---
        ax_sc = subplot(2, 2, 2);
        hold on;

        % Plot all units in grey
        for uj = 1:n_units
            if uj == ui, continue; end
            scatter(hw_vals(uj), ttp_vals(uj), 60, [0.75 0.75 0.75], ...
                'filled', 'MarkerEdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.3);
        end

        % Highlight current unit
        scatter(hw_vals(ui), ttp_vals(ui), 180, get_unit_color(ui), ...
            'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

        % Classification lines
        xline(pyr_hw_cutoff, '--k', 'LineWidth', 1, 'Alpha', 0.4);
        yline(pyr_ttp_cutoff, '--k', 'LineWidth', 1, 'Alpha', 0.4);

        % Classify
        if hw_vals(ui) > pyr_hw_cutoff && ttp_vals(ui) > pyr_ttp_cutoff
            cell_type = 'Putative Pyramidal';
        else
            cell_type = 'Putative FS/Interneuron';
        end

        xlabel('Half-Width (ms)');
        ylabel('Trough-to-Peak Time (ms)');
        title(cell_type, 'FontSize', 10);
        xlim([0 hw_axis_max*1.3 + 0.1]);
        ylim([0 ttp_axis_max*1.3 + 0.2]);
        grid on; box on;

        % --- Panel C: main raster + narrow onset zoom raster ---
        ax_rr_slot = subplot(2, 2, 3);
        rr_pos = get(ax_rr_slot, 'Position');
        delete(ax_rr_slot);

        main_w_frac = 0.74;
        gap_w_frac = 0.04;
        zoom_w_frac = 0.22;
        main_pos = [rr_pos(1), rr_pos(2), rr_pos(3)*main_w_frac, rr_pos(4)];
        zoom_pos = [rr_pos(1) + rr_pos(3)*(main_w_frac + gap_w_frac), ...
                    rr_pos(2), rr_pos(3)*zoom_w_frac, rr_pos(4)];

        ax_rr = axes('Position', main_pos);
        hold(ax_rr, 'on');
        ax_rr_zoom = axes('Position', zoom_pos);
        hold(ax_rr_zoom, 'on');

        zoom_window = [-0.01 0.05];  % -10 to +50 ms
        spk = units(ui).spike_times_s;
        for trial_i = 1:n_primary
            t0 = primary_events(trial_i);
            rel_main = spk(spk >= t0 + opts.window(1) & spk < t0 + opts.window(2)) - t0;
            rel_zoom = spk(spk >= t0 + zoom_window(1) & spk < t0 + zoom_window(2)) - t0;

            if ~isempty(rel_main)
                scatter(ax_rr, rel_main, repmat(trial_i, numel(rel_main), 1), 10, ...
                    get_unit_color(ui), 'filled', 'MarkerFaceAlpha', 0.7);
            end
            if ~isempty(rel_zoom)
                scatter(ax_rr_zoom, rel_zoom, repmat(trial_i, numel(rel_zoom), 1), 10, ...
                    get_unit_color(ui), 'filled', 'MarkerFaceAlpha', 0.7);
            end
        end

        xlim(ax_rr, opts.window);
        ylim(ax_rr, [0 n_primary + 1]);
        set(ax_rr, 'YDir', 'reverse');
        xlabel(ax_rr, 'Time from event onset (s)');
        ylabel(ax_rr, 'Trial');
        title(ax_rr, sprintf('Raster — %s (%d trials)', primary_label, n_primary), ...
            'FontSize', 10);

        xlim(ax_rr_zoom, zoom_window);
        ylim(ax_rr_zoom, [0 n_primary + 1]);
        set(ax_rr_zoom, 'YDir', 'reverse', 'YTick', []);
        xlabel(ax_rr_zoom, 'Time (s)');
        title(ax_rr_zoom, '-10 to +50 ms', 'FontSize', 9);

        % --- Panel D/E: split lower-right quadrant (z-score + 5V rate) ---
        ax_slot = subplot(2, 2, 4);
        slot_pos = get(ax_slot, 'Position');
        delete(ax_slot);

        % Z-score panel takes half the original vertical space.
        z_pos = [slot_pos(1), slot_pos(2) + 0.53*slot_pos(4), ...
                 slot_pos(3), 0.47*slot_pos(4)];
        fr_pos = [slot_pos(1), slot_pos(2), slot_pos(3), 0.43*slot_pos(4)];

        ax_hh = axes('Position', z_pos);
        hold(ax_hh, 'on');

        tag_z = [];
        tag_ylabels = {};
        for ti = 1:numel(unique_tags)
            tag = unique_tags(ti);
            evt = all_event_times(all_event_tags == tag);
            n_t = numel(evt);

            z_row = compute_zscore_row(spk, evt, heat_edges, opts.heatmap_bin, bl_mask);
            tag_z(end+1, :) = max(min(z_row, 3), -3);  %#ok

            if tag_info.isKey(int32(tag))
                tl = tag_info(int32(tag));
            else
                tl = sprintf('tag %d', tag);
            end
            tag_ylabels{end+1} = sprintf('%s (n=%d)', tl, n_t);  %#ok
        end

        imagesc(ax_hh, heat_centers, 1:numel(unique_tags), tag_z);
        set(ax_hh, 'CLim', [-3 3]);
        colormap(ax_hh, redblue(256));
        cb = colorbar(ax_hh);
        cb.Label.String = 'Z-score';
        cb.FontSize = 7;
        set(ax_hh, 'YTick', 1:numel(unique_tags), 'YTickLabel', tag_ylabels, ...
            'XTickLabel', [], 'FontSize', 7);
        ylabel(ax_hh, 'Condition');
        title(ax_hh, 'Z-scored rate by condition', 'FontSize', 9);

        % 5V-only firing rate panel (absolute Hz, unit-specific y-axis max)
        ax_fr = axes('Position', fr_pos);
        hold(ax_fr, 'on');

        if isempty(fiveV_events)
            text(ax_fr, mean(opts.window), 0.5, 'No 5V trials found', ...
                'HorizontalAlignment', 'center', 'FontSize', 8, ...
                'Color', [0.3 0.3 0.3]);
            xlim(ax_fr, opts.window);
            ylim(ax_fr, [0 1]);
        else
            n_fiveV = numel(fiveV_events);
            trial_rates_5v = zeros(n_fiveV, numel(psth_centers));

            for trial_i = 1:n_fiveV
                t0 = fiveV_events(trial_i);
                rel = spk(spk >= t0 + psth_edges(1) & spk < t0 + psth_edges(end)) - t0;
                trial_rates_5v(trial_i, :) = histcounts(rel, psth_edges) / opts.bin_size;
            end

            mean_rate_5v = mean(trial_rates_5v, 1);
            unit_max_rate = max(trial_rates_5v(:));

            plot(ax_fr, psth_centers, mean_rate_5v, 'Color', get_unit_color(ui), ...
                'LineWidth', 1.6);
            xlim(ax_fr, opts.window);
            if unit_max_rate < 1e-6
                ylim(ax_fr, [0 1]);
            else
                ylim(ax_fr, [0 unit_max_rate]);
            end
        end

        xlabel(ax_fr, 'Time from event onset (s)');
        ylabel(ax_fr, 'Firing rate (Hz)');
        title(ax_fr, sprintf('5V firing rate (mean; unit peak y-max, n=%d)', numel(fiveV_events)), ...
            'FontSize', 9);
        grid(ax_fr, 'on');
        box(ax_fr, 'on');

        sgtitle(sprintf('%s — %s', units(ui).label, cell_type), ...
            'FontSize', 13, 'FontWeight', 'bold');

        figFile = fullfile(unitDir, sprintf('unit_Ch%d_Cl%d.png', ...
            units(ui).channel, units(ui).cluster));
        saveas(fig_u, figFile);
        fprintf('  Saved: %s\n', figFile);
        close(fig_u);
    end

    %% ================================================================
    %  FIGURE 4: Per-unit waveform overlay around LED onset
    %  ================================================================
    waveDir = fullfile(outDir, 'waveform_led_comparison');
    if ~exist(waveDir, 'dir'), mkdir(waveDir); end

    fprintf('\n=== LED waveform comparisons (%d LED trials) ===\n', numel(led_events));

    for ui = 1:n_units
        fig_w = figure('Position', [80 80 760 430], 'Color', 'w', 'Visible', 'off');
        ax_w = axes('Parent', fig_w);
        hold(ax_w, 'on');

        if ~isempty(units(ui).waveforms)
            spk = units(ui).spike_times_s;
            pre_mask = spike_mask_in_windows(spk, led_events, [-1 0]);
            early_mask = spike_mask_in_windows(spk, led_events, [0 0.1]);

            n_pre = sum(pre_mask);
            n_early = sum(early_mask);

            n_samples = size(units(ui).waveforms, 2);
            wf_time_ms = (0:n_samples-1) / Fs * 1000;

            legend_h = gobjects(0);
            legend_txt = {};

            if n_pre > 0
                pre_avg = mean(units(ui).waveforms(pre_mask, :), 1);
                h1 = plot(ax_w, wf_time_ms, pre_avg, 'Color', [0.2 0.45 0.9], ...
                    'LineWidth', 2);
                legend_h(end+1) = h1; %#ok<AGROW>
                legend_txt{end+1} = sprintf('Pre-LED (-1 to 0 s), n=%d', n_pre); %#ok<AGROW>
            end

            if n_early > 0
                early_avg = mean(units(ui).waveforms(early_mask, :), 1);
                h2 = plot(ax_w, wf_time_ms, early_avg, 'Color', [0.9 0.25 0.2], ...
                    'LineWidth', 2);
                legend_h(end+1) = h2; %#ok<AGROW>
                legend_txt{end+1} = sprintf('Early LED (0 to 0.1 s), n=%d', n_early); %#ok<AGROW>
            end

            if isempty(legend_h)
                text(ax_w, mean(wf_time_ms), 0, 'No spikes in requested windows', ...
                    'HorizontalAlignment', 'center', 'Color', [0.35 0.35 0.35], ...
                    'FontSize', 9);
            else
                legend(ax_w, legend_h, legend_txt, 'Location', 'best', 'FontSize', 9);
            end

            xlabel(ax_w, 'Time (ms)');
            ylabel(ax_w, 'Amplitude (mV)');
            grid(ax_w, 'on');
            box(ax_w, 'on');
        else
            text(ax_w, 0.5, 0.5, 'No waveforms available for this unit', ...
                'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                'Color', [0.35 0.35 0.35], 'FontSize', 10);
            axis(ax_w, 'off');
        end

        title(ax_w, sprintf('%s\nWaveforms around LED onset', units(ui).label), 'FontSize', 11);

        figFile = fullfile(waveDir, sprintf('unit_waveforms_Ch%d_Cl%d.png', ...
            units(ui).channel, units(ui).cluster));
        saveas(fig_w, figFile);
        fprintf('  Saved: %s\n', figFile);
        close(fig_w);
    end

    fprintf('\nAll figures saved to: %s\n', outDir);
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

%% ================================================================
%  HELPER FUNCTIONS
%  ================================================================

function z_row = compute_zscore_row(spike_s, event_times, heat_edges, bin_width, bl_mask)
    % Compute z-scored firing rate for one unit / one condition
    n_trials = numel(event_times);
    n_bins = numel(heat_edges) - 1;
    trial_rates = zeros(n_trials, n_bins);

    for trial_i = 1:n_trials
        t0 = event_times(trial_i);
        rel = spike_s(spike_s >= t0 + heat_edges(1) & spike_s < t0 + heat_edges(end)) - t0;
        trial_rates(trial_i, :) = histcounts(rel, heat_edges) / bin_width;
    end

    bl_vals = trial_rates(:, bl_mask);
    mu = mean(bl_vals(:));
    sd = std(bl_vals(:));

    if sd < 1e-6
        z_row = zeros(1, n_bins);
    else
        z_row = (mean(trial_rates, 1) - mu) / sd;
    end
end

function mask = spike_mask_in_windows(spike_s, event_times, win_s)
    % Logical mask of spikes that fall within any event-aligned time window.
    mask = false(size(spike_s));
    if isempty(spike_s) || isempty(event_times)
        return;
    end

    for ei = 1:numel(event_times)
        t0 = event_times(ei);
        mask = mask | (spike_s >= t0 + win_s(1) & spike_s < t0 + win_s(2));
    end
end

function [half_width_ms, trough_to_peak_ms] = compute_waveform_metrics(avg_wf, Fs)
    % Half-width: duration at half-minimum amplitude
    % Trough-to-peak: time from trough (min) to subsequent peak (max)

    [~, trough_idx] = min(avg_wf);
    [~, peak_rel] = max(avg_wf(trough_idx:end));
    peak_idx = trough_idx + peak_rel - 1;

    trough_to_peak_ms = (peak_idx - trough_idx) / Fs * 1000;

    % Half-amplitude crossing for half-width
    half_amp = avg_wf(trough_idx) / 2;

    % Left crossing: last point before trough that's above half-amp
    left_cross = find(avg_wf(1:trough_idx) >= half_amp, 1, 'last');
    if isempty(left_cross), left_cross = 1; end

    % Right crossing: first point after trough that's above half-amp
    right_pts = find(avg_wf(trough_idx:end) >= half_amp, 1, 'first');
    if isempty(right_pts)
        right_cross = numel(avg_wf);
    else
        right_cross = trough_idx + right_pts - 1;
    end

    half_width_ms = (right_cross - left_cross) / Fs * 1000;
end

function col = get_unit_color(idx)
    colors = [
        0.12 0.47 0.71
        1.00 0.50 0.05
        0.17 0.63 0.17
        0.84 0.15 0.16
        0.58 0.40 0.74
        0.55 0.34 0.29
        0.89 0.47 0.76
        0.50 0.50 0.50
        0.74 0.74 0.13
        0.09 0.75 0.81
    ];
    col = colors(mod(idx-1, size(colors,1)) + 1, :);
end

function cmap = redblue(n)
    if nargin < 1, n = 256; end
    half = floor(n/2);
    r1 = linspace(0.23, 1, half)';
    g1 = linspace(0.30, 1, half)';
    b1 = linspace(0.75, 1, half)';
    r2 = linspace(1, 0.75, n-half)';
    g2 = linspace(1, 0.15, n-half)';
    b2 = linspace(1, 0.15, n-half)';
    cmap = [r1 g1 b1; r2 g2 b2];
end
