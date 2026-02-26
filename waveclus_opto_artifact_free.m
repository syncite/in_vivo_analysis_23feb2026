function waveclus_opto_artifact_free(master_mat_file, varargin)
% WAVECLUS_OPTO_ARTIFACT_FREE  Artifact-aware wave_clus workflow for LED experiments.
%
% This wrapper keeps wave_clus clustering templates clean by excluding spikes
% from LED windows during clustering, then optionally assigns LED-window spikes
% back to curated units using template distance.
% Input channel files (channel_*_filtered_CAR.mat) are not modified.
%
% Two modes:
%   1) 'prepare' (default)
%      - Detect spikes on full channel trace (Get_spikes)
%      - Split spikes into clean vs dirty using LED windows
%      - Save explicit *_clean_spikes.mat and *_dirty_spikes.mat files
%      - Run Do_clustering on clean spikes
%
%   2) 'assign'
%      - Load manually curated clean-only times_*.mat
%      - Assign non-template spikes using wave_clus-style force membership
%        (template_type / template_sdnum / template_k / template_k_min)
%      - Write combined output and (optionally) overwrite canonical times file
%
% Usage:
%   waveclus_opto_artifact_free('path/to/master.mat', 'mode', 'prepare');
%   % ... manual curation in wave_clus GUI ...
%   waveclus_opto_artifact_free('path/to/master.mat', 'mode', 'assign');
%
% Key options:
%   'mode'                  : 'prepare' | 'assign'            [prepare]
%   'channels'              : channel list (empty = all)      []
%   'par'                   : wave_clus parameter struct       []
%   'led_tags'              : event tags considered LED trials [32385]
%   'led_duration_s'        : scalar or per-tag durations (s)   1.0
%   'post_led_buffer_s'     : fixed post-offset exclusion (s)   0.1
%   'time_reference_mode'   : optional override for time mode   [index_ms+firstAD]
%   'assign_unlabeled_clean': optional override                 [true]
%   'clean_distance_threshold_sd': optional override            [2.5]
%   'distance_threshold_sd' : optional override                 [2.5]
%   'save_assignment_png'   : optional override                 [true]
%   'overwrite_main_times'  : in assign mode, replace times_*.mat [true]
%   'run_clustering'        : in prepare mode, run Do_clustering [true]
%   'debug'                 : print Get_spikes path diagnostics [false]
%
% Hard-coded defaults (to keep commands simple):
%   - time reference mode: index_ms+firstAD
%   - assign unlabeled clean spikes in step2: true
%   - clean reassignment threshold: 2.5 SD
%   - dirty reassignment threshold: 2.5 SD
%   - save step2 PNG summary: true
%
% Files created per channel:
%   channel_<N>_..._spikes.mat       - full spikes (from Get_spikes)
%   channel_<N>_..._clean_spikes.mat - clean-only spikes
%   channel_<N>_..._dirty_spikes.mat - dirty-only spikes
%   channel_<N>_..._led_split.mat    - masks + metadata
%   times_..._with_led_assignments.mat (assign mode)
%   times_..._clean_only.mat (backup of curated clean file if overwritten)

    p = inputParser;
    addOptional(p, 'master_mat_file', '', @ischar);
    addParameter(p, 'mode', 'prepare', @ischar);
    addParameter(p, 'channels', [], @isnumeric);
    addParameter(p, 'par', [], @(x) isempty(x) || isstruct(x));
    addParameter(p, 'led_tags', 32385, @isnumeric);
    addParameter(p, 'led_duration_s', 1.0, @(x) isnumeric(x) && all(x > 0));
    addParameter(p, 'post_led_buffer_s', 0.1, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'time_reference_mode', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'assign_unlabeled_clean', [], @(x) isempty(x) || islogical(x));
    addParameter(p, 'clean_distance_threshold_sd', [], ...
        @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    addParameter(p, 'distance_threshold_sd', [], ...
        @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    addParameter(p, 'save_assignment_png', [], @(x) isempty(x) || islogical(x));
    addParameter(p, 'overwrite_main_times', true, @islogical);
    addParameter(p, 'run_clustering', true, @islogical);
    addParameter(p, 'debug', false, @islogical);
    parse(p, master_mat_file, varargin{:});
    opts = p.Results;

    opts.mode = lower(strtrim(opts.mode));
    if ~ismember(opts.mode, {'prepare', 'assign'})
        error('mode must be ''prepare'' or ''assign''.');
    end

    % Defaults for reproducible two-pass workflow (can be overridden).
    if isempty(opts.time_reference_mode)
        opts.time_reference_mode = 'index_ms+firstAD';
    end
    opts.time_reference_mode = normalize_time_reference_mode(opts.time_reference_mode);
    if isempty(opts.assign_unlabeled_clean)
        opts.assign_unlabeled_clean = true;
    end
    if isempty(opts.distance_threshold_sd)
        opts.distance_threshold_sd = 2.5;
    end
    if isempty(opts.clean_distance_threshold_sd)
        opts.clean_distance_threshold_sd = opts.distance_threshold_sd;
    end
    if isempty(opts.save_assignment_png)
        opts.save_assignment_png = true;
    end

    if isscalar(opts.led_duration_s)
        opts.led_duration_s = repmat(opts.led_duration_s, size(opts.led_tags));
    end
    if numel(opts.led_duration_s) ~= numel(opts.led_tags)
        error('led_duration_s must be scalar or one duration per led_tags entry.');
    end

    if isempty(opts.par) && strcmp(opts.mode, 'prepare')
        if exist('set_parameters_custom', 'file') == 2
            opts.par = set_parameters_custom();
        elseif exist('set_parameters', 'file') == 2
            opts.par = set_parameters();
        else
            error(['No parameter struct provided and no set_parameters_custom()/' ...
                'set_parameters() found on path.']);
        end
    end

    if isempty(opts.master_mat_file)
        [f, d] = uigetfile('*.mat', 'Select master .mat file (from extract_and_filter)');
        if isequal(f, 0), return; end
        opts.master_mat_file = fullfile(d, f);
    end

    [master, masterPath, masterName] = load_master(opts.master_mat_file);
    channelDir = resolve_channel_dir(masterPath, masterName);

    if isempty(opts.channels)
        channels = 1:master.nChannels;
    else
        channels = opts.channels(:)';
    end

    led_windows = build_led_windows(master.Eventstime(:), master.EventTag(:), ...
        opts.led_tags(:), opts.led_duration_s, opts.post_led_buffer_s);

    if isempty(led_windows)
        warning('No LED windows built from the requested tags. All spikes will be treated as clean.');
    end

    fprintf('\n=== waveclus_opto_artifact_free (%s mode) ===\n', opts.mode);
    fprintf('Master: %s\n', opts.master_mat_file);
    fprintf('Channel dir: %s\n', channelDir);
    fprintf('Channels: %s\n', mat2str(channels));
    fprintf('Time reference mode: %s\n', opts.time_reference_mode);
    fprintf('Step2 assignment: clean unlabeled=%d (thr=%.2f SD), dirty thr=%.2f SD\n', ...
        opts.assign_unlabeled_clean, opts.clean_distance_threshold_sd, opts.distance_threshold_sd);
    if numel(unique(opts.led_duration_s(:))) == 1
        fprintf(['Masking rule: [LED onset, LED onset + %.3f + %.3f] s ' ...
            '(duration + fixed post-offset buffer)\n'], opts.led_duration_s(1), opts.post_led_buffer_s);
    else
        fprintf(['Masking rule: [LED onset, LED onset + duration(tag) + %.3f] s ' ...
            '(fixed post-offset buffer)\n'], opts.post_led_buffer_s);
    end
    fprintf('LED windows: %d merged intervals\n', size(led_windows, 1));

    switch opts.mode
        case 'prepare'
            run_prepare(master, channelDir, channels, led_windows, opts);
        case 'assign'
            run_assign(master, channelDir, channels, opts);
    end

    fprintf('\nDone.\n');
end

function run_prepare(master, channelDir, channels, led_windows, opts)
    if exist('Get_spikes', 'file') ~= 2
        error('Get_spikes not found on MATLAB path. Add wave_clus to path first.');
    end
    if opts.run_clustering && exist('Do_clustering', 'file') ~= 2
        error('Do_clustering not found on MATLAB path. Add wave_clus to path first.');
    end

    sr = get_sampling_rate(master, opts.par);

    for ch = channels
        chFile = resolve_channel_file(channelDir, ch);
        if isempty(chFile)
            warning('Channel %d: no input channel file found, skipping.', ch);
            continue;
        end

        fprintf('\n[prepare] Channel %d\n', ch);
        fprintf('  Input: %s\n', chFile);

        ensure_sr_in_mat_file(chFile, sr, opts.debug);

        if opts.debug
            chDir = fileparts(chFile);
            fprintf('  [debug] MATLAB pwd before Get_spikes: %s\n', pwd);
            fprintf('  [debug] Channel file dir: %s\n', chDir);
            print_spike_file_listing(chDir, 'before');
            if ~strcmpi(chDir, pwd)
                print_spike_file_listing(pwd, 'before');
            end
        end

        run_get_spikes(chFile, opts.par, opts.debug);

        if opts.debug
            chDir = fileparts(chFile);
            fprintf('  [debug] MATLAB pwd after Get_spikes: %s\n', pwd);
            print_spike_file_listing(chDir, 'after');
            if ~strcmpi(chDir, pwd)
                print_spike_file_listing(pwd, 'after');
            end
        end
        spkFileFull = resolve_spike_file_after_get_spikes(chFile);

        fullData = load(spkFileFull);
        if ~isfield(fullData, 'spikes')
            error('Spike file missing ''spikes'': %s', spkFileFull);
        end
        if ~isfield(fullData, 'sr')
            ensure_sr_in_mat_file(spkFileFull, sr, opts.debug);
            fullData.sr = sr;
        end
        [index_raw, index_field] = get_index_field(fullData, spkFileFull);

        first_sample_s = get_first_sample_seconds(master, ch);
        [spike_s, ref_mode] = infer_spike_times_seconds( ...
            index_raw, sr, first_sample_s, master.Eventstime(:), opts.time_reference_mode);

        dirty_mask = in_any_window(spike_s, led_windows);
        clean_mask = ~dirty_mask;

        n_total = numel(index_raw);
        n_clean = sum(clean_mask);
        n_dirty = sum(dirty_mask);
        fprintf('  Spikes: total=%d | clean=%d | dirty=%d | time ref=%s\n', ...
            n_total, n_clean, n_dirty, ref_mode);

        if n_clean < 20
            warning('Channel %d has only %d clean spikes; clustering may fail.', ch, n_clean);
        end

        spkFileClean = strrep(spkFileFull, '_spikes.mat', '_clean_spikes.mat');
        spkFileDirty = strrep(spkFileFull, '_spikes.mat', '_dirty_spikes.mat');
        splitFile = strrep(spkFileFull, '_spikes.mat', '_led_split.mat');

        cleanData = subset_spike_struct(fullData, clean_mask, n_total);
        dirtyData = subset_spike_struct(fullData, dirty_mask, n_total);
        cleanData.sr = sr;
        dirtyData.sr = sr;

        save(spkFileClean, '-struct', 'cleanData', '-v7.3');
        save(spkFileDirty, '-struct', 'dirtyData', '-v7.3');

        split = struct();
        split.channel = ch;
        split.created_at = datestr(now, 'yyyy-mm-dd HH:MM:SS');
        split.channel_file = chFile;
        split.spike_file_full = spkFileFull;
        split.spike_file_all = spkFileFull; % backwards-compatible alias
        split.spike_file_clean = spkFileClean;
        split.spike_file_dirty = spkFileDirty;
        split.times_file_clean = times_file_from_spike_file(spkFileClean);
        split.times_file_main = times_file_from_spike_file(spkFileFull);
        split.times_file = split.times_file_main; % backwards-compatible alias
        split.index_field = index_field;
        split.sr = sr;
        split.first_sample_s = first_sample_s;
        split.spike_time_reference = ref_mode;
        split.clean_mask = clean_mask(:);
        split.dirty_mask = dirty_mask(:);
        split.led_windows_s = led_windows;
        split.led_duration_s = opts.led_duration_s;
        split.post_led_buffer_s = opts.post_led_buffer_s;
        split.led_tags = opts.led_tags(:);
        split.n_total = n_total;
        split.n_clean = n_clean;
        split.n_dirty = n_dirty;

        save(splitFile, '-struct', 'split');
        fprintf('  Saved clean spikes: %s\n', spkFileClean);
        fprintf('  Saved dirty spikes: %s\n', spkFileDirty);
        fprintf('  Saved split metadata: %s\n', splitFile);

        if opts.run_clustering
            run_do_clustering(spkFileClean, opts.par, opts.debug);
            tFileCleanExpected = times_file_from_spike_file(spkFileClean);
            tFileClean = resolve_times_file_after_do_clustering(spkFileClean);
            if ~isempty(tFileClean) && ~strcmpi(tFileClean, tFileCleanExpected)
                copyfile(tFileClean, tFileCleanExpected);
                if opts.debug
                    fprintf('  [debug] Non-canonical times file copied to expected name:\n');
                    fprintf('    source: %s\n', tFileClean);
                    fprintf('    target: %s\n', tFileCleanExpected);
                end
                tFileClean = tFileCleanExpected;
            elseif isempty(tFileClean)
                tFileClean = '';
            end

            tFileMain = times_file_from_spike_file(spkFileFull);
            if exist(tFileClean, 'file')
                copyfile(tFileClean, tFileMain);
                fprintf('  Clean clustering file: %s\n', tFileClean);
                fprintf('  Main times file for GUI: %s\n', tFileMain);
            else
                warning('Do_clustering finished but expected times file not found: %s', tFileCleanExpected);
            end
        end
    end
end

function run_assign(master, channelDir, channels, opts)
    for ch = channels
        chFile = resolve_channel_file(channelDir, ch);
        if isempty(chFile)
            warning('Channel %d: no input channel file found, skipping.', ch);
            continue;
        end

        spkFile = spike_file_from_channel_file(chFile);
        splitFile = strrep(spkFile, '_spikes.mat', '_led_split.mat');
        if ~exist(splitFile, 'file')
            warning('Channel %d: split file not found (%s). Run prepare mode first.', ch, splitFile);
            continue;
        end

        split = load(splitFile);
        if isfield(split, 'spike_file_full')
            spkFileFull = split.spike_file_full;
        elseif isfield(split, 'spike_file_all')
            spkFileFull = split.spike_file_all;
        else
            spkFileFull = '';
        end

        if isfield(split, 'spike_file_clean')
            spkFileClean = split.spike_file_clean;
        else
            spkFileClean = strrep(spkFile, '_spikes.mat', '_clean_spikes.mat');
            if ~exist(spkFileClean, 'file')
                legacy = strrep(spkFile, '_spikes.mat', '_spikes_clean.mat');
                if exist(legacy, 'file')
                    spkFileClean = legacy;
                end
            end
        end

        if isfield(split, 'spike_file_dirty')
            spkFileDirty = split.spike_file_dirty;
        else
            spkFileDirty = strrep(spkFile, '_spikes.mat', '_dirty_spikes.mat');
            if ~exist(spkFileDirty, 'file')
                legacy = strrep(spkFile, '_spikes.mat', '_spikes_dirty.mat');
                if exist(legacy, 'file')
                    spkFileDirty = legacy;
                end
            end
        end

        if isempty(spkFileFull) || ~exist(spkFileFull, 'file')
            warning('Channel %d: full spike file missing (%s).', ch, spkFileFull);
            continue;
        end
        if ~exist(spkFileClean, 'file')
            warning('Channel %d: clean spike file missing (%s).', ch, spkFileClean);
            continue;
        end
        if ~exist(spkFileDirty, 'file')
            warning('Channel %d: dirty spike file missing (%s).', ch, spkFileDirty);
            continue;
        end

        if isfield(split, 'times_file_main') && ~isempty(split.times_file_main)
            timesFile = split.times_file_main;
        elseif isfield(split, 'times_file') && ~isempty(split.times_file)
            timesFile = split.times_file;
        else
            timesFile = times_file_from_spike_file(spkFile);
        end
        if ~exist(timesFile, 'file')
            warning('Channel %d: curated times file not found (%s).', ch, timesFile);
            continue;
        end

        fullData = load(spkFileFull);
        if ~isfield(fullData, 'spikes')
            warning('Channel %d: full spike file has no spikes, skipping.', ch);
            continue;
        end
        [full_index_raw, ~] = get_index_field(fullData, spkFileFull);

        cleanData = load(spkFileClean);
        dirtyData = load(spkFileDirty);
        if ~isfield(cleanData, 'spikes') || ~isfield(dirtyData, 'spikes')
            warning('Channel %d: clean/dirty spike files missing spikes field.', ch);
            continue;
        end
        [clean_index, ~] = get_index_field(cleanData, spkFileClean);
        [dirty_index, ~] = get_index_field(dirtyData, spkFileDirty);

        if numel(clean_index) + numel(dirty_index) ~= numel(full_index_raw)
            warning(['Channel %d: clean+dirty spike counts (%d) do not equal full spike count (%d). ' ...
                'Proceeding with clean/dirty files only.'], ...
                ch, numel(clean_index) + numel(dirty_index), numel(full_index_raw));
        end

        if isfield(split, 'clean_mask') && isfield(split, 'dirty_mask')
            if numel(full_index_raw) == numel(split.clean_mask) && ...
                    numel(full_index_raw) == numel(split.dirty_mask)
                clean_mask = logical(split.clean_mask(:));
                dirty_mask = logical(split.dirty_mask(:));
            else
                clean_mask = false(numel(full_index_raw), 1);
                dirty_mask = false(numel(full_index_raw), 1);
            end
        else
            clean_mask = false(numel(full_index_raw), 1);
            dirty_mask = false(numel(full_index_raw), 1);
        end

        if ~any(clean_mask) && ~any(dirty_mask)
            [~, loc_clean] = ismember(round(clean_index * 1000), round(full_index_raw * 1000));
            [~, loc_dirty] = ismember(round(dirty_index * 1000), round(full_index_raw * 1000));
            clean_mask(loc_clean(loc_clean > 0)) = true;
            dirty_mask(loc_dirty(loc_dirty > 0)) = true;
        end
        if any(clean_mask & dirty_mask)
            overlap = clean_mask & dirty_mask;
            dirty_mask(overlap) = false;
        end
        if ~all(clean_mask | dirty_mask)
            unresolved = ~(clean_mask | dirty_mask);
            clean_mask(unresolved) = false;
            dirty_mask(unresolved) = false;
        end

        if sum(clean_mask) ~= numel(clean_index) || sum(dirty_mask) ~= numel(dirty_index)
            warning(['Channel %d: mask reconstruction mismatch (mask clean=%d dirty=%d, ' ...
                'files clean=%d dirty=%d). Using file-derived spikes/times for assignment.'], ...
                ch, sum(clean_mask), sum(dirty_mask), numel(clean_index), numel(dirty_index));
        end

        if size(cleanData.spikes,1) ~= numel(clean_index) || ...
                size(dirtyData.spikes,1) ~= numel(dirty_index)
            warning('Channel %d: index/spike row mismatch in clean or dirty file.', ch);
            continue;
        end

        timesData = load(timesFile);
        if ~isfield(timesData, 'cluster_class')
            warning('Channel %d: %s has no cluster_class, skipping.', ch, timesFile);
            continue;
        end

        cc = double(timesData.cluster_class);
        if isempty(cc)
            warning('Channel %d: empty cluster_class in %s, skipping.', ch, timesFile);
            continue;
        end

        clean_spikes = double(cleanData.spikes);
        dirty_spikes = double(dirtyData.spikes);
        clean_index = double(clean_index(:));
        dirty_index = double(dirty_index(:));

        if isfield(split, 'sr')
            sr = split.sr;
        else
            sr = [];
        end
        if isempty(sr) || ~isfinite(sr) || sr <= 0
            sr = get_sampling_rate(master, []);
        end

        if isfield(split, 'first_sample_s')
            first_sample_s = split.first_sample_s;
        else
            first_sample_s = get_first_sample_seconds(master, ch);
        end
        [cc_to_clean, mode_id, mode_label, n_match] = map_cc_to_clean_indices( ...
            clean_index, cc(:,2), sr, first_sample_s);

        if n_match < round(0.8 * size(cc,1))
            warning(['Channel %d: low cc-to-clean matching (%d/%d) using %s. ' ...
                'Check time reference consistency.'], ch, n_match, size(cc,1), mode_label);
        end

        cc_original = cc;
        clean_labels = zeros(numel(clean_index), 1);
        for i = 1:size(cc,1)
            idx = cc_to_clean(i);
            if idx > 0 && idx <= numel(clean_labels)
                clean_labels(idx) = cc(i,1);
            end
        end

        n_clean_unassigned_before = sum(clean_labels == 0);
        clean_labels_for_templates = clean_labels;
        clean_reassigned_labels = zeros(size(clean_labels));
        clean_reassigned_best_dist = inf(size(clean_labels));
        clean_unassigned_mask = false(size(clean_labels));
        if opts.assign_unlabeled_clean
            clean_unassigned_mask = clean_labels == 0;
        end

        % Assign all non-template spikes equally in one pass:
        %   targets = [clean cluster-0 (optional); dirty spikes]
        % using the same curated clean templates and one threshold.
        n_clean_targets = sum(clean_unassigned_mask);
        target_spikes = [clean_spikes(clean_unassigned_mask, :); dirty_spikes];
        unified_thr = opts.distance_threshold_sd;
        if opts.assign_unlabeled_clean && abs(opts.clean_distance_threshold_sd - opts.distance_threshold_sd) > 1e-9
            warning(['Channel %d: using unified threshold %.2f SD for all non-template spikes ' ...
                '(clean/dirty), ignoring separate clean threshold %.2f SD.'], ...
                ch, unified_thr, opts.clean_distance_threshold_sd);
        end
        force_par = resolve_force_params(timesData, unified_thr);

        [target_labels_all, target_dist_all, template_stats, force_used] = ...
            assign_spikes_waveclus_force(clean_spikes, clean_labels, target_spikes, force_par);

        if n_clean_targets > 0
            idx_unassigned = find(clean_unassigned_mask);
            clean_reassigned_labels(idx_unassigned) = target_labels_all(1:n_clean_targets);
            clean_reassigned_best_dist(idx_unassigned) = target_dist_all(1:n_clean_targets);
            clean_labels_for_templates(idx_unassigned) = target_labels_all(1:n_clean_targets);
        end

        dirty_labels = target_labels_all((n_clean_targets + 1):end);
        dirty_best_dist = target_dist_all((n_clean_targets + 1):end);

        cc_clean_reassigned = cc;
        n_cc_clean_reassigned = 0;
        for i = 1:size(cc_clean_reassigned, 1)
            idx = cc_to_clean(i);
            if idx > 0 && cc_clean_reassigned(i,1) == 0 && clean_labels_for_templates(idx) > 0
                cc_clean_reassigned(i,1) = clean_labels_for_templates(idx);
                n_cc_clean_reassigned = n_cc_clean_reassigned + 1;
            end
        end
        n_clean_unassigned_after = sum(clean_labels_for_templates == 0);

        dirty_time_ms = transform_index_to_ms(double(dirty_index), sr, first_sample_s, mode_id);
        dirty_cc = [double(dirty_labels(:)) double(dirty_time_ms(:))];

        cc_combined = [cc_clean_reassigned; dirty_cc];
        [~, sort_idx] = sort(cc_combined(:,2), 'ascend');
        cc_combined = cc_combined(sort_idx, :);

        n_clean_rows = size(cc,1);
        cc_spikes = [];
        if isfield(timesData, 'spikes') && size(timesData.spikes,1) == n_clean_rows
            cc_spikes = double(timesData.spikes);
        else
            cc_spikes = nan(n_clean_rows, size(fullData.spikes,2));
            valid = cc_to_clean > 0;
            cc_spikes(valid, :) = clean_spikes(cc_to_clean(valid), :);
            if any(~valid)
                cc_spikes = [];
            end
        end

        outData = strip_row_aligned_fields(timesData, n_clean_rows, {'cluster_class', 'spikes'});
        outData.cluster_class = cc_combined;
        if ~isempty(cc_spikes)
            combined_spikes = [cc_spikes; dirty_spikes];
            combined_spikes = combined_spikes(sort_idx, :);
            outData.spikes = combined_spikes;
        else
            if isfield(outData, 'spikes')
                outData = rmfield(outData, 'spikes');
            end
        end

        led_assign_info = struct();
        led_assign_info.created_at = datestr(now, 'yyyy-mm-dd HH:MM:SS');
        led_assign_info.distance_threshold_sd = opts.distance_threshold_sd;
        led_assign_info.assign_unlabeled_clean = opts.assign_unlabeled_clean;
        led_assign_info.clean_distance_threshold_sd = opts.clean_distance_threshold_sd;
        led_assign_info.template_mode = 'wave_clus_force_membership';
        led_assign_info.force_template_type = force_used.template_type;
        led_assign_info.force_feature = force_used.force_feature;
        led_assign_info.force_template_sdnum = force_used.template_sdnum;
        led_assign_info.force_template_k = force_used.template_k;
        led_assign_info.force_template_k_min = force_used.template_k_min;
        led_assign_info.time_transform_mode = mode_label;
        led_assign_info.cc_match_count = n_match;
        led_assign_info.cc_total = size(cc,1);
        led_assign_info.n_clean = numel(clean_index);
        led_assign_info.n_clean_unassigned_before = n_clean_unassigned_before;
        led_assign_info.n_clean_reassigned = sum(clean_reassigned_labels > 0);
        led_assign_info.n_cc_rows_reassigned = n_cc_clean_reassigned;
        led_assign_info.n_clean_unassigned_after = n_clean_unassigned_after;
        led_assign_info.clean_reassigned_best_distance = clean_reassigned_best_dist(:);
        led_assign_info.n_dirty = numel(dirty_index);
        led_assign_info.n_dirty_assigned = sum(dirty_labels > 0);
        led_assign_info.n_dirty_unassigned = sum(dirty_labels == 0);
        led_assign_info.template_stats = template_stats;
        led_assign_info.dirty_best_distance = dirty_best_dist(:);
        outData.led_assignment = led_assign_info;

        combinedFile = strrep(timesFile, '.mat', '_with_led_assignments.mat');
        save(combinedFile, '-struct', 'outData', '-v7.3');

        if opts.overwrite_main_times
            cleanBackup = strrep(timesFile, '.mat', '_clean_only.mat');
            if ~exist(cleanBackup, 'file')
                copyfile(timesFile, cleanBackup);
            end
            save(timesFile, '-struct', 'outData', '-v7.3');
            fprintf('\n[assign] Channel %d\n', ch);
            fprintf('  Overwrote main times file: %s\n', timesFile);
            fprintf('  Clean-only backup: %s\n', cleanBackup);
        else
            fprintf('\n[assign] Channel %d\n', ch);
            fprintf('  Wrote combined file: %s\n', combinedFile);
        end

        reportFile = strrep(timesFile, '.mat', '_dirty_assignment.mat');
        assignment_report = struct(); %#ok<NASGU>
        assignment_report.channel = ch;
        assignment_report.created_at = datestr(now, 'yyyy-mm-dd HH:MM:SS');
        assignment_report.times_file_source = timesFile;
        assignment_report.times_file_combined = combinedFile;
        assignment_report.dirty_mask = dirty_mask;
        assignment_report.dirty_index = dirty_index(:);
        assignment_report.dirty_labels = dirty_labels(:);
        assignment_report.dirty_best_distance = dirty_best_dist(:);
        assignment_report.distance_threshold_sd = opts.distance_threshold_sd;
        assignment_report.assign_unlabeled_clean = opts.assign_unlabeled_clean;
        assignment_report.clean_distance_threshold_sd = opts.clean_distance_threshold_sd;
        assignment_report.clean_labels_initial = clean_labels(:);
        assignment_report.clean_labels_after_reassign = clean_labels_for_templates(:);
        assignment_report.clean_reassigned_labels = clean_reassigned_labels(:);
        assignment_report.clean_reassigned_best_distance = clean_reassigned_best_dist(:);
        assignment_report.template_stats = template_stats;
        assignment_report.force_used = force_used;
        save(reportFile, 'assignment_report');

        if opts.save_assignment_png
            pngFile = strrep(timesFile, '.mat', '_assignment_summary.png');
            try
                save_assignment_summary_png(ch, pngFile, cc_original(:,1), ...
                    cc_clean_reassigned(:,1), dirty_best_dist, dirty_labels, ...
                    opts.distance_threshold_sd, n_clean_unassigned_before, ...
                    n_clean_unassigned_after, sum(clean_reassigned_labels > 0));
                fprintf('  Assignment PNG: %s\n', pngFile);
            catch MEpng
                warning('Channel %d: failed to save assignment PNG (%s).', ch, MEpng.message);
            end
        end

        fprintf('  Clean cluster0: before=%d | reassigned=%d | after=%d\n', ...
            n_clean_unassigned_before, sum(clean_reassigned_labels > 0), n_clean_unassigned_after);
        fprintf('  Dirty spikes: %d | assigned: %d | unassigned: %d\n', ...
            numel(dirty_index), sum(dirty_labels > 0), sum(dirty_labels == 0));
        fprintf('  Assignment report: %s\n', reportFile);
    end
end

function [master, masterPath, masterName] = load_master(masterFile)
    if ~exist(masterFile, 'file')
        error('Master file not found: %s', masterFile);
    end
    info = whos('-file', masterFile);
    names = {info.name};

    required_vars = {'Eventstime', 'EventTag', 'nChannels', 'FsPlexon'};
    optional_vars = {'firstADsamples', 'firstADsample'};
    load_vars = [required_vars, optional_vars(ismember(optional_vars, names))];

    master = load(masterFile, load_vars{:});
    if ~isfield(master, 'Eventstime') || ~isfield(master, 'EventTag')
        error('Master file missing Eventstime/EventTag: %s', masterFile);
    end
    if ~isfield(master, 'nChannels') || isempty(master.nChannels)
        error('Master file missing nChannels: %s', masterFile);
    end

    [masterPath, masterName, ~] = fileparts(masterFile);
end

function ensure_sr_in_mat_file(matFile, sr, debug_mode)
    if nargin < 3
        debug_mode = false;
    end
    if isempty(sr) || ~isfinite(sr) || sr <= 0 || ~exist(matFile, 'file')
        return;
    end

    info = whos('-file', matFile);
    names = {info.name};
    if ismember('sr', names)
        return;
    end

    tmp = struct('sr', double(sr));
    save(matFile, '-struct', 'tmp', '-append');
    if debug_mode
        fprintf('  [debug] Added sr=%g to %s\n', sr, matFile);
    end
end

function channelDir = resolve_channel_dir(masterPath, masterName)
    c1 = fullfile(masterPath, [masterName '_channels']);
    c2 = fullfile(masterPath, masterName);
    if exist(c1, 'dir')
        channelDir = c1;
    elseif exist(c2, 'dir')
        channelDir = c2;
    else
        error('No channel directory found next to master file.');
    end
end

function chFile = resolve_channel_file(channelDir, ch)
    candidates = {
        fullfile(channelDir, sprintf('channel_%d_filtered_CAR.mat', ch))
        fullfile(channelDir, sprintf('channel_%d_filtered_CARfiltered.mat', ch))
        fullfile(channelDir, sprintf('channel_%d.mat', ch))
    };
    chFile = pick_first_existing(candidates);
end

function [index_raw, index_field] = get_index_field(S, context_file)
    if isfield(S, 'index')
        index_raw = double(S.index(:));
        index_field = 'index';
        return;
    end
    if isfield(S, 'indexes')
        index_raw = double(S.indexes(:));
        index_field = 'indexes';
        return;
    end
    error('No index/indexes field found in %s', context_file);
end

function out = subset_spike_struct(in, keep_mask, n_rows)
    out = in;
    fn = fieldnames(in);
    for i = 1:numel(fn)
        f = fn{i};
        v = in.(f);

        if isnumeric(v) || islogical(v)
            if isvector(v) && numel(v) == n_rows
                if isrow(v)
                    out.(f) = v(keep_mask(:)');
                else
                    out.(f) = v(keep_mask(:));
                end
            elseif ~isvector(v) && size(v,1) == n_rows
                out.(f) = v(keep_mask(:), :);
            end
        elseif iscell(v)
            if isvector(v) && numel(v) == n_rows
                out.(f) = v(keep_mask(:));
            elseif ~isvector(v) && size(v,1) == n_rows
                out.(f) = v(keep_mask(:), :);
            end
        end
    end
end

function out = strip_row_aligned_fields(in, n_rows, keep)
    out = in;
    fn = fieldnames(in);
    drop = {};
    for i = 1:numel(fn)
        f = fn{i};
        if any(strcmp(f, keep))
            continue;
        end
        v = in.(f);
        is_row_aligned = false;
        if isnumeric(v) || islogical(v)
            if isvector(v) && numel(v) == n_rows
                is_row_aligned = true;
            elseif ~isvector(v) && size(v,1) == n_rows
                is_row_aligned = true;
            end
        elseif iscell(v)
            if isvector(v) && numel(v) == n_rows
                is_row_aligned = true;
            elseif ~isvector(v) && size(v,1) == n_rows
                is_row_aligned = true;
            end
        end
        if is_row_aligned
            drop{end+1} = f; %#ok<AGROW>
        end
    end
    if ~isempty(drop)
        out = rmfield(out, drop);
    end
end

function spkFile = spike_file_from_channel_file(chFile)
    [d, n, ~] = fileparts(chFile);
    spkFile = fullfile(d, [n '_spikes.mat']);
end

function spkFile = resolve_spike_file_after_get_spikes(chFile)
    expected = spike_file_from_channel_file(chFile);
    if exist(expected, 'file')
        spkFile = expected;
        return;
    end

    [d, n, ~] = fileparts(chFile);
    patterns = {
        [n '_spikes.mat']
        [n '*_spikes.mat']
        [n '*spikes*.mat']
    };

    cand = {};
    for i = 1:numel(patterns)
        a = dir(fullfile(d, patterns{i}));
        b = dir(fullfile(pwd, patterns{i}));
        cand = [cand; to_paths(a, d); to_paths(b, pwd)]; %#ok<AGROW>
    end
    cand = unique(cand, 'stable');

    if isempty(cand)
        list_d = dir(fullfile(d, '*.mat'));
        preview = strjoin({list_d.name}, ', ');
        if numel(preview) > 300
            preview = [preview(1:300) ' ...'];
        end
        error(['Expected spike file not found after Get_spikes.\n' ...
            'Expected: %s\nDirectory scanned: %s\nFound MAT files: %s'], ...
            expected, d, preview);
    end

    % Choose the most recently modified candidate.
    chosen = cand{1};
    chosen_time = get_mtime(chosen);
    for i = 2:numel(cand)
        t = get_mtime(cand{i});
        if t > chosen_time
            chosen = cand{i};
            chosen_time = t;
        end
    end

    if ~strcmpi(chosen, expected)
        try
            copyfile(chosen, expected);
            fprintf('  Found non-canonical spikes file, copied to expected name:\n');
            fprintf('    source: %s\n', chosen);
            fprintf('    target: %s\n', expected);
            spkFile = expected;
            return;
        catch
            warning('Using non-canonical spike file path: %s', chosen);
            spkFile = chosen;
            return;
        end
    end

    spkFile = chosen;
end

function paths = to_paths(dir_struct, root_dir)
    if isempty(dir_struct)
        paths = {};
        return;
    end
    paths = cell(numel(dir_struct), 1);
    for i = 1:numel(dir_struct)
        paths{i} = fullfile(root_dir, dir_struct(i).name);
    end
end

function t = get_mtime(f)
    d = dir(f);
    if isempty(d)
        t = -inf;
    else
        t = d(1).datenum;
    end
end

function tFile = times_file_from_spike_file(spkFile)
    [d, n, ~] = fileparts(spkFile);
    if endsWith(n, '_spikes')
        base = n(1:end-7);
    else
        base = n;
    end
    tFile = fullfile(d, ['times_' base '.mat']);
end

function tFile = resolve_times_file_after_do_clustering(spkFile)
    expected = times_file_from_spike_file(spkFile);
    if exist(expected, 'file')
        tFile = expected;
        return;
    end

    [d, n, ~] = fileparts(spkFile);
    base_no_spikes = regexprep(n, '_spikes.*$', '');
    ch_token = regexp(base_no_spikes, 'channel_\d+', 'match', 'once');
    if isempty(ch_token)
        ch_token = base_no_spikes;
    end

    patterns = {
        ['times_' n '.mat']
        ['times_' n '*.mat']
        ['times_' base_no_spikes '.mat']
        ['times_' base_no_spikes '*.mat']
        ['times_*' ch_token '*.mat']
    };

    cand = {};
    for i = 1:numel(patterns)
        a = dir(fullfile(d, patterns{i}));
        b = dir(fullfile(pwd, patterns{i}));
        cand = [cand; to_paths(a, d); to_paths(b, pwd)]; %#ok<AGROW>
    end
    cand = unique(cand, 'stable');

    if isempty(cand)
        tFile = '';
        return;
    end

    chosen = cand{1};
    chosen_time = get_mtime(chosen);
    for i = 2:numel(cand)
        t = get_mtime(cand{i});
        if t > chosen_time
            chosen = cand{i};
            chosen_time = t;
        end
    end
    tFile = chosen;
end

function run_get_spikes(chFile, par, debug_mode)
    if nargin < 3
        debug_mode = false;
    end

    [chDir, chName, chExt] = fileparts(chFile);
    if isempty(chDir), chDir = pwd; end
    chBase = [chName chExt];
    oldDir = pwd;
    c = onCleanup(@() cd(oldDir)); %#ok<NASGU>
    cd(chDir);

    if isempty(par)
        if debug_mode
            fprintf('  [debug] Get_spikes call (cwd=%s): Get_spikes({%s})\n', chDir, chBase);
        end
        Get_spikes({chBase});
    else
        par_runtime = par;
        if isfield(par_runtime, 'cont_segment')
            par_runtime.cont_segment = false;
        end
        if isfield(par_runtime, 'segments_length') && (~isfinite(par_runtime.segments_length) || par_runtime.segments_length > 1e6)
            par_runtime.segments_length = 5;
        end

        try
            if debug_mode
                fprintf('  [debug] Get_spikes call (cwd=%s): Get_spikes({%s}, ''par'', par_runtime)\n', chDir, chBase);
            end
            Get_spikes({chBase}, 'par', par_runtime);
        catch ME1
            if debug_mode
                warning(['Get_spikes first attempt failed (%s). ' ...
                    'Retrying with single-file call style.'], ME1.message);
            end
            try
                if debug_mode
                    fprintf('  [debug] Get_spikes retry: Get_spikes(%s, ''par'', par_runtime)\n', chBase);
                end
                Get_spikes(chBase, 'par', par_runtime);
            catch ME2
                if contains(ME2.message, 'raw data is already fully loaded', 'IgnoreCase', true)
                    if debug_mode
                        warning(['Get_spikes failed due to segmented-read mode on fully loaded data. ' ...
                            'Retrying once with finite segment length for compatibility.']);
                    end
                    par_retry = par_runtime;
                    if isfield(par_retry, 'segments_length') && (~isfinite(par_retry.segments_length) || par_retry.segments_length > 1e6)
                        par_retry.segments_length = 5;
                    end
                    if debug_mode
                        fprintf('  [debug] Get_spikes retry2: segments_length=%g\n', par_retry.segments_length);
                        fprintf('  [debug] Get_spikes retry2 call: Get_spikes(%s, ''par'', par_retry)\n', chBase);
                    end
                    Get_spikes(chBase, 'par', par_retry);
                else
                    rethrow(ME2);
                end
            end
        end
    end
end

function print_spike_file_listing(dirPath, phase_label)
    if ~exist(dirPath, 'dir')
        fprintf('  [debug] (%s) dir missing: %s\n', phase_label, dirPath);
        return;
    end

    L = dir(fullfile(dirPath, '*spikes*.mat'));
    if isempty(L)
        fprintf('  [debug] (%s) no *spikes*.mat in %s\n', phase_label, dirPath);
        return;
    end

    fprintf('  [debug] (%s) spikes files in %s:\n', phase_label, dirPath);
    [~, order] = sort([L.datenum], 'descend');
    L = L(order);
    nShow = min(numel(L), 8);
    for i = 1:nShow
        fprintf('    %s | %s\n', datestr(L(i).datenum, 'yyyy-mm-dd HH:MM:SS'), L(i).name);
    end
end

function run_do_clustering(spkFile, par, debug_mode)
    if nargin < 3
        debug_mode = false;
    end

    [spkDir, spkName, spkExt] = fileparts(spkFile);
    if isempty(spkDir), spkDir = pwd; end
    spkBase = [spkName spkExt];
    oldDir = pwd;
    c = onCleanup(@() cd(oldDir)); %#ok<NASGU>
    cd(spkDir);

    if debug_mode
        fprintf('  [debug] Do_clustering call (cwd=%s): %s\n', spkDir, spkBase);
    end

    % Some wave_clus Windows builds emit this non-actionable warning due
    % backslash handling in sprintf format strings.
    ws = warning('query', 'MATLAB:dispatcher:invalidEscapeSequence');
    wc = onCleanup(@() warning(ws.state, 'MATLAB:dispatcher:invalidEscapeSequence')); %#ok<NASGU>
    warning('off', 'MATLAB:dispatcher:invalidEscapeSequence');

    if isempty(par)
        Do_clustering(spkBase);
    else
        Do_clustering(spkBase, 'par', par);
    end
end

function sr = get_sampling_rate(master, par)
    sr = [];
    if ~isempty(par) && isfield(par, 'sr') && isfinite(par.sr) && par.sr > 0
        sr = double(par.sr);
    elseif isfield(master, 'FsPlexon') && isfinite(master.FsPlexon) && master.FsPlexon > 0
        sr = double(master.FsPlexon);
    else
        sr = 40000;
        warning('Sampling rate unknown; defaulting to 40000 Hz.');
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

function led_windows = build_led_windows(event_s, event_tag, led_tags, led_duration_s, post_buffer_s)
    if isempty(event_s) || isempty(event_tag)
        led_windows = zeros(0, 2);
        return;
    end

    led_windows = zeros(0, 2);
    for i = 1:numel(led_tags)
        tag = led_tags(i);
        dur = led_duration_s(i);
        t = event_s(event_tag == tag);
        if isempty(t)
            continue;
        end
        w = [double(t(:)) double(t(:) + dur + post_buffer_s)];
        led_windows = [led_windows; w]; %#ok<AGROW>
    end

    led_windows = merge_windows(led_windows);
end

function merged = merge_windows(windows)
    if isempty(windows)
        merged = zeros(0, 2);
        return;
    end
    w = sortrows(double(windows), 1);
    merged = w(1, :);
    for i = 2:size(w, 1)
        if w(i,1) <= merged(end,2)
            merged(end,2) = max(merged(end,2), w(i,2));
        else
            merged(end+1, :) = w(i, :); %#ok<AGROW>
        end
    end
end

function [spike_s, mode] = infer_spike_times_seconds(index_raw, sr, first_sample_s, event_s, requested_mode)
    if nargin < 5 || isempty(requested_mode)
        requested_mode = 'auto';
    end
    requested_mode = normalize_time_reference_mode(requested_mode);

    idx = double(index_raw(:));
    cands = {
        idx / 1000, 'index_ms'
        idx / sr, 'index_samples'
        idx / 1000 + first_sample_s, 'index_ms+firstAD'
        idx / sr + first_sample_s, 'index_samples+firstAD'
    };

    if ~strcmp(requested_mode, 'auto')
        labels = cands(:,2);
        k = find(strcmp(labels, requested_mode), 1, 'first');
        if isempty(k)
            error('Unsupported time reference mode: %s', requested_mode);
        end
        spike_s = cands{k,1};
        mode = cands{k,2};
        return;
    end

    if isempty(event_s)
        spike_s = cands{1,1};
        mode = cands{1,2};
        return;
    end

    lo = min(event_s) - 2;
    hi = max(event_s) + 2;
    med_ev = median(event_s);

    scores = -inf(size(cands,1), 1);
    for k = 1:size(cands,1)
        t = cands{k,1};
        frac = mean(t >= lo & t <= hi);
        med_dist = abs(median(t) - med_ev);
        scores(k) = frac - 1e-3 * med_dist;
    end
    [~, best] = max(scores);
    spike_s = cands{best,1};
    mode = cands{best,2};
end

function mode = normalize_time_reference_mode(mode_in)
    mode = lower(strtrim(char(mode_in)));
    mode = strrep(mode, ' ', '');
    valid = {'auto', 'index_ms', 'index_samples', 'index_ms+firstad', 'index_samples+firstad'};
    if ~ismember(mode, valid)
        error(['time_reference_mode must be one of: auto, index_ms, index_samples, ' ...
            'index_ms+firstAD, index_samples+firstAD']);
    end
    if strcmp(mode, 'index_ms+firstad')
        mode = 'index_ms+firstAD';
    elseif strcmp(mode, 'index_samples+firstad')
        mode = 'index_samples+firstAD';
    end
end

function mask = in_any_window(t, windows)
    t = double(t(:));
    mask = false(size(t));
    if isempty(t) || isempty(windows)
        return;
    end
    for i = 1:size(windows,1)
        mask = mask | (t >= windows(i,1) & t <= windows(i,2));
    end
end

function [cc_to_clean, mode_id, mode_label, n_match] = map_cc_to_clean_indices( ...
    clean_index_raw, cc_time_ms, sr, first_sample_s)

    clean_index_raw = double(clean_index_raw(:));
    cc_time_ms = double(cc_time_ms(:));

    cand = {
        clean_index_raw, 'index_ms'
        clean_index_raw * (1000 / sr), 'index_samples_to_ms'
        clean_index_raw + first_sample_s * 1000, 'index_ms+firstAD_ms'
        clean_index_raw * (1000 / sr) + first_sample_s * 1000, 'index_samples_to_ms+firstAD_ms'
    };

    cc_keys = round(cc_time_ms * 1000); % micro-ms keys
    best_score = -1;
    best_loc = zeros(size(cc_keys));
    mode_id = 1;
    mode_label = cand{1,2};

    for k = 1:size(cand,1)
        clean_keys = round(cand{k,1} * 1000);
        [tf, loc] = ismember(cc_keys, clean_keys);
        score = sum(tf);
        if score > best_score
            best_score = score;
            best_loc = loc;
            mode_id = k;
            mode_label = cand{k,2};
        end
    end

    cc_to_clean = best_loc;
    n_match = sum(cc_to_clean > 0);
end

function time_ms = transform_index_to_ms(index_raw, sr, first_sample_s, mode_id)
    idx = double(index_raw(:));
    switch mode_id
        case 1
            time_ms = idx;
        case 2
            time_ms = idx * (1000 / sr);
        case 3
            time_ms = idx + first_sample_s * 1000;
        case 4
            time_ms = idx * (1000 / sr) + first_sample_s * 1000;
        otherwise
            time_ms = idx;
    end
end

function force_par = resolve_force_params(timesData, template_sdnum)
    force_par = struct();

    if isfield(timesData, 'par') && isstruct(timesData.par)
        src = timesData.par;
    else
        src = [];
        if exist('set_parameters_custom', 'file') == 2
            src = set_parameters_custom();
        elseif exist('set_parameters', 'file') == 2
            src = set_parameters();
        end
    end

    if isempty(src), src = struct(); end

    force_par.template_type = get_struct_field_or_default(src, 'template_type', 'center');
    force_par.template_k = get_struct_field_or_default(src, 'template_k', 10);
    force_par.template_k_min = get_struct_field_or_default(src, 'template_k_min', 10);
    force_par.force_feature = get_struct_field_or_default(src, 'force_feature', 'spk');
    force_par.features = get_struct_field_or_default(src, 'features', 'wav');
    force_par.scales = get_struct_field_or_default(src, 'scales', 4);
    force_par.min_inputs = get_struct_field_or_default(src, 'min_inputs', 10);
    force_par.max_inputs = get_struct_field_or_default(src, 'max_inputs', 0.75);

    force_par.template_type = lower(strtrim(char(force_par.template_type)));
    force_par.force_feature = lower(strtrim(char(force_par.force_feature)));
    if ~ismember(force_par.template_type, {'nn', 'center', 'ml', 'mahal'})
        warning('Unknown template_type "%s"; using center.', force_par.template_type);
        force_par.template_type = 'center';
    end
    if ~ismember(force_par.force_feature, {'spk', 'wav'})
        warning('Unknown force_feature "%s"; using spk.', force_par.force_feature);
        force_par.force_feature = 'spk';
    end

    force_par.template_k = max(1, round(double(force_par.template_k)));
    force_par.template_k_min = max(1, round(double(force_par.template_k_min)));
    force_par.template_sdnum = double(template_sdnum);
end

function [target_labels, target_best_dist, template_stats, force_used] = assign_spikes_waveclus_force( ...
    template_spikes, template_labels, target_spikes, force_par)

    template_spikes = double(template_spikes);
    target_spikes = double(target_spikes);
    template_labels = double(template_labels(:));

    target_labels = zeros(size(target_spikes,1), 1);
    target_best_dist = inf(size(target_spikes,1), 1);
    force_used = force_par;

    cids = unique(template_labels(template_labels > 0));
    template_stats = struct('cluster_id', {}, 'n_spikes', {});
    for i = 1:numel(cids)
        template_stats(i).cluster_id = cids(i);
        template_stats(i).n_spikes = sum(template_labels == cids(i));
    end
    if isempty(cids) || isempty(target_spikes)
        return;
    end

    in_mask = template_labels > 0;
    if ~any(in_mask)
        return;
    end

    f_in_all = template_spikes;
    f_out = target_spikes;

    if strcmp(force_par.force_feature, 'wav')
        if exist('wave_features', 'file') == 2
            try
                all_spk = [f_in_all; f_out];
                all_feat = wave_features(all_spk, force_par);
                n_in = size(f_in_all, 1);
                f_in_all = double(all_feat(1:n_in, :));
                f_out = double(all_feat((n_in+1):end, :));
            catch MEwf
                warning('wave_features failed (%s). Falling back to spike waveform forcing.', MEwf.message);
                force_used.force_feature = 'spk';
            end
        else
            warning('wave_features not found. Falling back to spike waveform forcing.');
            force_used.force_feature = 'spk';
        end
    end

    f_in = f_in_all(in_mask, :);
    class_in = template_labels(in_mask)';
    if isempty(f_in) || isempty(class_in)
        return;
    end

    class_out = force_membership_wc_local(f_in, class_in, f_out, force_used);
    target_labels = double(class_out(:));
    target_best_dist = force_best_distance_sd_units(f_in, class_in, f_out, force_used);
end

function class_out = force_membership_wc_local(f_in, class_in, f_out, par)
    nspk = size(f_out, 1);
    class_out = zeros(1, nspk);
    if isempty(f_in) || isempty(class_in) || nspk == 0
        return;
    end

    switch lower(par.template_type)
        case 'nn'
            sdnum = par.template_sdnum;
            k = par.template_k;
            k_min = par.template_k_min;
            sd = sqrt(sum(var(f_in, 1))) * ones(1, size(f_in, 1));
            for i = 1:nspk
                nn = nearest_neighbor_local(f_out(i,:), f_in, sdnum * sd, ...
                    Inf * ones(size(f_in)), Inf, k);
                if any(nn)
                    winner = mode(class_in(nn));
                    if nnz(class_in(nn) == winner) < k_min
                        class_out(i) = 0;
                    else
                        class_out(i) = winner;
                    end
                else
                    class_out(i) = 0;
                end
            end

        case 'center'
            [centers, sd] = build_templates_local(class_in, f_in);
            sdnum = par.template_sdnum;
            for i = 1:nspk
                class_out(i) = nearest_neighbor_local(f_out(i,:), centers, sdnum * sd);
            end

        case 'ml'
            [mu, inv_sigma] = fit_gaussian_local(f_in, class_in);
            for i = 1:nspk
                class_out(i) = ml_gaussian_local(f_out(i,:), mu, inv_sigma);
            end

        case 'mahal'
            classes = unique(class_in);
            mdistance = inf(numel(classes), nspk);
            maxdist = zeros(1, numel(classes));
            for ci = 1:numel(classes)
                c = classes(ci);
                fi = f_in(class_in == c, :);
                if size(fi,1) < 2 || exist('mahal', 'file') ~= 2
                    warning('mahal forcing unavailable for class %d; falling back to center forcing.', c);
                    p2 = par;
                    p2.template_type = 'center';
                    class_out = force_membership_wc_local(f_in, class_in, f_out, p2);
                    return;
                end
                mdistance(ci,:) = mahal(f_out, fi);
                maxdist(ci) = sqrt(mean(mahal(fi, fi)));
            end
            sdnum = par.template_sdnum;
            for i = 1:nspk
                [d, winner] = min(mdistance(:, i));
                if isfinite(d) && sqrt(d) < sdnum * maxdist(winner)
                    class_out(i) = classes(winner);
                end
            end

        otherwise
            warning('Unknown template_type "%s"; falling back to center forcing.', par.template_type);
            p2 = par;
            p2.template_type = 'center';
            class_out = force_membership_wc_local(f_in, class_in, f_out, p2);
    end
end

function best_sd = force_best_distance_sd_units(f_in, class_in, f_out, par)
    nspk = size(f_out, 1);
    best_sd = inf(nspk, 1);
    if isempty(f_in) || isempty(class_in) || nspk == 0
        return;
    end

    switch lower(par.template_type)
        case 'center'
            [centers, sd] = build_templates_local(class_in, f_in);
            sd = double(sd(:));
            for i = 1:nspk
                d = sqrt(sum((centers - f_out(i,:)).^2, 2));
                z = inf(size(d));
                ok = sd > 0 & isfinite(sd);
                z(ok) = d(ok) ./ sd(ok);
                if any(isfinite(z))
                    best_sd(i) = min(z);
                end
            end

        case 'nn'
            sd = sqrt(sum(var(f_in, 1)));
            if ~isfinite(sd) || sd <= 0
                return;
            end
            for i = 1:nspk
                d = sqrt(sum((f_in - f_out(i,:)).^2, 2));
                best_sd(i) = min(d) / sd;
            end

        case 'mahal'
            classes = unique(class_in);
            mdistance = inf(numel(classes), nspk);
            maxdist = zeros(1, numel(classes));
            for ci = 1:numel(classes)
                c = classes(ci);
                fi = f_in(class_in == c, :);
                if size(fi,1) < 2 || exist('mahal', 'file') ~= 2
                    return;
                end
                mdistance(ci,:) = mahal(f_out, fi);
                maxdist(ci) = sqrt(mean(mahal(fi, fi)));
            end
            for i = 1:nspk
                [d, winner] = min(mdistance(:, i));
                if isfinite(d) && maxdist(winner) > 0
                    best_sd(i) = sqrt(d) / maxdist(winner);
                end
            end

        case 'ml'
            best_sd = nan(nspk, 1);
    end
end

function index = nearest_neighbor_local(x, vectors, maxdist, varargin)
    if isempty(vectors)
        index = 0;
        return;
    end
    distances = sqrt(sum((ones(size(vectors,1),1) * x - vectors).^2, 2)');
    conforming = find(distances < maxdist);
    if ~isempty(varargin)
        pointdist = varargin{1};
        if numel(varargin) > 1
            pointlimit = varargin{2};
        else
            pointlimit = Inf;
        end
        pointwise_conforming = [];
        for i = 1:size(vectors,1)
            if sum(abs(x - vectors(i,:)) > pointdist(i,:)) < pointlimit
                pointwise_conforming = [pointwise_conforming i]; %#ok<AGROW>
            end
        end
        conforming = intersect(conforming, pointwise_conforming);
    end
    if isempty(conforming)
        index = 0;
    else
        if numel(varargin) > 2
            k = varargin{3};
            [~, ii] = sort(distances(conforming));
            ii = ii(1:min(numel(ii), k));
        else
            [~, ii] = min(distances(conforming));
        end
        index = conforming(ii);
    end
end

function [templates, maxdist, pointdist] = build_templates_local(classes, features)
    max_class = max(classes);
    feature_dim = size(features, 2);
    templates = zeros(max_class, feature_dim);
    maxdist = zeros(1, max_class);
    pointdist = zeros(max_class, feature_dim);
    for i = 1:max_class
        fi = features(classes == i, :);
        if isempty(fi), continue; end
        templates(i,:) = mean(fi, 1);
        maxdist(i) = sqrt(sum(var(fi, 1)));
        pointdist(i,:) = sqrt(var(fi, 1));
    end
end

function [mu, inv_sigma] = fit_gaussian_local(x, class)
    N = max(class);
    mu = zeros(N, size(x,2));
    inv_sigma = zeros(size(x,2), size(x,2), N);
    for i = 1:N
        xi = x(class == i, :);
        if size(xi,1) < 2
            inv_sigma(:,:,i) = eye(size(x,2));
            if ~isempty(xi)
                mu(i,:) = mean(xi, 1);
            end
        else
            mu(i,:) = mean(xi, 1);
            inv_sigma(:,:,i) = pinv(cov(xi));
        end
    end
end

function index = ml_gaussian_local(x, mu, inv_sigma)
    N = size(mu,1);
    if N == 0
        index = 0;
        return;
    end
    p = zeros(1, N);
    for i = 1:N
        q = (x - mu(i,:)) * inv_sigma(:,:,i) * (x - mu(i,:))';
        p(i) = sqrt(max(det(inv_sigma(:,:,i)), 0)) * exp(-0.5 * q);
    end
    [~, index] = max(p);
end

function v = get_struct_field_or_default(S, f, default_v)
    if isstruct(S) && isfield(S, f) && ~isempty(S.(f))
        v = S.(f);
    else
        v = default_v;
    end
end

function save_assignment_summary_png(ch, png_file, cc_before_labels, cc_after_clean_labels, ...
    dirty_best_dist, dirty_labels, dirty_threshold_sd, ...
    n_clean_unassigned_before, n_clean_unassigned_after, n_clean_reassigned)

    fig = figure('Color', 'w', 'Visible', 'off', 'Position', [80 80 1050 760]);
    tl = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, sprintf('Channel %d Step2 Assignment Summary', ch), 'FontWeight', 'bold');

    ax1 = nexttile(tl, 1);
    hold(ax1, 'on');
    ids = unique([cc_before_labels(:); cc_after_clean_labels(:)]);
    ids = ids(:)';
    if isempty(ids), ids = 0; end
    ids = sort(ids);
    x = 1:numel(ids);
    c_before = arrayfun(@(k) sum(cc_before_labels == k), ids);
    c_after = arrayfun(@(k) sum(cc_after_clean_labels == k), ids);
    b = bar(ax1, x, [c_before(:) c_after(:)], 'grouped'); %#ok<NASGU>
    set(ax1, 'XTick', x, 'XTickLabel', arrayfun(@num2str, ids, 'UniformOutput', false));
    xlabel(ax1, 'Cluster ID');
    ylabel(ax1, 'Spike count');
    title(ax1, 'Clean labels before vs after clean reassignment');
    legend(ax1, {'Before', 'After'}, 'Location', 'best');
    grid(ax1, 'on');
    box(ax1, 'on');

    ax2 = nexttile(tl, 2);
    hold(ax2, 'on');
    if isempty(dirty_best_dist)
        text(ax2, 0.5, 0.5, 'No dirty spikes', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center');
        axis(ax2, 'off');
    else
        histogram(ax2, dirty_best_dist, 50, 'FaceColor', [0.2 0.5 0.9], 'EdgeColor', 'none');
        xline(ax2, dirty_threshold_sd, '--r', sprintf('thr=%.2f', dirty_threshold_sd), 'LineWidth', 1.5);
        xlabel(ax2, 'Best normalized distance');
        ylabel(ax2, 'Count');
        title(ax2, 'Dirty spike template distances');
        grid(ax2, 'on');
        box(ax2, 'on');
    end

    ax3 = nexttile(tl, 3);
    hold(ax3, 'on');
    n_dirty = numel(dirty_labels);
    n_dirty_assigned = sum(dirty_labels > 0);
    n_dirty_unassigned = sum(dirty_labels == 0);
    bar(ax3, [1 2], [n_dirty_assigned n_dirty_unassigned], 0.6, ...
        'FaceColor', [0.3 0.7 0.3], 'EdgeColor', 'none');
    set(ax3, 'XTick', [1 2], 'XTickLabel', {'Dirty assigned', 'Dirty unassigned'});
    ylabel(ax3, 'Count');
    title(ax3, sprintf('Dirty assignment (n=%d)', n_dirty));
    grid(ax3, 'on');
    box(ax3, 'on');

    ax4 = nexttile(tl, 4);
    axis(ax4, 'off');
    txt = sprintf(['Clean cluster0 before: %d\n' ...
        'Clean reassigned: %d\n' ...
        'Clean cluster0 after: %d\n' ...
        'Dirty assigned: %d\n' ...
        'Dirty unassigned: %d\n' ...
        'Dirty threshold: %.2f SD'], ...
        n_clean_unassigned_before, n_clean_reassigned, n_clean_unassigned_after, ...
        n_dirty_assigned, n_dirty_unassigned, dirty_threshold_sd);
    text(ax4, 0.02, 0.98, txt, 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'FontName', 'Consolas', 'FontSize', 10);

    saveas(fig, png_file);
    close(fig);
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
