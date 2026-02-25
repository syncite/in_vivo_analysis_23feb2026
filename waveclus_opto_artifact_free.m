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
%      - Save explicit *_spikes_clean.mat and *_spikes_dirty.mat files
%      - Run Do_clustering on clean spikes
%
%   2) 'assign'
%      - Load manually curated clean-only times_*.mat
%      - Assign dirty spikes by normalized template distance
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
%   'distance_threshold_sd' : dirty assignment threshold        3
%   'overwrite_main_times'  : in assign mode, replace times_*.mat [true]
%   'run_clustering'        : in prepare mode, run Do_clustering [true]
%
% Files created per channel:
%   channel_<N>_..._spikes.mat       - full spikes (from Get_spikes)
%   channel_<N>_..._spikes_clean.mat - clean-only spikes
%   channel_<N>_..._spikes_dirty.mat - dirty-only spikes
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
    addParameter(p, 'distance_threshold_sd', 3, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'overwrite_main_times', true, @islogical);
    addParameter(p, 'run_clustering', true, @islogical);
    parse(p, master_mat_file, varargin{:});
    opts = p.Results;

    opts.mode = lower(strtrim(opts.mode));
    if ~ismember(opts.mode, {'prepare', 'assign'})
        error('mode must be ''prepare'' or ''assign''.');
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

        run_get_spikes(chFile, opts.par);
        spkFileFull = spike_file_from_channel_file(chFile);
        if ~exist(spkFileFull, 'file')
            error('Expected spike file not found: %s', spkFileFull);
        end

        fullData = load(spkFileFull);
        if ~isfield(fullData, 'spikes')
            error('Spike file missing ''spikes'': %s', spkFileFull);
        end
        [index_raw, index_field] = get_index_field(fullData, spkFileFull);

        first_sample_s = get_first_sample_seconds(master, ch);
        [spike_s, ref_mode] = infer_spike_times_seconds(index_raw, sr, first_sample_s, master.Eventstime(:));

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

        spkFileClean = strrep(spkFileFull, '_spikes.mat', '_spikes_clean.mat');
        spkFileDirty = strrep(spkFileFull, '_spikes.mat', '_spikes_dirty.mat');
        splitFile = strrep(spkFileFull, '_spikes.mat', '_led_split.mat');

        cleanData = subset_spike_struct(fullData, clean_mask, n_total);
        dirtyData = subset_spike_struct(fullData, dirty_mask, n_total);

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
            run_do_clustering(spkFileClean, opts.par);
            tFileClean = times_file_from_spike_file(spkFileClean);
            tFileMain = times_file_from_spike_file(spkFileFull);
            if exist(tFileClean, 'file')
                copyfile(tFileClean, tFileMain);
                fprintf('  Clean clustering file: %s\n', tFileClean);
                fprintf('  Main times file for GUI: %s\n', tFileMain);
            else
                warning('Do_clustering finished but expected times file not found: %s', tFileClean);
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
            spkFileClean = strrep(spkFile, '_spikes.mat', '_spikes_clean.mat');
        end

        if isfield(split, 'spike_file_dirty')
            spkFileDirty = split.spike_file_dirty;
        else
            spkFileDirty = strrep(spkFile, '_spikes.mat', '_spikes_dirty.mat');
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

        clean_labels = zeros(numel(clean_index), 1);
        for i = 1:size(cc,1)
            idx = cc_to_clean(i);
            if idx > 0 && idx <= numel(clean_labels)
                clean_labels(idx) = cc(i,1);
            end
        end

        [dirty_labels, dirty_best_dist, template_stats] = assign_dirty_spikes( ...
            clean_spikes, clean_labels, dirty_spikes, opts.distance_threshold_sd);

        dirty_time_ms = transform_index_to_ms(double(dirty_index), sr, first_sample_s, mode_id);
        dirty_cc = [double(dirty_labels(:)) double(dirty_time_ms(:))];

        cc_combined = [cc; dirty_cc];
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
        led_assign_info.template_mode = 'normalized_waveform_distance';
        led_assign_info.time_transform_mode = mode_label;
        led_assign_info.cc_match_count = n_match;
        led_assign_info.cc_total = size(cc,1);
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
        assignment_report.template_stats = template_stats;
        save(reportFile, 'assignment_report');

        fprintf('  Dirty spikes: %d | assigned: %d | unassigned: %d\n', ...
            numel(dirty_index), sum(dirty_labels > 0), sum(dirty_labels == 0));
        fprintf('  Assignment report: %s\n', reportFile);
    end
end

function [master, masterPath, masterName] = load_master(masterFile)
    if ~exist(masterFile, 'file')
        error('Master file not found: %s', masterFile);
    end
    master = load(masterFile, 'Eventstime', 'EventTag', 'nChannels', ...
        'FsPlexon', 'firstADsamples', 'firstADsample');
    if ~isfield(master, 'Eventstime') || ~isfield(master, 'EventTag')
        error('Master file missing Eventstime/EventTag: %s', masterFile);
    end
    if ~isfield(master, 'nChannels') || isempty(master.nChannels)
        error('Master file missing nChannels: %s', masterFile);
    end

    [masterPath, masterName, ~] = fileparts(masterFile);
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

function tFile = times_file_from_spike_file(spkFile)
    [d, n, ~] = fileparts(spkFile);
    if endsWith(n, '_spikes')
        base = n(1:end-7);
    else
        base = n;
    end
    tFile = fullfile(d, ['times_' base '.mat']);
end

function run_get_spikes(chFile, par)
    if isempty(par)
        Get_spikes({chFile});
    else
        Get_spikes({chFile}, 'par', par);
    end
end

function run_do_clustering(spkFile, par)
    if isempty(par)
        Do_clustering(spkFile);
    else
        Do_clustering(spkFile, 'par', par);
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

function [spike_s, mode] = infer_spike_times_seconds(index_raw, sr, first_sample_s, event_s)
    idx = double(index_raw(:));
    cands = {
        idx / 1000, 'index_ms'
        idx / sr, 'index_samples'
        idx / 1000 + first_sample_s, 'index_ms+firstAD'
        idx / sr + first_sample_s, 'index_samples+firstAD'
    };

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

function [dirty_labels, dirty_best_dist, template_stats] = assign_dirty_spikes( ...
    clean_spikes, clean_labels, dirty_spikes, threshold_sd)

    clean_labels = double(clean_labels(:));
    dirty_labels = zeros(size(dirty_spikes,1), 1);
    dirty_best_dist = inf(size(dirty_spikes,1), 1);

    cids = unique(clean_labels(clean_labels > 0));
    template_stats = struct('cluster_id', {}, 'n_spikes', {});
    if isempty(cids) || isempty(dirty_spikes)
        return;
    end

    K = numel(cids);
    D = inf(size(dirty_spikes,1), K);

    for ki = 1:K
        cid = cids(ki);
        mask = clean_labels == cid;
        w = clean_spikes(mask, :);
        mu = mean(w, 1);
        sd = std(w, 0, 1);
        sd(sd < 1e-6) = 1e-6;

        z = (dirty_spikes - mu) ./ sd;
        D(:, ki) = sqrt(mean(z.^2, 2));

        template_stats(ki).cluster_id = cid;
        template_stats(ki).n_spikes = sum(mask);
    end

    [dirty_best_dist, best_idx] = min(D, [], 2);
    pass = dirty_best_dist <= threshold_sd;
    dirty_labels(pass) = cids(best_idx(pass));
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
