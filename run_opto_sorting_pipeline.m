function run_opto_sorting_pipeline(varargin)
% RUN_OPTO_SORTING_PIPELINE  One-command interactive opto spike-sorting pipeline.
%
% MATLAB 2021a-compatible interactive workflow:
%   1) Pick .plx (raw) or master .mat (already filtered)
%   2) If .plx, run extract_and_filter
%   3) Prompt strictness (template generation + assignment)
%   4) Prompt channels (default all; supports 1:4 syntax)
%   5) Generate clean templates (prepare)
%   6) Display generated PNGs
%   7) Optional wave_clus GUI refinement
%   8) Assign remaining clean + dirty spikes (assign)
%   9) Display resulting PNGs
%
% Usage:
%   run_opto_sorting_pipeline
%   run_opto_sorting_pipeline('debug', true)
%   run_opto_sorting_pipeline('sdk_path', 'D:\path\to\Plexon SDK')

    p = inputParser;
    addParameter(p, 'sdk_path', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'debug', false, @islogical);
    parse(p, varargin{:});
    opts = p.Results;
    opts.sdk_path = char(opts.sdk_path);

    fprintf('\n=== Interactive Artifact-Free Spike Sorting Pipeline ===\n');
    fprintf('Select a .plx (raw) or master .mat (already filtered).\n\n');

    [fname, fpath] = uigetfile( ...
        {'*.plx;*.mat', 'Plexon or master MAT (*.plx, *.mat)'; ...
         '*.plx', 'Plexon raw (*.plx)'; ...
         '*.mat', 'Master MAT (*.mat)'}, ...
        'Select .plx (raw) or master .mat (already filtered)');
    if isequal(fname, 0)
        fprintf('Selection canceled.\n');
        return;
    end

    selected_file = fullfile(fpath, fname);
    [sel_dir, sel_base, sel_ext] = fileparts(selected_file);
    sel_ext = lower(sel_ext);

    switch sel_ext
        case '.plx'
            fprintf('Selected raw PLX: %s\n', selected_file);
            master_mat = run_extract_and_get_master(selected_file, opts);
        case '.mat'
            fprintf('Selected MAT: %s\n', selected_file);
            master_mat = selected_file;
        otherwise
            error('Unsupported file type: %s', sel_ext);
    end

    [is_master, missing_fields] = is_master_file(master_mat);
    if ~is_master
        error(['Selected MAT is not a master file from extract_and_filter.\n' ...
            'Missing fields: %s\nSelected: %s'], strjoin(missing_fields, ', '), master_mat);
    end

    [channel_dir, all_channels] = discover_filtered_channels(master_mat);
    if isempty(all_channels)
        error('No filtered channel files found in %s', channel_dir);
    end

    [~, master_name, ~] = fileparts(master_mat);
    ch_text = format_channel_compact(all_channels);
    proceed = ask_yes_no(sprintf([ ...
        '\n%s.mat file found along with filtered .mat files for channels %s.\n' ...
        'Do you wish to proceed with generating ''clean'' template? [y/n]: '], ...
        master_name, ch_text), true);
    if ~proceed
        fprintf('Stopped by user.\n');
        return;
    end

    template_choice = prompt_template_strictness();
    assign_choice = prompt_assignment_strictness();
    channels = prompt_channels(all_channels);
    run_plan = build_run_plan(template_choice, assign_choice);

    fprintf('\nSelected channels: %s\n', mat2str(channels));
    if numel(run_plan) == 1
        [~, template_label] = build_template_params(run_plan(1).template_mode);
        [clean_thr, dirty_thr, assign_label] = map_assignment_thresholds(run_plan(1).assign_mode);
        fprintf('Template strictness: %s\n', template_label);
        fprintf('Assignment strictness: %s (clean thr=%.2f SD, dirty thr=%.2f SD)\n', ...
            assign_label, clean_thr, dirty_thr);
        compare_root = '';
    else
        compare_root = fullfile(channel_dir, ['strictness_compare_' datestr(now, 'yyyymmdd_HHMMSS')]);
        if ~exist(compare_root, 'dir')
            mkdir(compare_root);
        end
        fprintf('Compare mode: running %d full pipeline passes.\n', numel(run_plan));
        for r = 1:numel(run_plan)
            [~, template_label] = build_template_params(run_plan(r).template_mode);
            [clean_thr, dirty_thr, assign_label] = map_assignment_thresholds(run_plan(r).assign_mode);
            fprintf('  Run %d -> template=%s | assign=%s (clean %.2f SD, dirty %.2f SD)\n', ...
                r, template_label, assign_label, clean_thr, dirty_thr);
        end
        fprintf('Compare outputs folder: %s\n', compare_root);
        fprintf('Compare mode runs unattended (no GUI pause between runs).\n');
    end

    for r = 1:numel(run_plan)
        run_start = now;
        [par, template_label] = build_template_params(run_plan(r).template_mode);
        [clean_thr, dirty_thr, assign_label] = map_assignment_thresholds(run_plan(r).assign_mode);
        run_tag = sprintf('run_%02d_template_%s_assign_%s', r, lower(template_label), lower(assign_label));

        fprintf('\n============================================================\n');
        fprintf('Pipeline run %d/%d\n', r, numel(run_plan));
        fprintf('Template strictness: %s | Assignment strictness: %s\n', template_label, assign_label);
        fprintf('============================================================\n');

        fprintf('\n=== Step 1: Generate Clean Templates ===\n');
        t_prepare_start = now;
        waveclus_opto_artifact_free(master_mat, ...
            'mode', 'prepare', ...
            'channels', channels, ...
            'par', par, ...
            'time_reference_mode', 'index_ms+firstAD', ...
            'assign_unlabeled_clean', true, ...
            'clean_distance_threshold_sd', clean_thr, ...
            'distance_threshold_sd', dirty_thr, ...
            'save_assignment_png', true, ...
            'debug', opts.debug);

        fprintf('\nDisplaying PNGs generated in step 1...\n');
        show_recent_pngs(channel_dir, t_prepare_start, sprintf('Step 1 PNGs (run %d)', r));

        if numel(run_plan) == 1
            while true
                resp = lower(strtrim(input([ ...
                    '\nAre you happy with all templates?\n' ...
                    '  y = proceed to assigning all spikes\n' ...
                    '  g = refine in wave_clus GUI first\n' ...
                    '  q = quit now\n' ...
                    'Choice [y/g/q]: '], 's')));
                if isempty(resp), resp = 'g'; end

                if strcmp(resp, 'y')
                    break;
                elseif strcmp(resp, 'q')
                    fprintf('Stopped by user before assignment step.\n');
                    return;
                elseif strcmp(resp, 'g')
                    print_gui_refine_instructions(channel_dir);
                    if exist('wave_clus', 'file') == 2
                        launch = ask_yes_no('Open wave_clus GUI now? [y/n]: ', true);
                        if launch
                            wave_clus;
                        end
                    else
                        fprintf('wave_clus function not found on path. Open GUI manually after adding wave_clus path.\n');
                    end
                    input('After you refine and save in GUI, press Enter to continue...', 's');
                else
                    fprintf('Invalid choice. Please type y, g, or q.\n');
                end
            end
        else
            fprintf('Compare mode: proceeding directly to assignment for this run.\n');
        end

        fprintf('\n=== Step 2: Assign Remaining Clean + Dirty Spikes ===\n');
        t_assign_start = now;
        waveclus_opto_artifact_free(master_mat, ...
            'mode', 'assign', ...
            'channels', channels, ...
            'time_reference_mode', 'index_ms+firstAD', ...
            'assign_unlabeled_clean', true, ...
            'clean_distance_threshold_sd', clean_thr, ...
            'distance_threshold_sd', dirty_thr, ...
            'save_assignment_png', true, ...
            'overwrite_main_times', true, ...
            'debug', opts.debug);

        fprintf('\nDisplaying PNGs generated in step 2...\n');
        show_recent_pngs(channel_dir, t_assign_start, sprintf('Step 2 PNGs (run %d)', r));

        if numel(run_plan) > 1
            run_out_dir = fullfile(compare_root, run_tag);
            archive_recent_outputs(channel_dir, run_start, run_out_dir);
            fprintf('Saved run outputs to: %s\n', run_out_dir);
        end
    end

    fprintf('\nPipeline complete.\n');
    fprintf('Master: %s\n', master_mat);
    fprintf('Channel dir: %s\n', channel_dir);
    if numel(run_plan) > 1
        fprintf('Compare outputs: %s\n', compare_root);
    end
end

function master_mat = run_extract_and_get_master(plx_file, opts)
    [plx_dir, plx_name, ~] = fileparts(plx_file);

    try
        if isempty(opts.sdk_path)
            extract_and_filter(plx_file);
        else
            extract_and_filter(plx_file, 'sdk_path', opts.sdk_path);
        end
    catch ME
        fprintf('\nInitial extract_and_filter failed:\n%s\n', ME.message);
        fprintf('Please choose Plexon Offline SDK folder to retry.\n');
        sdk_dir = uigetdir(plx_dir, 'Select Plexon Offline SDK folder');
        if isequal(sdk_dir, 0)
            rethrow(ME);
        end
        extract_and_filter(plx_file, 'sdk_path', sdk_dir);
    end

    master_mat = fullfile(plx_dir, [plx_name '.mat']);
    if ~exist(master_mat, 'file')
        error('Expected master MAT not found after extraction: %s', master_mat);
    end
end

function [tf, missing] = is_master_file(mat_file)
    required = {'Eventstime', 'EventTag', 'nChannels', 'FsPlexon'};
    info = whos('-file', mat_file);
    names = {info.name};
    missing = required(~ismember(required, names));
    tf = isempty(missing);
end

function [channel_dir, channels] = discover_filtered_channels(master_mat)
    [master_dir, master_name, ~] = fileparts(master_mat);
    c1 = fullfile(master_dir, [master_name '_channels']);
    c2 = fullfile(master_dir, master_name);

    if exist(c1, 'dir')
        channel_dir = c1;
    elseif exist(c2, 'dir')
        channel_dir = c2;
    else
        error('No channel directory found near %s', master_mat);
    end

    files = dir(fullfile(channel_dir, 'channel_*_filtered_CAR.mat'));
    if isempty(files)
        files = dir(fullfile(channel_dir, 'channel_*_filtered_CARfiltered.mat'));
    end

    channels = zeros(0, 1);
    for i = 1:numel(files)
        tk = regexp(files(i).name, '^channel_(\d+)_', 'tokens', 'once');
        if ~isempty(tk)
            channels(end+1, 1) = str2double(tk{1}); %#ok<AGROW>
        end
    end
    channels = unique(sort(channels(:)'));
end

function txt = format_channel_compact(channels)
    channels = unique(sort(channels(:)'));
    if isempty(channels)
        txt = '(none)';
        return;
    end
    if numel(channels) > 1 && all(diff(channels) == 1)
        txt = sprintf('%d-%d', channels(1), channels(end));
        return;
    end
    txt = mat2str(channels);
end

function channels = prompt_channels(all_channels)
    fprintf('\nChannel selection:\n');
    fprintf('  Available channels: %s\n', format_channel_compact(all_channels));
    fprintf('  Press Enter for all channels.\n');
    fprintf('  Example inputs: 1:4   or   [1 3 5 7]\n');
    raw = strtrim(input('Enter channels: ', 's'));

    if isempty(raw)
        channels = all_channels;
        return;
    end

    parsed = str2num(raw); %#ok<ST2NM>
    if isempty(parsed) || ~isnumeric(parsed)
        error('Could not parse channel input. Use syntax like 1:4 or [1 3 5].');
    end
    channels = unique(round(parsed(:)'));
    bad = channels(~ismember(channels, all_channels));
    if ~isempty(bad)
        error('Requested channels not found in filtered files: %s', mat2str(bad));
    end
end

function choice = prompt_template_strictness()
    fprintf('\nTemplate-generation strictness (default = 2):\n');
    fprintf('  1) Lenient  : stdmin=4.0, min_clus=20, force_auto=true,  template_sdnum=3.0\n');
    fprintf('  2) Default  : stdmin=4.5, min_clus=20, force_auto=false, template_sdnum=2.5\n');
    fprintf('  3) Strict   : stdmin=5.0, min_clus=20, force_auto=false, template_sdnum=2.0\n');
    fprintf('  4) Compare  : run full pipeline on 1/2/3 template levels\n');

    raw = strtrim(input('Choose [1/2/3/4] (Enter=2): ', 's'));
    if isempty(raw), raw = '2'; end
    choice = str2double(raw);
    if ~ismember(choice, [1 2 3 4])
        error('Invalid template strictness choice.');
    end
end

function choice = prompt_assignment_strictness()
    fprintf('\nAssignment strictness (default = 2):\n');
    fprintf('  1) Lenient  : clean/dirty threshold = 3.0 SD\n');
    fprintf('  2) Default  : clean/dirty threshold = 2.5 SD\n');
    fprintf('  3) Strict   : clean/dirty threshold = 2.0 SD\n');
    fprintf('  4) Compare  : run full pipeline on 1/2/3 assignment levels\n');

    raw = strtrim(input('Choose [1/2/3/4] (Enter=2): ', 's'));
    if isempty(raw), raw = '2'; end
    choice = str2double(raw);
    if ~ismember(choice, [1 2 3 4])
        error('Invalid assignment strictness choice.');
    end
end

function [par, label] = build_template_params(mode)
    if exist('set_parameters_custom', 'file') == 2
        par = set_parameters_custom();
    elseif exist('set_parameters', 'file') == 2
        par = set_parameters();
    else
        error('No set_parameters_custom or set_parameters found on path.');
    end

    % Keep min_clus fixed for all strictness levels.
    par.min_clus = 20;

    switch mode
        case 1
            par.stdmin = 4.0;
            par.force_auto = true;
            par.template_sdnum = 3.0;
            label = 'Lenient';
        case 2
            % Use defaults as defined by set_parameters_custom/set_parameters.
            label = 'Default';
        case 3
            par.stdmin = 5.0;
            par.force_auto = false;
            par.template_sdnum = 2.0;
            label = 'Strict';
        otherwise
            error('Unknown template strictness mode.');
    end
end

function [clean_thr, dirty_thr, label] = map_assignment_thresholds(mode)
    switch mode
        case 1
            clean_thr = 3.0;
            dirty_thr = 3.0;
            label = 'Lenient';
        case 2
            clean_thr = 2.5;
            dirty_thr = 2.5;
            label = 'Default';
        case 3
            clean_thr = 2.0;
            dirty_thr = 2.0;
            label = 'Strict';
        otherwise
            error('Unknown assignment strictness mode.');
    end
end

function run_plan = build_run_plan(template_choice, assign_choice)
    if template_choice == 4 && assign_choice == 4
        % Compare both dimensions as matched levels: (1,1), (2,2), (3,3).
        run_plan = struct('template_mode', {1, 2, 3}, 'assign_mode', {1, 2, 3});
        return;
    end

    if template_choice == 4
        run_plan = struct('template_mode', {1, 2, 3}, ...
            'assign_mode', {assign_choice, assign_choice, assign_choice});
        return;
    end

    if assign_choice == 4
        run_plan = struct('template_mode', {template_choice, template_choice, template_choice}, ...
            'assign_mode', {1, 2, 3});
        return;
    end

    run_plan = struct('template_mode', template_choice, 'assign_mode', assign_choice);
end

function tf = ask_yes_no(prompt_text, default_yes)
    resp = lower(strtrim(input(prompt_text, 's')));
    if isempty(resp)
        tf = default_yes;
        return;
    end
    tf = strcmp(resp, 'y') || strcmp(resp, 'yes');
end

function print_gui_refine_instructions(channel_dir)
    fprintf('\nRefine templates in wave_clus GUI (brief):\n');
    fprintf('  1) Run: wave_clus\n');
    fprintf('  2) Load: channel_<N>_filtered_CAR.mat from:\n');
    fprintf('     %s\n', channel_dir);
    fprintf('  3) Adjust clusters and save.\n');
end

function show_recent_pngs(channel_dir, since_datenum, header_text)
    all_png = dir(fullfile(channel_dir, '**', '*.png'));
    if isempty(all_png)
        fprintf('%s: no PNG files found in %s\n', header_text, channel_dir);
        return;
    end

    dt_eps = 2 / 86400; % 2 seconds tolerance
    recent = all_png([all_png.datenum] >= (since_datenum - dt_eps));
    if isempty(recent)
        fprintf('%s: no new PNG files detected since this step started.\n', header_text);
        return;
    end

    [~, ord] = sort([recent.datenum], 'ascend');
    recent = recent(ord);

    fprintf('%s (%d files):\n', header_text, numel(recent));
    for i = 1:numel(recent)
        png_file = fullfile(recent(i).folder, recent(i).name);
        fprintf('  %s\n', png_file);
        show_png_file(png_file);
    end
end

function show_png_file(png_file)
    try
        img = imread(png_file);
        fig = figure('Name', png_file, 'NumberTitle', 'off', 'Color', 'w');
        imshow(img, 'Parent', axes('Parent', fig), 'InitialMagnification', 'fit');
        title(strrep(png_file, '_', '\_'), 'Interpreter', 'none', 'FontSize', 9);
        drawnow;
    catch
        try
            open(png_file);
        catch
            fprintf('  Could not display: %s\n', png_file);
        end
    end
end

function archive_recent_outputs(channel_dir, since_datenum, out_dir)
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    pats = {'**/*.png', '**/times*.mat', '**/*_dirty_assignment.mat'};
    files = [];
    for i = 1:numel(pats)
        d = dir(fullfile(channel_dir, pats{i}));
        if ~isempty(d)
            files = [files; d(~[d.isdir])]; %#ok<AGROW>
        end
    end
    if isempty(files)
        fprintf('No files found to archive for this run.\n');
        return;
    end

    dt_eps = 2 / 86400; % 2 seconds tolerance
    files = files([files.datenum] >= (since_datenum - dt_eps));
    if isempty(files)
        fprintf('No new files detected to archive for this run.\n');
        return;
    end

    for i = 1:numel(files)
        src = fullfile(files(i).folder, files(i).name);
        rel = strrep(src, [channel_dir filesep], '');
        dst = fullfile(out_dir, rel);
        dst_dir = fileparts(dst);
        if ~exist(dst_dir, 'dir')
            mkdir(dst_dir);
        end
        copyfile(src, dst);
    end
end
