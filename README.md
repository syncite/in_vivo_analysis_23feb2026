# In Vivo Optogenetic Spike Sorting (Artifact-Free with wave_clus)

This is a step-by-step protocol for running your full analysis pipeline, including:
- PLX extraction + filtering
- clean-only clustering for wave_clus
- manual curation
- LED-period spike assignment
- downstream peri-event and waveform analysis

Scripts used:
- `extract_and_filter.m`
- `waveclus_opto_artifact_free.m`
- `run_opto_sorting_pipeline.m` (new one-command interactive entry point)
- `analyse_peri_event.m`
- `visualise_waveform.m`

## Quick Start (Single Command)

If you want one guided command in MATLAB:

```matlab
run_opto_sorting_pipeline
```

What it does interactively:
- asks you to pick `.plx` (raw) or master `.mat` (already filtered)
- if `.plx`, runs extraction/filtering automatically
- confirms detected filtered channel files and asks whether to proceed
- asks template-generation strictness (with compare option)
- asks assignment strictness (with compare option)
- asks channel selection (default all; supports input like `1:4`)
- runs clean-template generation (step 1)
- displays generated PNGs
- asks whether to refine in wave_clus GUI or proceed
- runs reassignment of remaining clean + dirty spikes (step 2)
- displays resulting PNGs

## 1) Start MATLAB and set paths

```matlab
cd('/Users/wjr/in_vivo_analysis_23feb2026');
addpath('/Users/wjr/in_vivo_analysis_23feb2026');
addpath(genpath('/ABSOLUTE/PATH/TO/wave_clus'));
```

## 2) Run PLX extraction + filtering

```matlab
plx = '/ABSOLUTE/PATH/TO/your_recording.plx';

extract_and_filter(plx, ...
    'sdk_path', '/ABSOLUTE/PATH/TO/Plexon_Offline_SDK_Bundle');
```

This creates:
- master file: `<plxname>.mat`
- channel folder: `<plxname>_channels/`
- per-channel input files for wave_clus: `channel_<N>_filtered_CAR.mat`

## 3) Define the master file path

```matlab
[plxDir, plxName, ~] = fileparts(plx);
master = fullfile(plxDir, [plxName '.mat']);
```

## 4) Prepare clean-only clustering input

This detects spikes on full data, splits into clean vs LED-contaminated windows, and clusters using clean spikes only.
The original channel input `.mat` files are not overwritten.

```matlab
waveclus_opto_artifact_free(master, 'mode', 'prepare');
```

Per channel, it writes:
- `channel_<N>_filtered_CAR_spikes.mat` (full spikes)
- `channel_<N>_filtered_CAR_clean_spikes.mat` (clean-only spikes)
- `channel_<N>_filtered_CAR_dirty_spikes.mat` (dirty-only spikes)
- `channel_<N>_filtered_CAR_led_split.mat` (metadata + masks)

Masking rule:
- from LED onset to `100 ms after LED offset`
- default is `1.0 s + 0.1 s = 1.1 s` per LED event
- default time reference is `index_ms+firstAD` (consistent across channels)

## 5) Review clusters first, then curate only if needed

After step 4, first review the wave_clus output figures (for example `fig2print` files in the channel folder).
If clusters already look good, you can skip GUI edits and continue to step 6.

If you want manual edits:

```matlab
wave_clus
```

In the GUI:
- Load each `channel_<N>_filtered_CAR.mat` from `<plxname>_channels`
- Curate clusters as usual
- Save results

Why this is clean-only:
- clustering was run from `*_clean_spikes.mat`
- that clean result is copied to canonical `times_channel_<N>_filtered_CAR.mat` for GUI compatibility

## 6) Reassign remaining clean + dirty spikes to curated templates

By default, step 6 will:
- reassign clean spikes still in cluster 0 using curated clean templates
- assign dirty spikes using the same template set
- leave spikes as 0 if they fail threshold
- save a step-2 PNG summary (`*_assignment_summary.png`)

These defaults are hard-coded in the script:
- time reference: `index_ms+firstAD`
- clean cluster-0 reassignment: ON
- clean threshold: `2.5 SD`
- dirty threshold: `2.5 SD`
- assignment PNG output: ON

```matlab
waveclus_opto_artifact_free(master, 'mode', 'assign');
```

With `overwrite_main_times = true`, the canonical `times_channel_<N>_filtered_CAR.mat` is replaced with combined clean+assigned output (and a clean-only backup is saved).

## 7) Optional second curation pass (after assignment)

```matlab
wave_clus
```

Open the same `channel_<N>_filtered_CAR.mat` files again and do final cleanup curation if needed.

## 8) Run downstream analysis scripts

```matlab
analyse_peri_event(master);
visualise_waveform(master);
```

## Important notes

- Run order is:
  1) `extract_and_filter`
  2) `waveclus_opto_artifact_free(..., 'mode', 'prepare')`
  3) optional wave_clus manual curation
  4) `waveclus_opto_artifact_free(..., 'mode', 'assign')`
  5) optional wave_clus recuration
  6) `analyse_peri_event` / `visualise_waveform`

- Do **not** repeatedly run `assign` on already combined files.
- If you need to rerun assignment, restore clean-only files first:
  - `times_channel_<N>_filtered_CAR_clean_only.mat`

## Optional custom LED settings

If your LED tags or pulse durations differ:

```matlab
waveclus_opto_artifact_free(master, 'mode', 'prepare', ...
    'led_tags', [32385 32321 32417], ...
    'led_duration_s', [1.0 0.1 0.01], ...
    'post_led_buffer_s', 0.1);
```
