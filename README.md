# Alternative_virmen_alignment

MATLAB pipeline for aligning multimodal neuroscience experiment data: Digidata/WaveSurfer sync signals, two-photon imaging frames, ViRMEn behavioral events, speaker stimuli, and ball-velocity sensors.

This repository is a research data-alignment work sample. It is **not** a turnkey production system and is **not** fully reproducible without institutional raw/processed data, mapped server drives, and a few external MATLAB dependencies.

---

## Start here

**Canonical entrypoint:** [`run_multimodal_alignment.m`](run_multimodal_alignment.m)

1. Edit the configuration blocks at the top of `run_multimodal_alignment.m` (session identity, channels, thresholds, run flags, calibration path, code folder).
2. Ensure required inputs exist on disk (see [Expected inputs](#expected-inputs)).
3. Run `run_multimodal_alignment` from MATLAB.
4. Optionally inspect saved outputs under `ProcessedData/<mouse>/<date>/`.

Supporting canonical scripts:

| Script | Role |
|--------|------|
| `get_frame_times_imaging.m` | Detect imaging frame times from the slow-galvo sync channel; compare counts to TIFF frames |
| `run_velocity_alignment.m` | Extract, calibrate, validate, and save frame-aligned ball velocity |
| `align_virmen_data.m` | Join ViRMEn behavior, stimuli, and neural activity into a trialized `imaging` structure |

Older drivers and alternate implementations live under [`legacy/`](legacy/) and are **not** part of the supported path.

---

## What this repository demonstrates

For a Data Ingestion Coordinator / research-data role, this codebase shows:

- **Ingestion of heterogeneous streams** — ABF/H5 sync files, TIFF-backed acquisitions, ViRMEn `.mat` behavior, preprocessed fluorescence (`dff` / `deconv`)
- **Synchronization across clocks** — map galvo frame times and ViRMEn iteration pulses onto a shared Digidata/WaveSurfer timebase
- **Validation checks** — acquisition-count matching, frame-count QC, channel-conflict checks, velocity shape/frame checks, trial “good_trial” marking
- **Harmonization** — classify speaker events into experimental conditions; project events into imaging frame indices; produce analysis-ready trial structures
- **Configurable session parameters** — paths, channels, sync filename filters, detection thresholds, and run flags grouped at the top of the wrapper

---

## High-level data flow

```
Raw sync (.abf / .h5) + imaging TSeries folders
        │
        ▼
get_frame_times_imaging  →  alignment_info.mat
        │
        ├──► run_velocity_alignment  →  corrected_velocity.mat
        │
        ├──► speaker detection (find_spkr_output_task_simple)
        │         + ViRMEn iteration decode
        │         + align_virmen_data
        │              → imaging.mat, alignment_variables.mat, vr_sound_frames*.mat
        │
        └──► passive speaker detection + binarize_passive_sounds
                  + align_passive_imagingst_updated_noise
                       → passive_frames.mat, imaging_st.mat
```

In short: **raw sync and imaging files → frame detection → stimulus and ViRMEn event alignment → frame-aligned imaging structure → calibrated ball velocity → saved outputs.**

---

## Configuration (edit in `run_multimodal_alignment.m`)

### Session

- `config.mouse`, `config.date`, `config.server`, `config.experimenter`
- `config.virmen_file_base` — ViRMEn filename stem (without `_Cell` / `_1` / `_2` suffixes)
- `config.condition_file` — `conditions_per_speaker` lookup `.mat`
- `config.calibration_file` — ball-sensor `calibration_info.mat`
- `config.code_folder` — local clone of this repository (added to the MATLAB path; `legacy/` is then removed from the path)

### Acquisition / task / passive

- `acquisition.galvo_channel`, `acquisition.virmen_channel`
- `acquisition.vr_sync_string`, `acquisition.passive_sync_string` — substrings used to select sync files
- `speaker.channel_numbers`, `speaker.ids`, `speaker.multiple_speakers_per_trial`
- `task.*` — stim flag and ITI / minimum trial timing
- `passive_align_info.before_frames`, `passive_align_info.after_frames`

### Sound detection

- `task_sound.*` / `passive_sound.*` — distances, duration, smoothing, detection thresholds, ITI-tone version, optional per-file overrides

### Run flags (`run_options`)

| Flag | Effect |
|------|--------|
| `process_task` | Run VR/task alignment branch |
| `process_passive` | Run passive alignment branch |
| `recompute_alignment` | Recompute `alignment_info` even if a cached file exists |
| `plot_task_alignment_qc` | Plot random aligned trials |
| `interactive_qc` | Enable interactive `pause` inside speaker detection |
| `calculate_velocity` | Run velocity extraction/calibration |
| `plot_velocity_qc` | Plot corrected pitch/roll/yaw |

Task and passive branches are independent: set either flag to `false` to skip that branch. Velocity runs when `calculate_velocity` is `true` (after frame alignment, before task/passive processing).

---

## Expected inputs

Institutional data are **not** included in this repository. Typical layout assumed by the wrapper:

```
<server>/<experimenter>/RawData/<mouse>/
    wavesurfer/<date>/          % *.abf or *.h5 sync files (names contain VR / passive)
    <date>/                     % *TSeries* folders with *Ch2* TIFFs
    virmen/<virmen_file_base>*.mat

<server>/<experimenter>/ProcessedData/<mouse>/<date>/
    dff.mat
    deconv/deconv.mat

config.condition_file           % must contain conditions_per_speaker
config.calibration_file         % must contain calibration_info (if calculate_velocity)
```

ViRMEn loading accepts one of: `*_2.mat` + `*_Cell_2.mat`, `*_1.mat` + `*_Cell_1.mat`, or `.mat` + `*_Cell.mat`.

Passive imaging alignment additionally expects (under the processed session folder): `z_dff.mat` and `corrected_velocity.mat` (the latter is produced when `calculate_velocity` is enabled).

---

## Generated outputs

Under `ProcessedData/<mouse>/<date>/`:

| Output | Contents |
|--------|----------|
| `alignment_info.mat` | Per-acquisition `frame_times`, sync IDs, sampling rate, frame-count QC fields |
| `corrected_velocity.mat` | `raw_velocity`, `corrected_velocity`, `calibration_info`, velocity `qc` |

Under `.../VR/`:

| Output | Contents |
|--------|----------|
| `alignment_variables.mat` | Session `info`, `sound_info`, `task_info`, trial/sound intermediate structures |
| `imaging.mat` | Trialized neural + behavior structure |
| `vr_sound_frames.mat` / `vr_sound_frames_updated.mat` | Stimulus frames; updated version restricted to good imaging trials |

Under `.../passive/`:

| Output | Contents |
|--------|----------|
| `alignment_variables.mat` | Passive sound intermediates |
| `passive_frames.mat` | Passive stimulus frame indices |
| `imaging_st.mat` | Passive trial-aligned snippets |

---

## Setup and usage

1. Clone this repository and set `config.code_folder` to that local path.
2. Map or mount institutional servers so `config.server` and related paths resolve.
3. Install/configure required MATLAB toolboxes and external helpers (below).
4. Point `config.condition_file` and `config.calibration_file` at your lookup/calibration mats.
5. Set channels and sync strings to match the recording.
6. Choose `run_options` (task / passive / velocity / plots).
7. Run `run_multimodal_alignment`.

The wrapper validates required folders/files and channel conflicts before loading data. It does **not** use `uigetdir` for path selection.

---

## MATLAB toolboxes and external dependencies

**In-repo**

- `abfload.m` — Axon Binary File loader (third-party utility included here)

**MATLAB / toolboxes used by the canonical path**

- Base MATLAB
- Signal Processing Toolbox — `findpeaks` (frame detection)
- Statistics and Machine Learning Toolbox — `zscore` in `align_virmen_data`

**External (not in this repository)**

- WaveSurfer MATLAB API — `ws.loadDataFile` when sync files are `.h5` instead of `.abf`
- Preprocessed imaging products (`dff`, `deconv`, and for passive `z_dff`) from the lab imaging pipeline
- Lab-specific `condition_per_speaker.mat` and `calibration_info.mat`

---

## Quality control (what actually exists)

Honest summary of checks in the canonical path:

- **Acquisition pairing** — number of `*TSeries*` folders must equal number of sync files; pairings are printed
- **Frame-count validation** — detected/filled galvo peaks compared to `*Ch2*` TIFF counts; mismatch warns; extras truncated to TIFF count
- **Photostimulation gaps** — long galvo intervals filled at the estimated frame rate; filled indices stored in `bad_frames` / QC metadata
- **Missing-frame warnings** — warning when final frame times &lt; expected TIFF count
- **Config validation** — required paths/files; galvo / ViRMEn / speaker channel conflicts; speaker ID count match
- **Velocity checks** — expects 3×N array (pitch, roll, yaw); warns if length ≠ total imaging frames; counts NaN frames
- **Trial checks** — `align_virmen_data` marks in-bounds trials with `good_trial = 1` and clears the flag for “weird” trials; `fix_vr_sound_frames` keeps sounds only for good trials
- **Optional visual QC** — galvo peak plots, random trial overlays, velocity traces; speaker-detection `pause` only if `interactive_qc` is true

These checks improve confidence but do **not** guarantee a perfect alignment. Some heuristic recoveries (for example unpaired speaker events or inferred conditions) can still produce saved outputs that need human review.

---

## Canonical vs legacy

| Location | Status |
|----------|--------|
| `run_multimodal_alignment.m`, `get_frame_times_imaging.m`, `run_velocity_alignment.m`, `align_virmen_data.m`, and helpers they call at the repo root / `align_sounds/` | **Canonical** |
| `legacy/` (including `legacy/old_code/` and `legacy/alternative_code_for_missing_iterations/`) | **Legacy** — older drivers, alternate matchers, and recovery variants |

The wrapper adds `config.code_folder` with `genpath`, then **removes** `legacy/` from the MATLAB path so legacy duplicates cannot shadow canonical functions.

Some root-level `.m` files (for example `match_trialsperfile*.m`, `shift_sync_data.m`, `find_reward*.m`, `redo_imaging.m`) are retained from earlier workflows but are **not** called by `run_multimodal_alignment.m`. Treat them as non-entrypoint helpers unless you are extending the pipeline.

---

## Limitations and known assumptions

- Institutional raw/processed data and server drive letters are required; example paths in the config block are lab-specific placeholders.
- Assumes two imaging channels and counts `*Ch2*` TIFFs for expected frame counts.
- Assumes TSeries folders and sync files sort into matching order.
- Reward-event detection is **not** enabled in the canonical wrapper (`reward_loc_pure_frames` is empty).
- Several helpers still use `cd` internally when reading sync directories; the main wrapper itself does not.
- `binarize_passive_sounds` still inserts a short `pause(2)` during plotting.
- Outputs overwrite same-named `.mat` files in the processed session folders on re-run.
- No unit test suite or CI is included.

---

## Repository structure

```
.
├── run_multimodal_alignment.m      # canonical entrypoint
├── get_frame_times_imaging.m       # frame-clock detection + QC
├── run_velocity_alignment.m        # velocity extract / calibrate / save
├── get_velocity.m
├── correct_velocity.m
├── align_virmen_data.m             # trialized multimodal join
├── abfload.m
├── align_sounds/                   # speaker detection + passive alignment helpers
├── legacy/                         # older / alternate workflows (not on path)
└── README.md
```

---

## License / data note

Large research datasets are not distributed with this repository. Contact the author for access policies if you need example data for evaluation.
