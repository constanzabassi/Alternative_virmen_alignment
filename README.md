# Multimodal Imaging and Behavioral Alignment

A MATLAB workflow for synchronizing behavioral, neural imaging, stimulus, and acquisition data collected during virtual-reality neuroscience experiments.

Large-scale neuroscience experiments often generate multiple independent data streams, including behavioral events, two-photon imaging recordings, audio stimulus outputs, ball-motion signals, and hardware synchronization channels. Timing differences, missing events, dropped frames, and acquisition inconsistencies can prevent reliable downstream analysis.

This repository contains a research workflow developed to align these data streams in a common imaging-frame coordinate system and generate analysis-ready outputs for downstream statistical modeling and machine-learning analyses.

This is a research code sample rather than a turnkey production system. Full execution requires institutional data, mapped server locations, and external MATLAB dependencies that are not distributed with the repository.

## Start here

The canonical entrypoint is `run_multimodal_alignment.m`.

Edit the configuration block at the top of the wrapper to specify:

- session ID and date
- local data paths
- galvo and ViRMEn channels
- speaker channels and IDs
- synchronization strings
- sound-detection thresholds
- velocity calibration path
- task, passive, and plotting options

## High-level workflow

Raw synchronization and imaging files  
→ imaging-frame detection  
→ stimulus and ViRMEn event alignment  
→ frame-aligned neural and behavioral data  
→ calibrated ball velocity and saved outputs

Main functions:

- `get_frame_times_imaging.m` detects imaging-frame times and compares detected frames with TIFF counts.
- `align_virmen_data.m` aligns behavioral, stimulus, and neural data into imaging-frame coordinates.
- `run_velocity_alignment.m` extracts, calibrates, validates, and saves frame-aligned ball velocity.
- `run_multimodal_alignment.m` coordinates the full workflow.

## Expected inputs

The associated research data are not included in this repository.

Typical inputs include:

- Digidata `.abf` or WaveSurfer `.h5` synchronization files
- two-photon imaging TSeries folders containing TIFF files
- ViRMEn behavioral `.mat` files (virtual reality behavioral data)
- preprocessed `dff.mat` and `deconv.mat` files
- speaker-condition lookup data
- microscope-specific ball-calibration data

## Generated outputs

Depending on the enabled processing stages, the workflow produces:

- `alignment_info.mat`
- `corrected_velocity.mat`
- `VR/alignment_variables.mat`
- `VR/imaging.mat`
- `VR/vr_sound_frames.mat`
- `VR/vr_sound_frames_updated.mat`
- `passive/alignment_variables.mat`
- `passive/passive_frames.mat`
- `passive/imaging_st.mat`

## Quality control

The canonical workflow includes:

- validation that the number of TSeries folders matches the number of synchronization files
- printed imaging/synchronization file pairings for review
- filtering of detected galvo peaks by scanning amplitude
- comparison of detected frame counts with TIFF counts
- identification and reconstruction of photostimulation-related frame gaps
- warnings when fewer frame times than TIFFs remain
- saved expected, initial, final, and filled-frame counts
- validation of channel conflicts and speaker mappings
- validation of velocity shape and total frame count
- counts of frames containing missing velocity values
- optional plots for frame detection, frame intervals, trial alignment, and corrected velocity

Some thresholds and alignment decisions are experiment-specific and should be reviewed using the saved QC fields, warnings, and visualizations.

## Setup

1. Clone the repository.
2. Open MATLAB in the repository directory.
3. Edit the configuration values at the top of `run_multimodal_alignment.m`.
4. Confirm that required institutional drives and input files are available.
5. Configure channels, speaker mappings, synchronization strings, thresholds, and run options.
6. Run `run_multimodal_alignment`.
7. Review warnings, QC plots, and saved QC fields before downstream analysis.

## Dependencies

Included in this repository:

- `abfload.m` for reading ABF files

MATLAB requirements:

- base MATLAB
- Signal Processing Toolbox for `findpeaks`

External dependencies may include:

- WaveSurfer MATLAB API for `ws.loadDataFile`
- institutional raw and processed data
- speaker-condition lookup files
- microscope-specific ball-calibration files

## Canonical and legacy code

`run_multimodal_alignment.m` is the supported entrypoint.

The `legacy/` folder contains historical and alternative workflows used for sessions with unusual or incomplete synchronization signals. These files are retained for provenance but are not recommended for new runs.

## Known assumptions and limitations

- TSeries folders and synchronization files are assumed to sort into corresponding order.
- The frame-count workflow assumes two imaging channels and counts `*Ch2*` TIFF files.
- Detection thresholds may require adjustment for different acquisition configurations.
- Reward-event detection is not enabled in the canonical wrapper.
- Re-running the workflow may overwrite existing output files with the same names.
- Full reproduction requires institutional data and external dependencies.
