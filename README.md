# Alternative_virmen_alignment

A MATLAB-based pipeline for synchronizing behavioral, neural imaging, and stimulus data streams collected during virtual reality neuroscience experiments.

## Overview

Large-scale neuroscience experiments often generate multiple independent data streams, including behavioral events, neural imaging recordings, audio stimulus outputs, and hardware synchronization signals. Small timing errors, missing events, or acquisition inconsistencies can prevent reliable downstream analysis.

This project implements a robust alignment and quality-control framework that detects stimulus events, reconstructs behavioral timing, reconciles independent acquisition clocks, and maps experimental events directly into imaging frame space. The pipeline was developed to process and validate data across dozens of experiments containing thousands of behavioral trials.

## Key Features

### Signal Processing and Event Detection

* Processes raw speaker output channels to identify stimulus events.
* Applies baseline correction, thresholding, nonlinear compression, and envelope-based processing to improve event detection robustness.
* Detects stimulus onset and offset times while rejecting noise and spurious events.

### Trial Reconstruction and Classification

* Groups detected events into behavioral trials using timing constraints.
* Classifies trials according to experimental conditions and active stimulus channels.
* Handles incomplete or truncated recordings through automated recovery logic.

### Multimodal Data Alignment

* Reconstructs behavioral iteration timing from synchronization pulse trains.
* Aligns behavioral, stimulus, and acquisition systems with independent clocks.
* Estimates and corrects timing offsets between data streams.
* Projects events into neural imaging frame coordinates for downstream analysis.

### Quality Control and Validation

* Performs automated consistency checks on event ordering, timing, durations, and alignment accuracy.
* Identifies synchronization failures, anomalous trials, and recording irregularities.
* Generates visualization tools for rapid human validation and troubleshooting.
* Includes handling for common real-world acquisition failures and edge cases.

## Applications

This pipeline was used to generate analysis-ready datasets from multimodal neuroscience experiments, enabling reliable integration of behavioral measurements, sensory stimuli, and neural population recordings for downstream statistical modeling and machine learning analyses.

## Technologies

* MATLAB
* Signal Processing
* Time-Series Analysis
* Data Validation and Quality Control
* Multimodal Data Integration
* Experimental Data Synchronization

