# Hammer Neuromodulation Lab (HNL)

This repository contains MATLAB scripts used for processing data from **Medtronic Percept™ PC Neurostimulators** with BrainSense™ technology, as part of research conducted in the **Hammer Neuromodulation Lab** on patients with **Parkinson's Disease (PD)**.

## Overview

These tools were developed to streamline the organization and processing of JSON data recorded from deep brain stimulation (DBS) devices. The scripts help researchers efficiently parse, label, and structure patient neural data for analysis.

---

## Scripts

### `Chan_LoadJSON.m`

Organizes raw JSON files exported from Medtronic's Percept device:

* Automatically renames files using a standardized format:
  `PatientLastName_Date_Diagnosis.json`
* Sorts and moves each file into a structured directory system based on patient identity.
* Outputs organized files into a selected Box project folder for further analysis or archival.

### `rawsignalandstimulationintensity2.m`

Processes a single JSON file from a Percept recording session and outputs a clean MATLAB table with:

| Time (s) | Left LFP Signal | Right LFP Signal | Left Stimulation Intensity | Right Stimulation Intensity |
| -------- | --------------- | ---------------- | -------------------------- | --------------------------- |

* Synchronizes signal and stimulation values.
* Useful for time-series visualization, feature extraction, or integration into larger analyses.
* Provides the processed data for signalvsstimuation scripts.

### `signalvsstimulation_part1.m`

Coregisters high-rate LFP signal data with low-rate stimulation intensity and plots them together:

* Reads a JSON file and decodes time-domain (TD) and LFP streams.
* Builds precise timestamps for each TD sample and LFP packet event.
* Interpolates stimulation current (mA) onto each TD timestamp using previous-value interpolation.
* Outputs a MATLAB table `final` with columns:
  **Time**, **LeftRaw**, **LeftStim\_mA**, **RightRaw**, **RightStim\_mA**
* Generates a two-panel figure overlaying raw LFP signals (blue) and stimulation intensity (red) on dual y-axes for left and right channels.

### `signalvsstimulation_part2.m`

Extends the functionality of `signalvsstimulation_part1.m` by **coregistering the full power spectrum of raw LFP signals with stimulation intensity**, and calculating **gamma band power (60–90 Hz)** across time-segmented epochs.

**Key Features:**

* Reads and decodes a Percept JSON file to extract time-domain (TD) and LFP data streams.
* Synchronizes raw LFP signals with corresponding stimulation intensity using timestamp interpolation and packet sequence tracking.
* Segments the data into **10-second, non-overlapping windows** and checks each segment for time continuity (skips segments with dropped packets).
* Uses the **Welch power spectral density (PSD) estimation** method to compute power spectra from raw LFP signals.
* Calculates **average stimulation intensity** and **total gamma power** for both left and right hemispheres within each segment.
* Outputs a new table, `gammaData`, with the following columns:

| Segment | FreqL | PowerL | AvgStimL | FreqR | PowerR | AvgStimR |
|---------|-------|--------|----------|-------|--------|----------|

* Plots **gamma band power** and **average stimulation intensity** side-by-side across segments for both hemispheres, allowing visual correlation between stimulation level and neural oscillatory activity.

This script enables researchers to explore how adaptive stimulation may affect oscillatory activity in the gamma band—a signal often implicated in therapeutic responses in Parkinson’s Disease.

### `droppedpacketdetectionworking.m`

Checks for dropped data packets that may have occurred during syncing between the DBS implant and the tablet:

* Identifies missing or irregular time intervals in the data stream.
* Flags potential communication errors that could affect data integrity.
* Helps ensure accurate temporal alignment for downstream analyses.

---

## Requirements

* MATLAB R2021a or later
* JSON files exported directly from Medtronic Percept PC
* [JSONlab Toolbox](https://github.com/fangq/jsonlab) (if needed for parsing)

---

## Usage

1. Clone the repository or download the scripts.
2. Run `Chan_LoadJSON.m` to organize and label your patient JSON files.
3. Use `rawsignalandstimulationintensity2.m` to convert JSON data into MATLAB tables.
4. Execute `signalvsstimulation.m` to visualize coregistered LFP signal and stimulation intensity.
5. Run `droppedpacketdetectionworking.m` to check for time-series discontinuities or packet loss in the raw data.

---

## Author

**Andrew Chan**
Bioengineering, University of Pittsburgh
Researcher at the Hammer Neuromodulation Lab
