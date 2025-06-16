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

### `realtimestampingv2.m`

Generates accurate real-time timestamps for both Time Domain and LFP data in BrainSense JSON files:

* Anchors signals to UTC using the first packet timestamp and reconstructs sample-level timestamps based on sampling frequency and packet size.
* Detects and flags dropped packets and timing misalignments due to sequence number or tick discrepancies.
* Supports multiple LFP data structures (`BrainSenseLfp`, `LFPData`, `LFPmontageTimeDomain`) and assigns timestamp arrays back to MATLAB workspace for analysis.
* Ensures precise alignment across neural signal streams for time-resolved analyses.

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

### `addingstimrate_part1.m`

Extracts and aligns the stimulation rate from LFP data to each time‐domain sample:

* Prompts for a BrainSense JSON file and decodes its contents.
* Converts the raw `BrainSenseTimeDomain` struct array into a MATLAB table.
* Adds a new `StimRateHz` column, initialized to `NaN`.
* Parses all packet timestamps into `datetime` objects with UTC timezone.
* For each time‐domain row, finds the matching LFP packet within ±1.5 s.
* Splits combined LFP channel pairs (e.g. `'ZERO_TWO_LEFT,ZERO_TWO_RIGHT'`) and selects the side matching the time‐domain channel.
* Extracts the corresponding `RateInHertz` from `TherapySnapshot.Left` or `.Right` and assigns it to `StimRateHz`.
* Saves the updated `BrainSenseTimeDomain` table back into the base workspace for downstream analysis.

**Usage:**  
Run `addingstimrate_part1` in MATLAB to automatically load your JSON file and annotate `BrainSenseTimeDomain` with stimulation rates.

---

## Requirements

* MATLAB R2021a or later
* JSON files exported directly from Medtronic Percept PC

---

## Author

**Andrew Chan**

Bioengineering, University of Pittsburgh

Researcher at the Hammer Neuromodulation Lab, University of Pennsylvania
