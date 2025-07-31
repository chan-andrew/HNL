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

### `addingstimrate.m`

Extracts and aligns the stimulation rate from LFP data to each time‐domain sample:

* Prompts for a BrainSense JSON file and decodes its contents.
* Converts the raw `BrainSenseTimeDomain` struct array into a MATLAB table.
* Adds a new `StimRateHz` column, initialized to `NaN`.
* Parses all packet timestamps into `datetime` objects with UTC timezone.
* For each time‐domain row, finds the matching LFP packet within ±1.5 s.
* Splits combined LFP channel pairs (e.g. `'ZERO_TWO_LEFT,ZERO_TWO_RIGHT'`) and selects the side matching the time‐domain channel.
* Extracts the corresponding `RateInHertz` from `TherapySnapshot.Left` or `.Right` and assigns it to `StimRateHz`.
* Saves the updated `BrainSenseTimeDomain` table back into the base workspace for downstream analysis.

### `threecolumns.m`

Adds a structured column of nested tables to `BrainSenseTimeDomain`, containing **DateTime**, **Stimulation Amplitude**, and a placeholder for **Gamma Power**:

* Prompts for a BrainSense JSON file and decodes its contents.
* Converts `BrainSenseTimeDomain` to a MATLAB table if not already formatted.
* Parses all timestamps (`FirstPacketDateTime`) into `datetime` objects with UTC timezone.
* For each time-domain row, finds the corresponding LFP packet within ±1.5 seconds.
* Confirms that the LFP packet contains a matching channel (e.g. `'ZERO_TWO_LEFT'`).
* Extracts the stimulation rate (`RateInHertz`) and adds it to a new `StimRateHz` column.
* Calculates **average stimulation amplitude** (in mA) from LFP samples and stores it per row.
* Constructs a new column `ThreePart` containing 1×3 sub-tables with:
  - `DateTime`
  - `StimAmp` (average amplitude)
  - `GammaPower` (placeholder set to `NaN` for now)
* Outputs the updated `BrainSenseTimeDomain` table to the base workspace for downstream analysis.

* ### `timedomain.m`

Analyzes how **gamma-band power** in raw time-domain neural signals correlates with **stimulation amplitude** from a Medtronic Percept JSON file.

**Main Functions:**

* Parses both `BrainSenseTimeDomain` and `BrainSenseLfp` fields from the JSON file.
* Converts packet timestamps to `datetime` with proper UTC alignment.
* Matches each time-domain row to the closest LFP packet within ±1.5 seconds.
* Extracts stimulation **rate** and **amplitude**, synchronizing them with raw signal timestamps.
* Allows user to select a specific `(Channel, StimRateHz)` combination for focused analysis.

**Signal & Stimulation Visualization:**

* Aggregates raw signal samples (sampled at 250 Hz) and stimulation amplitude (sampled at 2 Hz).
* Plots signal and stimulation on **dual y-axes**, with:
  * **Left y-axis**: Raw neural signal (blue)
  * **Right y-axis**: Stimulation amplitude in mA (red)
* Highlights dropped packets (gaps >50 ms) with shaded regions on the plot.
* Aligns the y-axis zero lines visually for better interpretation.
* Prints detailed statistics about recording duration, signal/stimulation range, and sample counts.

**Gamma Power Analysis:**

* Segments the signal into **10-second epochs** and checks each for dropped packets or incomplete data.
* Computes **power spectral density (PSD)** using Welch’s method.
* Extracts **gamma-band power (60–70 Hz)** and aligns it with **average stimulation amplitude** per epoch.
* Produces a scatter plot of gamma power vs. stimulation amplitude:
  * Adds a **linear regression trendline**
  * Reports **R² value**, **correlation coefficient**, and **slope**
* Optionally generates a **2D heatmap** of gamma power vs stim amplitude density (requires ≥20 valid epochs).
* Saves results to the MATLAB workspace as:
  * `gamma_stim_results` (table)
  * `epoch_results` (array)
 
* ### `onespectrumforallpatients.m`

Comprehensive batch analysis tool that processes multiple BrainSense JSON files and automatically generates power spectral density analysis for **all channel/stimulation rate combinations**:

* **Revolutionary selective JSON extraction**: Processes large files (up to 500MB) by extracting only necessary sections without loading entire files into memory.
* **Multi-file batch processing**: Handles multiple JSON files simultaneously with automatic metadata extraction (lead location: STN/GPi, hemisphere: Left/Right).
* **Intelligent data standardization**: Concatenates data across files with different structures by automatically adding missing columns with appropriate default values.
* **Temporal alignment**: Matches TimeDomain and LFP data using optimized nearest-neighbor search with configurable time thresholds.
* **Automatic channel analysis**: Identifies and processes every unique channel/stimulation rate combination without user intervention.
* **Spectral analysis**: Segments data into 10-second epochs using Welch's method (1-second windows, 50% overlap) and calculates mean power spectral density with confidence intervals.
* **Comprehensive visualization**: Generates multi-panel subplot grid showing all combinations with highlighted gamma bands (stimulation_rate/2 ± 1 Hz).
* **Memory optimized**: Uses chunked processing and immediate cleanup for efficient handling of large datasets.

**Outputs:**
* Multi-panel figure with spectral analysis for all channel/rate combinations
* `all_spectrum_results`: Detailed spectral data for each combination  
* `channelRateCombos`: List of all processed combinations
* Concatenated TimeDomain and LFP data tables

---

## Requirements

* MATLAB R2021a or later
* JSON files exported directly from Medtronic Percept PC

---

## Author

**Andrew Chan**

Bioengineering, University of Pittsburgh

Researcher at the Hammer Neuromodulation Lab, University of Pennsylvania
