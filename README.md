# Hammer Neuromodulation Lab (HNL)

This repository contains MATLAB scripts used for processing data from **Medtronic Percept™ PC Neurostimulators** with BrainSense™ technology, as part of research conducted in the **Hammer Neuromodulation Lab** on patients with **Parkinson's Disease (PD)**.

## Overview

These tools were developed to streamline the organization and processing of JSON data recorded from deep brain stimulation (DBS) devices. The scripts help researchers efficiently parse, label, and structure patient neural data for analysis.

---

## Scripts

### `Chan_LoadJSON.m`

Organizes raw JSON files exported from Medtronic's Percept device:

- Automatically renames files using a standardized format:  
  `PatientLastName_Date_Diagnosis.json`
- Sorts and moves each file into a structured directory system based on patient identity.
- Outputs organized files into a selected Box project folder for further analysis or archival.

### `rawsignalandstimulationintensity2.m`

Processes a single JSON file from a Percept recording session and outputs a clean MATLAB table with:

| Time (s) | Left LFP Signal | Right LFP Signal | Left Stimulation Intensity | Right Stimulation Intensity |
|----------|------------------|-------------------|-----------------------------|------------------------------|

- Synchronizes signal and stimulation values.
- Useful for time-series visualization, feature extraction, or integration into larger analyses.

### `droppedpacketdetectionworking.m`

Checks for dropped data packets that may have occurred during syncing between the DBS implant and the tablet:

- Identifies missing or irregular time intervals in the data stream.
- Flags potential communication errors that could affect data integrity.
- Helps ensure accurate temporal alignment for downstream analyses.

---

## Requirements

- MATLAB R2021a or later
- JSON files exported directly from Medtronic Percept PC
- [JSONlab Toolbox](https://github.com/fangq/jsonlab) (if needed for parsing)

---

## Usage

1. Clone the repository or download the scripts.
2. Run `Chan_LoadJSON.m` to organize and label your patient JSON files.
3. Use `rawsignalandstimulationintensity2.m` to convert JSON data into MATLAB tables.
4. Run `droppedpacketdetectionworking.m` to check for time-series discontinuities or packet loss in the raw data.

---

## Author

**Andrew Chan**  
Bioengineering, University of Pittsburgh  
Researcher at the Hammer Neuromodulation Lab
