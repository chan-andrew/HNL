# 🧠 Hammer Neuromodulation Lab (HNL)
## MATLAB Tools for Medtronic Percept™ PC Data Processing

> **Advanced MATLAB toolkit for processing neural data from Medtronic Percept™ PC Neurostimulators with BrainSense™ technology, specifically designed for Parkinson's Disease research.**

---

## 📋 Table of Contents
- [Overview](#overview)
- [Installation & Requirements](#installation--requirements)
- [Core Scripts](#core-scripts)
  - [Data Organization](#data-organization)
  - [Signal Processing](#signal-processing)
  - [Spectral Analysis](#spectral-analysis)
  - [Quality Control](#quality-control)
- [Data Output Formats](#data-output-formats)
- [Usage Examples](#usage-examples)
- [Contributing](#contributing)
- [Author](#author)

---

## 🎯 Overview

This repository provides a comprehensive suite of MATLAB scripts designed to streamline the processing and analysis of JSON data recorded from deep brain stimulation (DBS) devices. Our tools enable researchers to efficiently parse, organize, and analyze patient neural data with a focus on:

- **Data Organization**: Automated file management and standardized naming conventions
- **Signal Processing**: High-precision timestamp reconstruction and signal alignment
- **Spectral Analysis**: Advanced frequency domain analysis with stimulation correlation
- **Quality Assurance**: Comprehensive packet loss detection and data validation

---

## 🔧 Installation & Requirements

### System Requirements
- **MATLAB**: R2021a or later
- **Input Format**: JSON files exported directly from Medtronic Percept™ PC
- **Memory**: Minimum 8GB RAM recommended for large datasets
- **Storage**: Variable depending on dataset size

### Setup
1. Clone or download this repository
2. Add the script directory to your MATLAB path
3. Ensure JSON files are accessible from MATLAB workspace

---

## 📁 Core Scripts

### 🗂️ Data Organization

#### `Chan_LoadJSON.m`
**Purpose**: Automated organization of raw JSON exports

**Features**:
- Standardized file renaming: `PatientLastName_Date_Diagnosis.json`
- Intelligent directory structuring based on patient identity
- Integration with Box project folders for collaborative research
- Batch processing capabilities

---

### 🔄 Signal Processing

#### `rawsignalandstimulationintensity2.m`
**Purpose**: Single-file processing with synchronized output

**Output Table Structure**:
| Column | Description | Units |
|--------|-------------|-------|
| `Time` | Temporal reference | seconds |
| `Left LFP Signal` | Left hemisphere neural activity | µV |
| `Right LFP Signal` | Right hemisphere neural activity | µV |
| `Left Stimulation Intensity` | Left hemisphere stimulation | mA |
| `Right Stimulation Intensity` | Right hemisphere stimulation | mA |

**Key Capabilities**:
- Precise signal-stimulation synchronization
- Time-series visualization ready output
- Foundation for downstream analysis pipelines

#### `realtimestampingv2.m`
**Purpose**: High-precision timestamp reconstruction

**Advanced Features**:
- UTC-anchored timestamp generation using first packet reference
- Sample-level timestamp reconstruction based on sampling frequency
- Multi-format LFP data structure support:
  - `BrainSenseLfp`
  - `LFPData` 
  - `LFPmontageTimeDomain`
- Comprehensive packet validation and error flagging
- Automatic detection of sequence number discrepancies

#### `addingstimrate.m`
**Purpose**: Stimulation rate extraction and alignment

**Process Flow**:
1. JSON file parsing and structure conversion
2. Temporal packet matching (±1.5s tolerance)
3. Channel-specific rate extraction from therapy snapshots
4. Integration into existing data structures

---

### 📊 Spectral Analysis

#### `signalvsstimulation_part1.m`
**Purpose**: Basic signal-stimulation coregistration with visualization

**Output Specifications**:
- **Data Table**: `final` with synchronized timestamps and stimulation data
- **Visualization**: Dual y-axis plots (neural signal in blue, stimulation in red)
- **Interpolation**: Previous-value method for stimulation current alignment

#### `signalvsstimulation_part2.m`
**Purpose**: Advanced spectral analysis with gamma band focus

**Enhanced Capabilities**:
- **Segmentation**: 10-second non-overlapping analysis windows
- **Quality Control**: Automatic dropped packet detection per segment
- **Spectral Method**: Welch power spectral density estimation
- **Target Analysis**: Gamma band power (60-90 Hz) quantification

**Output Table Structure**:
| Column | Description |
|--------|-------------|
| `Segment` | Analysis window identifier |
| `FreqL/FreqR` | Frequency vectors for left/right hemispheres |
| `PowerL/PowerR` | Power spectral density values |
| `AvgStimL/AvgStimR` | Average stimulation intensity per segment |

#### `timedomain.m`
**Purpose**: Comprehensive gamma-stimulation correlation analysis

**Analysis Pipeline**:
1. **Data Matching**: Precise TimeDomain-LFP packet alignment
2. **User Selection**: Interactive channel/rate combination selection
3. **Visualization**: 
   - Dual y-axis signal plotting with dropped packet highlighting
   - Zero-line alignment for improved interpretation
4. **Statistical Analysis**:
   - Linear regression with R² reporting
   - Correlation coefficient calculation
   - Optional 2D density heatmaps (≥20 epochs required)

**Output Variables**:
- `gamma_stim_results`: Tabulated analysis results
- `epoch_results`: Detailed epoch-by-epoch data

#### `onespectrumforallpatients.m`
**Purpose**: Revolutionary batch processing with intelligent optimization

**Advanced Architecture**:
- **Memory Optimization**: Selective JSON extraction without full file loading
- **Scalability**: Handles files up to 500MB efficiently
- **Metadata Integration**: Automatic lead location (STN/GPi) and hemisphere detection
- **Data Standardization**: Intelligent column matching across diverse file structures
- **Comprehensive Analysis**: Processes all channel/rate combinations automatically

**Key Innovations**:
- Chunked processing with immediate memory cleanup
- Optimized nearest-neighbor temporal alignment
- Configurable time threshold matching
- Multi-panel visualization with gamma band highlighting

**Output Components**:
- `all_spectrum_results`: Complete spectral dataset
- `channelRateCombos`: Processed combination registry
- Concatenated and standardized data tables

---

### 🔍 Quality Control

#### `droppedpacketdetectionworking.m`
**Purpose**: Communication error identification and data integrity validation

**Detection Methods**:
- Irregular time interval identification
- Missing packet sequence analysis
- Communication error flagging between implant and tablet
- Temporal alignment accuracy verification

#### `threecolumns.m`
**Purpose**: Structured data enhancement with nested organization

**Data Structure Enhancement**:
- Nested table creation within `BrainSenseTimeDomain`
- Three-component sub-tables per row:
  - `DateTime`: Precise temporal reference
  - `StimAmp`: Calculated average stimulation amplitude
  - `GammaPower`: Placeholder for future analysis integration

---

## 📈 Data Output Formats

### Standard Table Structure
All processing scripts generate standardized MATLAB tables optimized for:
- **Temporal Analysis**: Precise timestamp alignment across all data streams
- **Statistical Processing**: Ready integration with MATLAB's statistical toolboxes
- **Visualization**: Direct compatibility with plotting functions
- **Export Compatibility**: Easy conversion to CSV, Excel, or other formats

### Key Data Types
- **Neural Signals**: High-resolution LFP data (typically 250 Hz sampling)
- **Stimulation Parameters**: Rate, amplitude, and timing synchronized to neural data
- **Spectral Data**: Power spectral density across frequency ranges of interest
- **Quality Metrics**: Packet loss statistics and data integrity assessments

---

## 🚀 Usage Examples

### Basic Processing Workflow
```matlab
% 1. Organize raw files
Chan_LoadJSON();

% 2. Process individual file
rawsignalandstimulationintensity2();

% 3. Generate timestamps
realtimestampingv2();

% 4. Analyze signal-stimulation relationships
signalvsstimulation_part1();
```

### Advanced Spectral Analysis
```matlab
% Comprehensive gamma-stimulation correlation
timedomain();

% Batch processing across multiple patients
onespectrumforallpatients();
```

---

## 🤝 Contributing

We welcome contributions from the research community! Please consider:
- **Bug Reports**: Submit issues with detailed descriptions and sample data
- **Feature Requests**: Propose new analysis methods or visualization improvements  
- **Code Contributions**: Follow MATLAB best practices and include comprehensive documentation
- **Validation Studies**: Share results from different patient populations or device configurations

---

## 👨‍🔬 Author

**Andrew Chan**  
*Bioengineering, University of Pittsburgh*  
*Researcher, Hammer Neuromodulation Lab*  
*University of Pennsylvania*

---

**🧠 Advancing Parkinson's Disease Research Through Innovative Data Analysis 🧠**

*For questions, support, or collaboration opportunities, please contact me*
