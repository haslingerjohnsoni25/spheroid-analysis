# Spheroid Fragmentation & Entropy Analysis Script

## Overview
This Python script analyzes spheroid images from cancer research experiments to quantify structural metrics such as fragmentation, entropy (complexity), number of internal holes, and total area. It is optimized for evaluating treatment effects (e.g., immune-based cytotoxicity) using time-series data of control and experimental groups.

The script uses OpenCV and skimage for image processing, with SciPy, seaborn, and pandas for data analysis and visualization. Results are saved in Excel format along with annotated image overlays and line plots of entropy trends.

## Features
- Histogram equalization, Gaussian blurring, and adaptive thresholding
- Morphological cleanup and watershed segmentation
- Entropy measurement using local texture (rank-based filter)
- Hole detection via contour hierarchy
- Area computation based on calibrated µm² conversion
- Batch processing of nested well folders for two experimental groups
- Excel export and entropy vs. time plots for group comparison

## Requirements
Install dependencies using pip:
```bash
pip install opencv-python numpy pandas matplotlib seaborn scikit-image scikit-learn statsmodels termcolor
```

## Input
- Two top-level folders: one for the Control group and one for the Experimental group.
- Each folder must contain subfolders representing different wells.
- Each subfolder should contain spheroid image files (`.jpg`, `.jpeg`, `.png`).

## Output
- `fragmentation_results.xlsx`: A table with columns:
  - `Image`: Well identifier
  - `Session`: Timestamp (from filename)
  - `Group`: Control or Experimental
  - `Area_um2`: Total spheroid mask area in µm²
  - `Fragments`: Number of segmented regions
  - `Entropy`: Mean local texture entropy
  - `Holes`: Internal holes within spheroids
- Annotated images saved as `<original_name>_frag_detected.png` with visual overlays
- `entropy_fragments_plot.png`: Line plot showing entropy trends by group

## How to Use
1. Modify the following directory paths in the script to match your setup:
```python
control_dir = "/path/to/prevention_control"
experimental_dir = "/path/to/prevention_experimental"
```
2. Run the script from a terminal:
```bash
python spheroid_analysis_nested_v1.1.py
```

## Calibration
The following conversion is used to compute area in micrometers:
```python
um_per_px = 200.0 / 298.412  # µm per pixel
area_conversion = um_per_px ** 2  # µm² per pixel²
```
You can adjust this ratio based on your imaging system.

## Example Application
This script was developed to support image-based assessment of immune-mediated cytotoxicity in a 3D MyC-CaP prostate tumor spheroid model. It helps quantify the impact of a peptide-based vaccine by analyzing spheroid degradation and complexity over time.