# Image Analysis Pipeline for Synthetic Patterning Experiments

This repository contains Fiji/ImageJ macros and Python scripts used for the analysis of microscopy images associated with the preprint:

**Swedlund B. et al.**
*A novel reaction-diffusion architecture for engineering self-organized patterns in mammalian cells*
https://www.biorxiv.org/content/10.64898/2026.05.24.727552v1

The workflow is designed for multi-channel fluorescence microscopy images acquired on a Zeiss microscope and exported as `.czi` files. The analysis pipeline consists of:

1. Conversion and TIF/JPEG export of raw microscopy images
2. Generation of binary masks from fluorescence channels
3. Quantification of area fraction and analyze particles
4. Spatial pattern analysis using autocorrelation, cross-correlation, and overlap metrics

---

## Fiji Macros

### 1. `czi_file_merge_adjust_channels.ijm`

Converts multi-channel `.czi` files into merged TIFF and JPEG images using the Bio-Formats importer. The script performs background subtraction, intensity normalization, channel merging, and image export.

Key parameters:

* **channelOrder** and **channelColors**: define channel mapping and display colors.
* **backgroundRolling**: radius of the rolling-ball background subtraction.
* **lower/upper percentiles**: determine display scaling for fluorescent channels.
* **threshPos**: intensity threshold used to determine whether a pixel is considered activated
* **pctSwitch**: fraction of positive pixels used to switch between display scaling modes (thresholds 1 or 2).

---

### 2. `tif_making_binary_mask.ijm`

Generates binary masks from fluorescence images through thresholding, Gaussian smoothing, and morphological operations.

Key parameters:

* **batch_mode**: process a single image or all images within a folder.
* **channelsToProcess**: channels for which masks should be generated.
* **min_t_arr / max_t_arr**: critical threshold values defining what is considered a positive pixel.
* **cropping options**: crop images to circular or rectangular regions of interest.
* **dilate_erode_iterations**: helps to remove small objects and to smoothe mask boundaries.
* **first_gaussian_blur** and **second_gaussian_blur**: determine the degree of smoothing. Important parameters, to be changed depending on the resolution of the images and the wanted granularity.

---

### 3.1 `mask_area_fraction.ijm`

Measures the fraction of image area occupied by a binary mask and exports the results to a CSV file.

---

### 3.2 `mask_analyze_particles.ijm`

Runs Fiji's **Analyze Particles** workflow on binary masks and exports particle measurements, including:

* Area
* Centroid position
* Shape descriptors

Optional cropping can be applied before analysis.

---

## Python Analysis Scripts

### 4.1 `run_autocorrelation.py`

Computes the spatial autocorrelation of a single image channel and generates radial autocorrelation profiles.


---

### 4.2 `run_cross_correlation.py`

Computes the spatial cross-correlation between two image channels (e.g., GFP and mCherry).


---

### 4.3 `run_overlap_analysis.py`

Computes overlap and colocalization metrics between two binary masks.

Reported metrics include:

* Fractional area occupied by each mask
* Fractional overlap area
* Manders M1 and M2 coefficients
* Dice–Sørensen coefficient
* Cohen's kappa coefficient
* Expected overlap under random placement
* Excess overlap and normalized colocalization coefficients


---

## Helper Modules

### `analysis.py`

Core analysis functions used throughout the repository, including:

* Image preprocessing
* Autocorrelation and cross-correlation calculations
* Radial profile extraction
* Colocalization and overlap metrics

---

### `viz.py`

Visualization utilities used to generate figures for:

* Autocorrelation analysis
* Cross-correlation analysis
* Radial profile summaries

These functions produce publication-ready plots directly from the outputs of `analysis.py`.

---

## Dependencies

### Fiji

Required plugins:

* Bio-Formats Importer

### Python

Major dependencies:

* numpy
* scipy
* pandas
* matplotlib
* h5py
* scikit-image
