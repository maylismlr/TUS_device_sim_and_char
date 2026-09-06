# TUS_device_sim_and_char

In this repository, you will find the code that was developed and used to obtain the results in the Master's Thesis "Characterization and Simulation of a 256-elements Transcranial Ultrasound Stimulation Device". This thesis was written by Maylis Müller as part of the Master's degree in Life Sciences Engineering at EPFL.

## Introduction / Background and Motivation

 work mentioned above.

=> Background

=> Context of the project

=> BASICALLY I WILL ADD MY ABSTRACT

## Requirements

- MATLAB (we used R2025a)
- k-Wave (v1.4, cite)
- k-Wave CUDA binaries for GPU-accelerated three-dimensional simulations cite

## Description

The project's architecture is described below:

```text
 -tu_muellerm/
  ├── convergence_testing/             
  │   ├── cvg_free-field/
  │   │   ├── CFL_number/
  │   │   ├── dx/
  │   │   └── PML_size/
  │   └── cvg_with_skull/
  │       ├── CFL_number/
  │       └── PML_size/
  │
  ├── delay_calibration/
  │   ├── delay_calib_results/      
  │   └── pressure_fields_corrected/  
  │
  ├── element_size_calibration/             
  │   ├── broad_sweep/
  │   ├── narrow_sweep/
  │   ├── narrower_sweep/
  │   ├── element_size_calibration_figures_39/ 
  │   └── element_size_calibration_figures_nom_size/
  │
  ├── intensity_calibration/          
  │   └── input_calib_pressure_fields/
  │
  ├── little_TUS_anim_for_fun/           
  │
  ├── PAC/                         # Phase aberration correction
  │   └── PAC_plots/
  │
  ├── skull_models/             
  │   ├── python_files/            # To map acoustic properties to HU values for the skulls, and to plot simulation pressure fields
  │   └── skull_model_results/
  │
  ├── skull_simulation_videos/     # Contains short 2D videos showing ultrasound wave attenuation through muCT skull        
  │
  ├── to_reproduce_paper_figures/  # To reproduce Figures 6 and 7 from Tretbar et al.
  │   └── simulation_plots/
  │       ├── Figure_6/
  │       └── Figure_7/
  │
  └── utils/                       # Helper functions
```

### Folder contents

Each main folder corresponds to a specific part of the simulation, processing, or analysis workflow. The subfolders and their most relevant scripts are described below.

#### `convergence_testing`

Contains the code used to assess numerical convergence of the acoustic simulations. For free-field simulations, convergence was evaluated with respect to the spatial grid spacing, CFL number, and perfectly matched layer (PML) size. For simulations including a skull model, convergence was evaluated with respect to the CFL number and PML size.

The scripts used to analyse and plot the convergence results are located at the top level of the corresponding folders, while the resulting figures are stored in the associated subfolders.

#### `delay_calibration` 

Contains the code and results used to calibrate the focal-depth input of the simulated TRUST array against the experimentally recorded focal positions. The resulting calibration is stored in \texttt{delay\_calibration.mat} and is used to convert a requested device focal depth to the corresponding corrected simulation input.

**Main files**

* `delay_calibration.m` — Generates \texttt{delay\_calibration.mat}, which contains the mapping used to determine the corrected axial focal input.
* `get_rec_max_p_loc_Fraunhofer.m` — Extracts the spatial peak-pressure locations from the experimental TRUST recordings.
* `get_validation_delays` - Generates the corrected simulation pressure fields used to validate the focal-depth calibration.

**NB**: The hydrophone axial pressure recordings used are not distributed with the repository.

#### `element_size_calibration` 

Contains the simulations and analysis used to determine the effective radiating element size of the TRUST array by comparison with experimentally measured axial pressure profiles.

**Main files**

* `fit_element_size_profile.m` — actual code to run to get element size that minimizes L_2 error from chosen interval.

**NB**: The hydrophone axial pressure recordings used (same as above) are not distributed with the repository.

#### `intensity_calibration` 

Contains the code used to evaluate the pressure-amplitude calibration implemented in \texttt{trust\_p0.m}. The analysis verifies the relationship between the requested device voltage, focal depth, and source pressure used in the simulations based on the available TÜV SÜD hydrophone measurements.

The folder also contains element-size-calibrated free-field simulation pressure fields obtained at several focal depths.


#### `little_TUS_anim_for_fun` 

Contains a simplified two-dimensional k-Wave simulation used to generate an illustrative animation of pulsed transcranial ultrasound stimulation. The corresponding animation is also included in the folder.

This material is provided for visualisation purposes and is not used in the quantitative analyses presented in the thesis.

#### `PAC`

Contains the implementation and analysis of phase-aberration correction (PAC). The code estimates element-wise phase corrections from skull-dependent acoustic propagation and uses the resulting delays to perform corrected TRUST simulations.

The folder also contains scripts for extracting PAC-related metrics, comparing corrected and uncorrected simulations, and generating the corresponding figures. Generated figures are stored in \texttt{PAC\_plots/}.

**Main files**

* `tFUS_3D_[skull_model_name]_PAC.m` — Runs the simulations with PAC and plots figures.


#### `skull_models` 

Contains the code used to generate the acoustic skull models and perform transcranial ultrasound simulations through these models.

**Main subfolder**

* `python_files/`- Contains the Python scripts used to map image-derived skull values to acoustic properties and to visualise the resulting acoustic-property maps and simulation pressure fields.

**Main files**

* `tFUS_3D_[skull_model_name]`- Runs a three-dimensional TRUST simulation using the corresponding skull model as the heterogeneous acoustic medium.

#### `skull_simulation_videos`

Contains short two-dimensional animations illustrating ultrasound propagation and attenuation through the $\mu$CT-derived skull model.

#### `to_reproduce_paper_figures`

Contains simulations and plotting scripts used to reproduce selected results from Tretbar et al., in particular Figures~6 and~7 describing the acoustic field of the TRUST device.

The generated simulation figures are stored in:
* ```simulation_plots/Figure_6/```
* ```simulation_plots/Figure_7/```

#### `utils`

Contains utility functions shared across multiple parts of the repository, including TRUST array and source construction, signal generation, pressure calibration, focal-field analysis, phase-aberration correction, and data visualisation.

**Main files**

* `convert_fig_to_png.m` — Converts saved MATLAB `.fig` figures to high-resolution `.png` files.
* `build_karray.m` — Builds the 256-element TRUST transducer geometry as a `kWaveArray` using the specified element dimensions and pitch.
* `create_signal.m` — Generates the TRUST temporal excitation signal, including tone bursts, signal envelopes, pulse timing, and pulse trains.
* `create_source.m` — Creates the TRUST k-Wave source, applying geometric focusing delays and, optionally, element-wise phase-aberration correction delays, and distributing the excitation signal over the array elements using `build_karray.m` and `create_signal.m`.
* `focal_vol_dB.m` — Computes the connected focal volume at a specified dB threshold using both voxel counting and a marching-cubes isosurface estimate.
* `build_afp_source_bli.m` — Builds band-limited source amplitude and phase matrices for `acousticFieldPropagator`, avoiding grid-dependent rasterisation of the array elements. Used for element-size calibration.
* `get_fwhm.m` — Computes the full width at half maximum (FWHM) of a pressure profile using interpolated half-maximum crossings.
* `trust_p0.m` — Computes the calibrated acoustic source-pressure amplitude corresponding to a requested TRUST excitation voltage and commanded focal depth.
* `trust_reference_data.m` — Stores the hydrophone measurements and simulated reference data used for the TRUST pressure-amplitude calibration.
* `get_DSC.m` — Computes the Dice similarity coefficient between two binary focal-region masks.
* `get_focal_mask.m` — Extracts the connected thresholded focal region containing a specified pressure peak.
* `estimate_pac_phase_reversal.m` — Estimates element-wise skull-induced phase shifts and PAC delays using reciprocal water and skull simulations.
* `focal_metrics.m` — Extracts the peak pressure and focal position from a 3D pressure field, including sub-grid focal-position refinement. Used primarily for the convergence analysis.
