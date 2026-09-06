# TUS_device_sim_and_char

In this repository, you will find the code that was developed and used to obtain the results in the Master's Thesis "Characterization and Simulation of a 256-elements Transcranial Ultrasound Stimulation Device". This thesis was written by Maylis Müller, MSc student at EPFL, to obtain the Life Sciences Engineering diploma.

## Introduction / Background and Motivation

 work mentioned above.

=> Background

=> Context of the project

=> BASICALLY I WILL ADD MY ABSTRACT

## Requirements

- MATLAB (we used R2025a)
- k-Wave (v1.4, cite)
- k-Wave CUDA binaries for GPU-accelerated 3-D simulations

=> other packages ? NO

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

Contains the code used to verify convergence for simulation parameters, with and without skull in the medium. For water-only, convergence was assessed for dx, CFL number and PML size. Only PML size and CFL number were evaluated with the skull. The script to plot is at the top of the folder, and the resulting figures in each subfolder.

#### `delay_calibration` 

Contains the code to obtain .mat delay calibration file, as well as the code to get the .mat simulated pressure fields after calibration and the verification figures.

**Main files**

* `delay_calibration.m` — to get delay_calibration.mat, needed to infer corrected focal input.
* `get_rec_max_p_loc_Fraunhofer.m` — to get the device's recorded spatial peak.
* `get_validation_delays` - to get the corrected simulation pressure profiles.

#### `element_size_calibration` 

**Main files**

* `fit_element_size_profile.m` — actual code to run to get element size that minimizes L_2 error from chosen interval.

**NB**: The hydrophone axial pressure recordings used for the element-size calibration are not divulged here.


#### `intensity_calibration` 

Contains the code used to verify that the function `trust_p0.m` located in utils changes the input pressure to match TÜV SÜD recordings, and check that our rationale works. Also contains element-size calibrated simulation pressure fields at different foci.


#### `little_TUS_anim_for_fun` 

Contains the code used to make video of "TUS" pulsing in 2D with k-Wave, and the video itself.

#### `PAC`

Contains the code used to get corrected delays depending on the skull in the medium, to then run a simulation with these corrected inputs, as well as to plot the wanted figures and extract wanted metrics. Also contains a subfolder with all the figures.


#### `skull_models` 

Contains the code used to map the HU values from each skull model to their acoustic properties (and plot), and to then run the simulations with the models in the medium and save pressure outputs (and plot). 

**Main subfolder**

* `python_files/`- to obtain acoustic properties

**Main files**

* `tFUS_3D_[skull_model_name]`- to run the simulation with corresponding skull model

#### `utils`

Contains utility functions shared across multiple parts of the repository, including TRUST source construction, calibration, focal-field analysis, phase-aberration correction, and data visualisation.

**Main files**

* `convert_fig_to_png.m` — Converts saved MATLAB `.fig` figures to high-resolution `.png` files.
* `build_karray.m` — Builds the 256-element TRUST transducer geometry as a `kWaveArray` using the specified element dimensions and pitch.
* `create_signal.m` — Generates the TRUST temporal excitation signal, including tone bursts, envelopes, pulse timing, and pulse trains.
* `create_source.m` — Creates the TRUST k-Wave source, applying geometric or pac-corrected focusing delays and distributing the excitation signal over the array elements (using `build_karray.m` and `create_signal.m`).
* `focal_vol_dB.m` — Computes the connected focal volume at a configurable dB threshold, including voxel and marching-cubes volume estimates.
* `build_afp_source_bli.m` — Builds band-limited source amplitude and phase matrices for `acousticFieldPropagator`, avoiding grid-dependent rasterisation of the array elements. Used for element-size calibration.
* `create_source_pac.m` — Creates a TRUST source combining geometric focusing delays with element-wise phase-aberration correction delays.
* `get_fwhm.m` — Computes the full width at half maximum (FWHM) of a pressure profile using interpolated half-maximum crossings.
* `trust_p0.m` — Computes the calibrated source pressure amplitude from the requested TRUST voltage and focal depth.
* `trust_reference_data.m` — Stores the hydrophone measurements and simulated reference data used for TRUST amplitude calibration.
* `get_DSC.m` — Computes the Dice similarity coefficient between two binary focal-region masks.
* `get_focal_mask.m` — Extracts the connected thresholded focal region containing a specified pressure peak.
* `estimate_pac_phase_reversal.m` — Estimates element-wise skull-induced phase shifts and PAC delays using reciprocal water and skull simulations.
* `focal_metrics.m` — Extracts peak pressure and focal position from a 3D pressure field, including sub-grid focal-position refinement, used in convergence study.
