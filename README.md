# TUS_device_sim_and_char
Code for the Master's Thesis "Characterization and Simulation of a 256-elements Transcranial Ultrasound Stimulation Device"
DON'T FORGET TO REMOVE SENSITIVE DATA
## Introduction / Background and Motivation

In this repository, you will find the code (only MATLAB ??) that was written and used to obtain the results in Maylis Müller's Master's Thesis. 

=> Background

=> Context of the project

## Requirements

- MATLAB (we used R2025a)
- k-Wave (v1.4, cite)
- k-Wave CUDA binaries for GPU-accelerated 3-D simulations

=> other packages ?

## Description

The project's architecture is described below:

```text
 -tu_muellerm_cloud/
  ├── skulls_FUN/
  │   ├── sensor_data/      # Description
  │         ├── el_size_calib/ 
                      ├── pac/ # PAC DATA, pac_ at beginning means file with actual delays; if _PAC at the end means sensor data after                            │         # applying PAC
                      ├── unsteered/ # 
  │                   └── steered/
  │                  
  ├── folder_02/             # Description
  │   ├── subfolder_01/
  │   ├── subfolder_02/
  │   └── subfolder_03/
  ├── folder_03/
  │   ├── subfolder_01/
  │   ├── subfolder_02/
  │   └── subfolder_03/
  
 -tu_muellerm/
  ├── folder_01/
  │   ├── subfolder_01/      # Description
  │   ├── subfolder_02/
  │   └── subfolder_03/
  ├── folder_02/             # Description
  │   ├── subfolder_01/
  │   ├── subfolder_02/
  │   └── subfolder_03/
  ├── folder_03/
  │   ├── subfolder_01/
  │   ├── subfolder_02/
  │   └── subfolder_03/
  ├── folder_04/             # Description
  │   ├── subfolder_01/
  │   ├── subfolder_02/
  │   └── subfolder_03/
  ├── folder_05/             # Description
  │   ├── subfolder_01/
  │   ├── subfolder_02/
  │   └── subfolder_03/
  ├── folder_06/             # Description
  ├── folder_07/             # Description
  ├── folder_08/             # Description
  ├── folder_09/             # Description
  └── utils/             # Description (kind of straigth forward)
```

### Folder contents

Each main folder corresponds to a specific part of the simulation, processing, or analysis workflow. The most relevant scripts and subfolders are described below.

#### `folder_01` — [Folder title]

Contains the code used to [main purpose of the folder].

**Main files**

* `script_01.m` — [Short description of what the script does].
* `script_02.m` — [Short description].
* `script_03.m` — [Short description].

**Subfolders**

* `subfolder_01/` — [Purpose of the subfolder].
* `subfolder_02/` — [Purpose of the subfolder].
* `subfolder_03/` — [Purpose of the subfolder].

#### `folder_02` — [Folder title]

Contains the code used to [main purpose of the folder].

**Main files**

* `script_01.m` — [Short description].
* `script_02.m` — [Short description].

**Subfolders**

* `subfolder_01/` — [Purpose of the subfolder].
* `subfolder_02/` — [Purpose of the subfolder].
* `subfolder_03/` — [Purpose of the subfolder].

#### `folder_03` — [Folder title]

Contains the code used to [main purpose of the folder].

**Main files**

* `script_01.m` — [Short description].
* `script_02.m` — [Short description].

**Subfolders**

* `subfolder_01/` — [Purpose of the subfolder].
* `subfolder_02/` — [Purpose of the subfolder].
* `subfolder_03/` — [Purpose of the subfolder].

#### `folder_04` — [Folder title]

Contains the code used to [main purpose of the folder].

**Main files**

* `script_01.m` — [Short description].
* `script_02.m` — [Short description].

**Subfolders**

* `subfolder_01/` — [Purpose of the subfolder].
* `subfolder_02/` — [Purpose of the subfolder].
* `subfolder_03/` — [Purpose of the subfolder].

#### `folder_05` — [Folder title]

Contains the code used to [main purpose of the folder].

**Main files**

* `script_01.m` — [Short description].
* `script_02.m` — [Short description].

**Subfolders**

* `subfolder_01/` — [Purpose of the subfolder].
* `subfolder_02/` — [Purpose of the subfolder].
* `subfolder_03/` — [Purpose of the subfolder].

#### `folder_06` — [Folder title]

Contains the code used to [main purpose of the folder].

**Main files**

* `script_01.m` — [Short description].
* `script_02.m` — [Short description].

#### `folder_07` — [Folder title]

Contains the code used to [main purpose of the folder].

**Main files**

* `script_01.m` — [Short description].
* `script_02.m` — [Short description].

#### `folder_08` — [Folder title]

Contains the code used to [main purpose of the folder].

**Main files**

* `script_01.m` — [Short description].
* `script_02.m` — [Short description].

#### `folder_09` — [Folder title]

Contains the code used to [main purpose of the folder].

**Main files**

* `script_01.m` — [Short description].
* `script_02.m` — [Short description].

#### `utils`

Contains utility functions shared across multiple parts of the repository.

**Main files**

* `function_01.m` — [Short description].
* `function_02.m` — [Short description].
* `function_03.m` — [Short description].


