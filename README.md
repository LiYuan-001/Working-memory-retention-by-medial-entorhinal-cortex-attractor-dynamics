# Working-memory-retention-by-medial-entorhinal-cortex-attractor-dynamics

Code accompanying the paper on working memory retention by medial entorhinal cortex attractor dynamics.

## Overview

This repository contains MATLAB code used for the analyses and figure generation in the paper.

Because of GitHub storage limitations, the raw and processed data files are not stored directly in this repository. All data files are provided through an open-access OneDrive folder.

## Data access

Due to space limitations, all data files are available in the following open-access OneDrive folder:

[Open-access OneDrive data folder](https://ucsdcloud-my.sharepoint.com/:f:/r/personal/lyuan_ucsd_edu/Documents/Working-memory-retention-by-medial-entorhinal-cortex-attractor-dynamics/Data?csf=1&web=1&e=rl2vmB)

Please download the required data files before running the code.

## Repository structure

- `Dependencies/`  
  Contains required external functions and supporting code.

- `Figure2_UMAP.m`  
  MATLAB script for analyses and plots related to Figure 2.

- `Figure6_seqNMF_existingRun.m`  
  MATLAB script for analyses and plots related to Figure 6.

- `Figure7_seqEvent_corr.m`  
  MATLAB script for analyses and plots related to Figure 7.

## Requirements

- MATLAB version 2020b or above
- Required toolboxes and third-party functions included in `Dependencies/` or otherwise noted in the code

## How to run

1. Clone or download this repository.
2. Download the data files from the open-access OneDrive folder above.
3. Add both the `Data` folder and the `Dependencies` folder to your local MATLAB path.
4. Run each MATLAB code file directly to reproduce the corresponding analyses and figures.

In MATLAB, you can add the folders to the path, for example:

```matlab
addpath(genpath('path_to/Dependencies'));
addpath(genpath('path_to/Data'));
