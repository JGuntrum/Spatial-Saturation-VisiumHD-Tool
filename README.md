# SpatialSatCheck



# Visium HD (2 µm) binomial subsampling

This repository contains a Jupyter notebook that generates binomially subsampled versions of **segmented Visium HD (2 µm)** spatial transcriptomics data using **10x Genomics Space Ranger** outputs. The subsampled datasets can be used to derive sequencing-depth–dependent metrics and to visualize how data quality changes across effective depth.

## Overview

The notebook implements a downsampling strategy that operates at the level of molecule/UMI read evidence to emulate lower sequencing depth. For a set of target depths (commonly expressed as reads per cell) and optional replicate draws, it produces:

- subsampled count matrices (saved as `.h5ad`)
- per-depth summary metrics and plots
- small manifest/index files describing generated outputs

## Inputs

- A **Space Ranger `outs/` directory** from the relevant Visium HD run (including `molecule_info.h5` and associated mapping files).
- A **reference AnnData (`.h5ad`)** representing the segmented 2 µm cells/bins to which molecules are assigned (used for consistent cell/gene indexing and metadata).

## Outputs

Typical outputs include:

- Multiple subsampled `.h5ad` files (one per target depth × replicate)
- Tables summarizing metrics per depth (and aggregated across replicates)
- Plots showing metric trends across depth
- A small schema/manifest describing the produced artifacts

## Installation (conda)

A `requirements.txt` file is provided. Create and activate a conda environment, then install the requirements inside it:

```bash
conda create -n visiumhd-subsample --file requirements.txt
conda activate visiumhd-subsample
