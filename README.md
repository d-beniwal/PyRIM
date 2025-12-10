# PyRIM
#Python implementation of Radial Integration routine from MIDAS

A Python tool for performing radial integration of 2D X-ray diffraction patterns into 1D or 2D polar coordinate systems. This tool is designed to work with MIDAS parameter files and supports various image formats.

## Features

- **Radial Integration**: Convert 2D diffraction patterns to 1D radial profiles or 2D polar coordinates
- **Detector Geometry Correction**: Handles detector rotations, distortions, and calibration parameters
- **Multiple Input Formats**: Supports TIFF, HDF5, and NumPy array files
- **Flexible Binning**: Customizable radial and azimuthal bin sizes
- **MIDAS Integration**: Compatible with MIDAS parameter files for detector calibration
- **Visualization**: Optional plotting of results

## What's missing?

The current implementation does not handle overlap of multiple bins on a pixel. In the MIDAS C/C++ code, this is handled through polygon intersection and shoelace area calculations.

## Requirements

```bash
pip install numpy matplotlib imageio h5py
```

# Run radial integration using below command.

```bash
python radial_integration.py -data test_data.npy -ps_filepath test_ps.txt
```

If you need radial integration as part of your bigger workflow, consider using the 'do_radial_integration' function in the 'radial_integration.py' file.

## Required Arguments

| Argument | Description |
|----------|-------------|
| `-data` | Path to the input image file (TIFF, HDF5, or NPY format) |
| `-ps_filepath` | Path to the MIDAS parameters text file |

## Optional Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `-output_dir` | string | `"out"` | Output directory for results |
| `-dR` | float | `0.5` | Radial bin size in pixels |
| `-dE` | float | `5` | Azimuthal bin size in degrees |
| `-h5datacontainer` | string | `"None"` | HDF5 dataset name (required for HDF5 files): path to data variable within hdf5 |
| `-plot_results` | string | `"yes"` | Show plots ("yes"/"no") |


# The MIDAS parameter file used as an input must have the following parameters defined in it

## Detector Geometry

| Parameter | Value | Description |
|-----------|-------|-------------|
| `px` | `200.0` | Pixel size in microns |
| `BC` | `1024.5 1024.5` | Beam center pixel coordinates (Y, Z) |
| `Lsd` | `820000` | Sample-to-detector distance in microns |
| `tx` | `0.124` | Rotation about X-axis (degrees) |
| `ty` | `-0.321` | Rotation about Y-axis (degrees) |
| `tz` | `0.021` | Rotation about Z-axis (degrees)|
| `p0` | `6.357880806e-06` | Distortion parameter for 2-fold angular distortion term |
| `p1` | `6.210839248e-05` | Distortion parameter for 4-fold angular distortion term |
| `p2` | `-0.00074891983` | Distortion parameter for radial distortion term |
| `p3` | `150.0` | Distortion parameter for 4-fold angular distortion term |
| `NrPixels` | `2048` | Number of pixels (square detector) |
| **OR** | | |
| `NrPixelsY` | `1650` | Horizontal pixels |
| `NrPixelsZ` | `1460` | Vertical pixels |