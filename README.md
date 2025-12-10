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

### Requirements

```bash
pip install numpy matplotlib imageio h5py
```

# Run radial integration using below command.

```bash
python radial_integration.py -data test_data.npy -ps_filepath test_ps.txt
```

If you need radial integration as part of your bigger workflow, consider using the 'do_radial_integration' function in the 'radial_integration.py' file.

Required Arguments:

-data: Path to the input image file (TIFF, HDF5, or NPY format)
-ps_filepath: Path to the MIDAS parameters text file


Optional Arguments:

-output_dir: Output directory (default: "out")
-dR: Radial bin size in pixels (default: 0.5)
-dE: Azimuthal bin size in degrees (default: 5)
-h5datacontainer: HDF5 dataset name (required for HDF5 files)
-plot_results: Show plots ("yes"/"no", default: "yes")


# Detector geometry
px 50.0                    # Pixel size in microns
BC 1024.5 1024.5         # Beam center coordinates (Y, Z)
Lsd 0.5                   # Sample-to-detector distance in meters

# Detector rotations (degrees)
tx 0.0                    # Rotation about X-axis
ty 0.0                    # Rotation about Y-axis  
tz 0.0                    # Rotation about Z-axis

# Distortion parameters
p0 0.0                    # 2-fold angular distortion
p1 0.0                    # 4-fold angular distortion
p2 0.0                    # Radial distortion
p3 0.0                    # 4-fold angular phase

# Detector size
NrPixels 2048            # Number of pixels (square detector)
# OR
NrPixelsY 2048           # Horizontal pixels
NrPixelsZ 2048           # Vertical pixels
