import numpy as np
import os
import argparse


def parse_midas_ps_file(ps_filepath):
    """
    Parse a MIDAS parameters file and return a dictionary of the data in key-value pairs
    Args:
        ps_filepath (str): Path to the MIDAS parameters txt file
    Returns:
        ps_params (dict): Dictionary of the parameters. All values are strings and need to\
             be converted to the appropriate data type when used.
    """

    ps_params = {}
    with open(ps_filepath, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            key = parts[0]
            value = " ".join(parts[1:])  # keep as a single string
            ps_params[key] = value

    return ps_params


def read_params_for_radint(ps_filepath):
    """
    Filter the parameters for the radial integration.
    Args:
        ps_filepath (str): Path to the MIDAS parameters txt file
    Returns:
        radint_params (dict): Dictionary of the parameters needed for the radial integration
    """

    ps_params = parse_midas_ps_file(ps_filepath) # parse the MIDAS parameters file

    # Create a dictionary to store the parameters needed for the radial integration
    radint_params = {}
    radint_params["px"] = float(ps_params["px"]) # pixel size in microns
    radint_params["BC"] = [float(i) for i in ps_params["BC"].split()] # First value is Y(horizontal px); second value is Z(vertical px)
    radint_params["Lsd"] = float(ps_params["Lsd"]) # distance from sample to detector in m

    # Read rotation angles (in degrees)
    radint_params["tx"] = float(ps_params["tx"]) # rotation about x-axis in degrees
    radint_params["ty"] = float(ps_params["ty"]) # rotation about y-axis in degrees
    radint_params["tz"] = float(ps_params["tz"]) # rotation about z-axis in degrees
    
    # Read distortion parameters
    radint_params["p0"] = float(ps_params["p0"])
    radint_params["p1"] = float(ps_params["p1"])
    radint_params["p2"] = float(ps_params["p2"])
    radint_params["p3"] = float(ps_params["p3"])

    # For nr. of pixels, first see if NrPixels is present & use that to assume square detector
    if "NrPixels" in ps_params.keys():
        radint_params["NrPixelsY"] = int(ps_params["NrPixels"])
        radint_params["NrPixelsZ"] = int(ps_params["NrPixels"])
    
    # If separate NrPixelsY & NrPixelsZ are present, use those to overwrite the square detector assumption
    if "NrPixelsY" in ps_params.keys():
        radint_params["NrPixelsY"] = int(ps_params["NrPixelsY"])
    if "NrPixelsZ" in ps_params.keys():
        radint_params["NrPixelsZ"] = int(ps_params["NrPixelsZ"])

    # If numPxY & numPxZ are present, use those to overwrite the NrPixelsY & NrPixelsZ
    if "numPxY" in ps_params.keys():
        radint_params["NrPixelsY"] = int(ps_params["numPxY"])
    if "numPxZ" in ps_params.keys():
        radint_params["NrPixelsZ"] = int(ps_params["numPxZ"])

    # Read bad pixel intensity values
    if "BadPxIntensity" in ps_params.keys():
        radint_params["BadPxIntensity"] = float(ps_params["BadPxIntensity"])

    # Read gap pixel intensity values
    if "GapIntensity" in ps_params.keys():
        radint_params["GapIntensity"] = float(ps_params["GapIntensity"])

    return radint_params


def create_rotation_matrices(tx, ty, tz):
    """
    Create the rotation matrices for the given rotation angles
    Args:
        tx (float): rotation about x-axis in degrees
        ty (float): rotation about y-axis in degrees
        tz (float): rotation about z-axis in degrees
    Returns:
        Rx (np.ndarray): rotation matrix about x-axis
        Ry (np.ndarray): rotation matrix about y-axis
        Rz (np.ndarray): rotation matrix about z-axis
    """

    Rx = np.array([
        [1, 0, 0],
        [0, np.cos(np.deg2rad(tx)), -np.sin(np.deg2rad(tx))],
        [0, np.sin(np.deg2rad(tx)), np.cos(np.deg2rad(tx))]
        ])

    Ry = np.array([
        [np.cos(np.deg2rad(ty)), 0, np.sin(np.deg2rad(ty))],
        [0, 1, 0],
        [-np.sin(np.deg2rad(ty)), 0, np.cos(np.deg2rad(ty))]
        ])

    Rz = np.array([
        [np.cos(np.deg2rad(tz)), -np.sin(np.deg2rad(tz)), 0],
        [np.sin(np.deg2rad(tz)), np.cos(np.deg2rad(tz)), 0],
        [0, 0, 1]
        ])

    return Rx, Ry, Rz


def create_distortion_map(RR, EE, p0, p1, p2, p3, px):
    """
    Create the distortion map for the given parameters
    Args:
        RR (np.ndarray): radial map (in microns)
        EE (np.ndarray): polar angle map (spans from 180 to -180 degrees in clockwise direction)
        p0 (float): distortion parameter
        p1 (float): distortion parameter
        p2 (float): distortion parameter
        p3 (float): distortion parameter
        px (float): pixel size in microns
    Returns:
        R_N (np.ndarray): normalized radial map in microns
    """

    RR_N = RR / np.max(RR) # Normalized radial map
    EE_T = 90 - EE # Transformed polar angle map

    # Distortion terms
    dist_gRE = p0 * (RR_N**2) * np.cos(np.deg2rad(2*EE_T)) # 2-fold angular
    dist_hRE = p1 * (RR_N**4) * np.cos(np.deg2rad(4*EE_T + p3)) # 4-fold angular
    dist_kRE = p2 * (RR_N**2) # radial distortion
    dist_fRE = 1 + dist_gRE + dist_hRE + dist_kRE # combined distortion factor

    return dist_fRE
    

def integrate_to_polar_bins(data_arr, RR, EE, dR, dE):
    """
    Integrate the 2D intesity to polar (radial-azimuthal) bins
    Args:
        data_arr (np.ndarray): data array with pixel intensity values
        RR (np.ndarray): radial map (in pixels)
        EE (np.ndarray): polar angle map (in degrees)
        dR (float): radial bin size (in pixels)
        dE (float): azimuthal bin size (in degrees)
    Returns:
        RE_bin_intensity (np.ndarray): 2D array with integrated intensity values. Row index is radial bin, column index is azimuthal bin
        R_bin_centers (np.ndarray): array of radial bin centers
        E_bin_centers (np.ndarray): array of azimuthal bin centers
    """

    R_bin_edges = np.arange(0, np.max(RR) + dR/2, dR)
    R_bin_centers = R_bin_edges[:-1] + dR/2

    E_bin_edges = np.arange(np.min(EE), np.max(EE) + dE/2, dE)
    E_bin_centers = E_bin_edges[:-1] + dE/2

    # The RE_bin_sum is a 2D array where the first dimension is the radial bins and the second dimension is the angular bins
    # So suppose I want lineout for only the angular range of 3rd Eta bin, then I can do: RE_bin_intensity[:,2]
    # If I want lineout for full azimuthal range, then I can do: np.sum(RE_bin_intensity, axis=1)
    
    RE_bin_intensity, _, _ = np.histogram2d(
        RR.ravel(), EE.ravel(), 
        bins=[R_bin_edges, E_bin_edges], 
        weights=data_arr.ravel()
    )

    return RE_bin_intensity, R_bin_centers, E_bin_centers


def do_radial_integration(data, ps_filepath, dR=0.5, dE=5, h5datacontainer="None"):
    """
    Perform radial integration on the data array
    Args:
        data (np.ndarray or str): data array with pixel intensity values or path to the tif/h5/npy file
        ps_filepath (str): path to the MIDAS parameters txt file
        dR (float): radial bin size (in pixels)
        dE (float): azimuthal (Eta) bin size (in degrees)
        h5datacontainer (str): name of the container in the h5 file (needed only if input 'data' is path to an h5 file)

    Returns:
        radial_lineout (np.ndarray): 1D array with integrated intensity values as a function of radial distance;\
            azimuthal bins are compressed into one
        RE_bin_intensity (np.ndarray): 2D array with integrated intensity values. Row index is radial bin, column index is azimuthal bin
        R_bin_centers (np.ndarray): array of radial bin centers
        E_bin_centers (np.ndarray): array of azimuthal bin centers
    """

    if isinstance(data, str):
        if data.endswith(".tif") or data.endswith(".tiff"):
            import imageio
            data_arr = imageio.imread(data)
        if data.endswith(".h5"):
            import h5py
            data_arr = h5py.File(data)[h5datacontainer]
        if data.endswith(".npy"):
            data_arr = np.load(data)

    if isinstance(data, np.ndarray):
        data_arr = data

    radint_params = read_params_for_radint(ps_filepath)

    # YY and ZZ are the pixel coordinates of the detector
    # These are relative to the top left corner pixel as origin
    YYpx, ZZpx = np.meshgrid(np.arange(0, radint_params["NrPixelsY"]), np.arange(0, radint_params["NrPixelsZ"]))

    # yyd and zzd are the physical coordinates of the detector in microns
    # These are relative to the beam center as origin
    yyd = (-YYpx + radint_params["BC"][0])*radint_params["px"]
    zzd = (ZZpx - radint_params["BC"][1])*radint_params["px"]

    # xyz_comb is the combined array of x, y, z physical coordinates of the detector in microns
    # The detector map here is flattened; first row is x, second row is y, third row is z coordinates of all pixels
    xyz_comb = np.array([
        np.zeros_like(yyd.flatten()),
        yyd.flatten(),
        zzd.flatten()
    ])

    # Create the rotation matrices
    Rx, Ry, Rz = create_rotation_matrices(radint_params["tx"], radint_params["ty"], radint_params["tz"])

    # xyz_comb_rot is the combined array of x, y, z physical coordinates of the detector in microns after rotation
    xyz_comb_rot = (Rx @ (Ry @ Rz)) @ xyz_comb

    # X, Y, Z are the cartesian lab coordinates of the detector in microns
    det_shape = (radint_params["NrPixelsZ"], radint_params["NrPixelsY"])
    XX = (radint_params["Lsd"] + xyz_comb_rot[0]).reshape(det_shape)
    YY = xyz_comb_rot[1].reshape(det_shape)
    ZZ = xyz_comb_rot[2].reshape(det_shape)

    RR = (radint_params["Lsd"]/XX) * np.sqrt(YY**2 + ZZ**2)
    EE = np.zeros_like(YY)
    EE[YY <= 0] = np.rad2deg(np.arccos(ZZ[YY<=0] / (np.sqrt(YY[YY<=0]**2 + ZZ[YY<=0]**2))))
    EE[YY > 0] = np.rad2deg(- np.arccos(ZZ[YY>0] / (np.sqrt(YY[YY>0]**2 + ZZ[YY>0]**2))))


    # Create the distortion map
    dist_fRE = create_distortion_map(RR, EE,
                                    radint_params["p0"], radint_params["p1"], radint_params["p2"], radint_params["p3"],
                                    radint_params["px"])

    RR_T = RR * dist_fRE / radint_params["px"]

    RE_bin_intensity, R_bin_centers, E_bin_centers = integrate_to_polar_bins(data_arr, RR_T, EE, dR, dE)

    # Just the radial lineout; compressing all azimuthal bins into one
    radial_lineout = np.sum(RE_bin_intensity, axis=1)

    return radial_lineout, RE_bin_intensity, R_bin_centers, E_bin_centers


def main():
    """
    Main function to perform radial integration
    """

    parser = argparse.ArgumentParser(
        description='Radial integration',
    )

    parser.add_argument('-data', type=str, required=True, help='Path to the tif/h5/npy file')
    parser.add_argument('-ps_filepath', type=str, required=True, help='Path to the MIDAS parameters txt file')
    parser.add_argument('-output_dir', type=str, required=False, default="out", help='Output directory')

    parser.add_argument('-dR', type=float, required=False, default=0.5, help='Radial bin size (in pixels)')
    parser.add_argument('-dE', type=float, required=False, default=5, help='Azimuthal (Eta) bin size (in degrees)')
    parser.add_argument('-h5datacontainer', type=str, required=False, default="None", help='Name of the container in the h5 file (needed only if input data is path to an h5 file)')
    parser.add_argument('-plot_results', type=str, required=False, default="yes", help='Plot the results')
    args = parser.parse_args()

    radial_lineout, RE_bin_intensity, R_bin_centers, E_bin_centers = do_radial_integration(args.data, args.ps_filepath, args.dR, args.dE, args.h5datacontainer)

    print(f"Radial lineout shape: {radial_lineout.shape}")
    print(f"RE_bin_intensity shape: {RE_bin_intensity.shape}")
    print(f"R_bin_centers shape: {R_bin_centers.shape}")
    print(f"E_bin_centers shape: {E_bin_centers.shape}")

    # Save the results
    os.makedirs(args.output_dir, exist_ok=True)
    np.save(os.path.join(args.output_dir, "radial_lineout.npy"), radial_lineout)
    np.save(os.path.join(args.output_dir, "RE_bin_intensity.npy"), RE_bin_intensity)
    np.save(os.path.join(args.output_dir, "R_bin_centers.npy"), R_bin_centers)
    np.save(os.path.join(args.output_dir, "E_bin_centers.npy"), E_bin_centers)

    print(f"Results saved to '{args.output_dir}' directory")

    if args.plot_results.lower() in ["yes", "true", "y", "t", "1"]:
        from matplotlib import pyplot as plt
        plt.figure(figsize=(10, 5))
        plt.plot(R_bin_centers, radial_lineout)
        plt.xlabel("Radial distance (pixels)")
        plt.ylabel("Integrated intensity")
        plt.title("Radial lineout")
        plt.show()


if __name__ == "__main__":
    main()