#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This file contains routines for manipulating SWOT SSH (Sea Surface Height) data.
Based on Jinbo Wang's SWOT-OpenToolkit code.

Functions:
1. find_swaths: Finds SWOT passes intersecting a geographical bounding box using orbit shapefiles.
2. download_passes: Downloads SWOT passes from the AVISO+ server using SFTP.
3. clean_incomplete_files: Deletes incomplete files from a local directory.

Author: Jinbo Wang (First version), Tatsu
Date: First version: 11.15.2024

Dependencies:
    - numpy
    - xarray
    - paramiko
    - geopandas
    - shapely
"""

# Import Python libraries
import geopandas as gpd  # For working with geospatial data
import shapely.geometry as geometry  # For defining bounding boxes and geometry
import paramiko  # For SFTP connections
import os  # For file system operations
import xarray as xr  # For working with multidimensional datasets like SWOT data
import traceback 

# Import local modules
import sys
sys.path.append('./')
import swot_utils
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function: find_swaths
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def find_swaths(sw_corner, ne_corner, path_to_sph_file="./orbit_data/sph_science_swath.zip"):
    """
    Finds SWOT passes that intersect a given geographical bounding box.

    Parameters
    ----------
    sw_corner : list of float
        [longitude, latitude] of the southwest corner of the bounding box.
    ne_corner : list of float
        [longitude, latitude] of the northeast corner of the bounding box.
    path_to_sph_file : str, optional
        Path to the SWOT shapefile that contains orbit information (default: "./orbit_data/sph_science_swath.zip").

    Returns
    -------
    pass_IDs_list : list of str
        List of pass IDs that intersect the bounding box.

    Notes
    -----
    - Uses GeoPandas to load the shapefile and filter passes based on intersection with the bounding box.
    - The bounding box is extended by 0.2° to account for potential nadir data overlaps.
    """
    try:
        # Load the SWOT shapefile containing orbital data
        gdf_karin = gpd.read_file(path_to_sph_file)
    except:
        print("Error: Ensure the shapefile exists in the correct directory.")
        return []

    # Define the bounding box as a Shapely geometry
    bbox = geometry.box(sw_corner[0], sw_corner[1], ne_corner[0], ne_corner[1])

    # Extend the bounding box by 0.2° for better coverage
    extended_bbox = geometry.box(
        bbox.bounds[0] - 0.2, bbox.bounds[1] - 0.2,
        bbox.bounds[2] + 0.2, bbox.bounds[3] + 0.2
    )

    # Filter the GeoDataFrame for passes that intersect the bounding box
    overlapping_segments = gdf_karin[gdf_karin.intersects(extended_bbox)]

    # Extract pass IDs and format them as 3-digit strings with leading zeros
    pass_IDs_list = []
    for foo in overlapping_segments["ID_PASS"].values:
        pass_IDs_list.append(str(foo).zfill(3))

    return pass_IDs_list
    

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function: download_passes
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def download_passes(pass_ID, cycle="001", remote_path="swot_products/l3_karin_nadir/l3_lr_ssh/v2_0_1/Unsmoothed",
                    save_path=f"~/scratch/SWOT_L3_v2_0_1/Unsmoothed/cycle_001", hostname="ftp-access.aviso.altimetry.fr",
                    port=2221, username="tdmonkman@uchicago.edu", password="2prSvl",
                    variables=None,
                    subset=False, lat_lims=False, lon_lims=False, trim_suffix="",**kwargs):
    """
    Downloads SWOT passes from the AVISO+ FTP server using SFTP.

    Parameters
    ----------
    pass_ID : str
        Numeric ID of the pass to download (e.g., "001", "002").
    cycle : str, optional
        Numeric ID of the cycle to download (default: "001").
    remote_path : str, optional
        Path to SWOT data on the AVISO+ server (default: location of "Unsmoothed" products).
    save_path : str, optional
        Local directory to save the downloaded files (default: "../../SWOT_L3/Unsmoothed/cycle_001").
    hostname : str, optional
        FTP server hostname (default: "ftp-access.aviso.altimetry.fr").
    port : int, optional
        Port for SFTP connection (default: 2221).
    username : str, optional
        Username for FTP server authentication.
    password : str, optional
        Password for FTP server authentication.
    subset : bool, optional
        Whether to subset swaths by latitude (default: False).
    lat_lims : list of float, optional
        [min_lat, max_lat] for subsetting if subset=True (default: False).
    trim_suffix : str, optional
        Suffix to add to trimmed files (default: "trimmed").

    Output
    ------
    Saves files to `save_path` and logs skipped passes in 'skipped_swaths.txt'.
    """
    # Extract SWOT release version from remote_path (hacky, should improve)
    try:
        version = int(remote_path.split("_lr_ssh")[0][-1])
    except:
        print(f"Invalid SWOT release version in path: {remote_path}")
        return

    # Define the remote file naming pattern based on product type
    if "Unsmoothed" in remote_path:
        target_remote_file = f"SWOT_L{version}_LR_SSH_Unsmoothed_{cycle}_{pass_ID}"
    elif "Expert" in remote_path:
        target_remote_file = f"SWOT_L{version}_LR_SSH_Expert_{cycle}_{pass_ID}"

    # Initialize SSH connection
    print(f"Attempting SSH connection for target remote file {target_remote_file}...")
    with paramiko.SSHClient() as ssh_client:
        ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh_client.connect(hostname, port, username, password, timeout=30)

        # Open SFTP client
        sftp_client = ssh_client.open_sftp()

        # Define the remote directory for the cycle
        remote_path_cycle = f"{remote_path}/cycle_{cycle}"
        try:
            sftp_client.stat(remote_path_cycle)  # Check if the directory exists
            print(f"Found cycle_{cycle} on remote server")
        except:
            print(f"Cannot find path for cycle_{cycle} on remote server.")
            return

        # List all files in the cycle directory on the server
        remote_files = sftp_client.listdir(remote_path_cycle)
        print(f"Looking for matches to {target_remote_file}...")

        # Loop through files and download matches
        for remote_file in remote_files:
            if target_remote_file not in remote_file:
                continue  # Skip non-matching files

            print(f"Found remote file: {remote_file}")
            local_file_path = f"{save_path}/{remote_file}"

            # Skip if file already exists locally
            if os.path.isfile(local_file_path):
                print(f"{local_file_path} already exists! Skipping for now...")
                continue

            try:
                # Download full or subset file
                if not subset:
                    sftp_client.get(f"{remote_path_cycle}/{remote_file}", local_file_path)
                    print(f"Downloaded full {remote_file}")
                else:
                    # Subset file by latitude
                    temp_file = f"./tmp_{trim_suffix}.nc"
                    sftp_client.get(f"{remote_path_cycle}/{remote_file}", temp_file)
                    # Open and trim datasets
                    with xr.open_dataset(temp_file) as swath:
                        if variables:
                            swath = swath[variables]
                        if "cross_track_distance" in swath.keys():
                            swath = swath.drop_vars("cross_track_distance")
                        if isinstance(lon_lims, list):
                            trimmed_swath = swot_utils.xr_subset(swath, lat_lims, lon_lims).load()          
                        else:
                            trimmed_swath = swot_utils.subset(swath, lat_bounds=lat_lims).load()
                            # The 'xr_subset' function returns a 'None' type if it can't find
                            # data in the lat bounds. NOE: I'm using a 'try' block here so 
                            # if a Nonetype is returned there may be other issues besides 
                            # data outside of the lat-lon range.
                        print("trimmed swath",trimmed_swath)
                        if trimmed_swath is None:
                            print(f"No valid data found in lats {lat_lims}")     
                            pass
                        else:
                            trimmed_filename = f"{remote_file[:-3]}{trim_suffix}.nc"
                            trimmed_swath.to_netcdf(f"{save_path}/{trimmed_filename}")

                    if os.path.isfile(temp_file):
                        try:
                            os.remove(temp_file)
                        except:
                            print(f"Having trouble removing {temp_file} from local directory")
                        
                    print(f"Downloaded and trimmed {remote_file}")
            except Exception as e:
                print(f"Failed to download {remote_file}. Error: {e}")
                print(traceback.print_exc(file=sys.stdout))

    return
    

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function: clean_incomplete_files
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def clean_incomplete_files(path, size=0.05):
    """
    Deletes incomplete files in a directory based on file size.

    Parameters
    ----------
    path : str
        Path to the directory containing files to check.
    size : float
        Minimum file size (in MB) to consider as valid.
        Setting to 50kb, since some of the Expert data is 
        quite lightweight.

    Notes
    -----
    - Deletes files smaller than the specified size.
    """
    files = os.listdir(path)
    for file in files:
        file_stats = os.stat(f"{path}/{file}")
        if file_stats.st_size / (1024 * 1024) < size:  # Convert bytes to MB
            os.remove(f"{path}/{file}")
            print(f"Deleted possible incomplete file: {path}/{file}")

    return



