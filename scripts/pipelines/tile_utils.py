import os
import sys
import re
import math
from .pipeline_transform_vrt_gdal import gdal_build_vrt

import settings

config = settings.get_config()

USGS_dem_regex = settings.get_config_params().USGS_dem_files_regex_match
ESA_dem_files_regex_match = settings.get_config_params().ESA_dem_files_regex_match

coastlines = config["NA_coastlines"]
regex_match = ESA_dem_files_regex_match
regex_match = USGS_dem_regex

def generate_neighbouring_coordinates(center: tuple, delta=0.15):
    """
    Generate coordinates for a given center point, considering eight directions.

    Parameters:
    - center (tuple): Tuple of two numbers (lat, lon) representing the central point.
    - delta (float): Distance from the central point to the boundary in each direction.

    Returns:
    - list of tuples: List containing four numbers (x_min, y_min, x_max, y_max) for each direction.
    """
    lat, lon = center
    lat = int(lat)
    lon = int(lon)

    left = (lon - delta, lat, lon, lat + 1)
    down = (lon, lat - delta, lon + 1, lat)
    right = (lon + 1, lat, lon + 1 + delta, lat + 1)
    up = (lon, lat + 1, lon + 1, lat + 1 + +delta)
    down_right = (lon + 1, lat - delta, lon + 1 + delta, lat)
    down_left = (lon - delta, lat - delta, lon, lat)
    up_right = (lon + 1, lat + 1, lon + 1 + delta, lat + 1 + delta)
    up_left = (lon - delta, lat + 1, lon, lat + 1 + delta)

    return left, down, right, up, down_right, down_left, up_right, up_left


def replace_coordinates(input_string, lat_value, lon_value):
    """
    Replace latitude and longitude values in a string based on a specified regex pattern.

    Args:
        input_string (str): The input string containing latitude and longitude placeholders.
        lat_value (str): The new latitude value to replace the matched latitude placeholder.
        lon_value (str): The new longitude value to replace the matched longitude placeholder.

    Returns:
        str: A new string with replaced latitude and longitude values.

    Example:
        >>> input_string = "Copernicus_DSM_03_N45_00_E015_00_DEM"
        >>> lat_value = "1"
        >>> lon_value = "1"
        >>> replace_coordinates(input_string, lat_value, lon_value)
        'Copernicus_DSM_03_N1_00_E1_00_DEM'
    """
    regex_pattern = regex_match

    def replacement(match):
        """
        Replacement function for re.sub.

        Args:
            match (re.Match): The regex match object.

        Returns:
            str: The replacement string with adjusted latitude and longitude values.
        """
        return f"_{match.group(1)}{lat_value:02d}_00_{match.group(3)}{lon_value:03d}"

    result = re.sub(regex_pattern, replacement, input_string)
    return result


def extract_coordinates_from_tile_name(file_name: str):
    """
    Extract latitude and longitude coordinates from a given tile name.

    Parameters:
    - file_name (str): The tile name containing latitude and longitude information.
    - regex_match (str): Regular expression pattern to match latitude and longitude components.
      Default pattern assumes a format like 'Copernicus_DSM_03_N45_00_E015_00_DEM.tif'.

    Returns:
    - tuple: A tuple containing the following components:
        - lat_direction (str): Latitude direction ('S' for South, 'N' for North).
        - lat_value (str): Latitude value extracted from the tile name.
        - lon_direction (str): Longitude direction ('W' for West, 'E' for East).
        - lon_value (str): Longitude value extracted from the tile name.
    """
    file_name = os.path.basename(file_name)
    geocell_id_match = re.search(regex_match, file_name)

    if geocell_id_match:
        lat_direction = geocell_id_match.group(1)
        lat_value = geocell_id_match.group(2)
        lon_direction = geocell_id_match.group(3)
        lon_value = geocell_id_match.group(4)

        return lat_direction, lat_value, lon_direction, lon_value
    

def gdalwarp_string(ul_lon, ul_lat, lr_lon, lr_lat, input_path, output_path):
    return f"gdalwarp -overwrite -r near -of GTiff -te {ul_lon} {ul_lat} {lr_lon} {lr_lat} {input_path} {output_path}"


def expand_tile(input_tile_path: str, tmp_folder:str, expansion = 0.007):
    """
    Expands the geographical coverage of a given input tile by including neighbouring tiles.

    Parameters:
    - input_tile_path (str): Path to the input tile file.
    - tmp_folder (str): Path to the temporary folder where intermediate and output files will be stored.
    - expansion (float, optional): The distance, in geographical degrees, by which to expand the tile.
                                   Defaults to 0.007 degrees.

    Returns:
    - str: Path to the created Virtual Raster Tile (VRT) file representing the expanded tile.
    """
    lat_direction,lat_value,lon_direction,lon_value = extract_coordinates_from_tile_name(input_tile_path)
    central_tile_lat_lon = (lat_value, lon_value)
    neighbouring_coords = generate_neighbouring_coordinates(central_tile_lat_lon, delta=expansion)

    vrt_paths = [input_tile_path]

    # neighbouring tiles clipping part
    for i in range(0, 8):
        ul_lon, ul_lat, lr_lon, lr_lat = neighbouring_coords[i]
        gdalwarp_tile_name = f"gdalwarp_{math.floor(ul_lat)}_{math.floor(ul_lon)}.tif"
        input_warp = replace_coordinates(
            input_tile_path, math.floor(ul_lat), math.floor(ul_lon)
        )

        if os.path.exists(input_warp):
            output_path = os.path.join(tmp_folder, gdalwarp_tile_name)
            command = gdalwarp_string(
                ul_lon, ul_lat, lr_lon, lr_lat, input_warp, output_path
            )
            os.system(command)

            vrt_paths.append(output_path)

    # create vrt
    vrt_file_name = f"expanded_tile{lat_value}_{lon_value}.vrt"
    out_vrt = os.path.join(tmp_folder, vrt_file_name)
    gdal_build_vrt(vrt_paths, out_vrt)

    return out_vrt


def collect_neighbouring_coastal_flood_files(central_file:str, neighbours = 1, include_central = False, coastal_flood_sufix = "coastal_flood"):
    """
    Collects the paths of neighboring coastal flood files given the central file.

    Parameters:
    - central_file (str): Path to the central coastal flood file.
    - neighbours (int): The number of surrounding layers to include.
                       If 1, returns 8 neighbors; if 2, returns 24 neighbors (8 initial + 16 additional)

    Returns:
    - list: A list containing the paths of the 8 neighboring coastal flood files, if they exist.
            The order of files in the list corresponds to the surrounding tiles
            (e.g., top-left, top-center, top-right, left, center, right, bottom-left, bottom-center, bottom-right).
    """
    # naming the files
    if "DSM_10" in central_file:
        dsm = "10"    
    elif "DSM_03" in central_file:
        dsm = "03"
    elif "DSM_30" in central_file:
        dsm = "30"

    if include_central == True:
        neighbour_files = [central_file]
    else:
        neighbour_files = []

    file_name = os.path.basename(central_file)
    tile_folder = os.path.dirname(central_file)
    if "USGS" in central_file:
        usgs_files = os.listdir(tile_folder)
    # find file coords
    geocell_id_match = re.search(regex_match, file_name)

    if geocell_id_match:
        lat_direction = geocell_id_match.group(1)
        lat_value = int(geocell_id_match.group(2))
        lon_direction = geocell_id_match.group(3)
        lon_value = int(geocell_id_match.group(4))

        # find neighbour files
        for i in range(lon_value - 1, lon_value + 2):
            for j in range(lat_value - 1, lat_value + 2):
                # Construct the file name for each surrounding tile

                # england edge case is: W001
                if lon_direction == "W" and lon_value == 1:
                    for k in ["W", "E"]:
                        tile_name = f"Copernicus_DSM_{dsm}_{lat_direction}{j:02d}_00_{k}{i:03d}_00_DEM_{coastal_flood_sufix}.shp"
                        tile_path = os.path.join(tile_folder, tile_name)
                        if os.path.exists(tile_path) and tile_path != central_file:
                            neighbour_files.append(tile_path)
                
                elif lon_direction == "E" and lon_value == 0:
                    for k in ["W", "E"]:
                        tile_name = f"Copernicus_DSM_{dsm}_{lat_direction}{j:02d}_00_{k}{i:03d}_00_DEM_{coastal_flood_sufix}.shp"
                        tile_path = os.path.join(tile_folder, tile_name)
                        if os.path.exists(tile_path) and tile_path != central_file:
                            neighbour_files.append(tile_path)

                else:
                    if "Copernicus" in central_file:
                        tile_name = f"Copernicus_DSM_{dsm}_{lat_direction}{j:02d}_00_{lon_direction}{i:03d}_00_DEM_{coastal_flood_sufix}.shp"
                        tile_path = os.path.join(tile_folder, tile_name)
                        if os.path.exists(tile_path) and tile_path != central_file:
                            neighbour_files.append(tile_path)

                    # TODO: use function for constructing tile names
                    elif "USGS" in central_file:
                        tile_name = f"USGS_13_{lat_direction.lower()}{j:02d}{lon_direction.lower()}{i:03d}_"
                        try:
                            matching_tile = [file for file in usgs_files if tile_name in file and ".shp" in file][0] #this will prolly fail if there is no real file
                            tile_path = os.path.join(tile_folder, matching_tile)
                            if os.path.exists(tile_path) and tile_path != central_file:
                                neighbour_files.append(tile_path)
                        except Exception as e:
                            print(f"skipping {tile_name} because Exception: {e}") # if concatenating nodata file

    return neighbour_files


def build_expanded_tif(input_file, coastal_dataset, tmpdirname):
    """
    Builds an expanded GeoTIFF file from an input file and neighboring coastal flood files.

    Parameters:
    - input_file (str): Path to the input GeoTIFF file.
    - coastal_dataset (list): List of coastal flood files.
    - tmpdirname (str): Path to the temporary directory for storing intermediate files.

    Returns:
    - str: Path to the resulting expanded GeoTIFF file.
    """
    neighbours_list = collect_neighbouring_coastal_flood_files(input_file)
    common_files = common_files_between_lists(neighbours_list, coastal_dataset)

    # build expanded tif/vrt with neighboring dem tiles
    # vrt
    out_vrt = f"{tmpdirname}/{os.path.splitext(os.path.basename(input_file))[0]}_flooding.vrt"
    gdal_build_vrt(common_files, out_vrt)

    # make tif from vrt
    expanded_tif = f"{tmpdirname}/{os.path.splitext(os.path.basename(input_file))[0]}_expanded.tif"
    gdal_translate_cmd = f"gdal_translate -of GTiff -b 1 {out_vrt} {expanded_tif}"
    os.system(gdal_translate_cmd)

    return expanded_tif


def common_files_between_lists(list1, list2, keep_order = True):
    """
    Finds common file names between two lists.

    Parameters:
    - list1 (list): First list of file names.
    - list2 (list): Second list of file names.
    - keep_order (bool): keeps the order of the first list if true

    Returns:
    - list: A list containing file names that are present in both input lists.
    """
    set1 = frozenset(list1)
    set2 = frozenset(list2)
    
    if keep_order == False:
        common_files = list(set1.intersection(set2))
    else:
        common_files = [x for x in list1 if x in set2]    

    return common_files