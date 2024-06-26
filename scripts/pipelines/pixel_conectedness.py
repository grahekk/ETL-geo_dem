import os
import sys
import geopandas as gpd
import scipy
import numpy as np
import rasterio
from rasterio.features import geometry_mask
from scipy import ndimage
from skimage import measure
from shapely.geometry import shape, LineString
from shapely.geometry import shape
import tempfile
import warnings
import shutil
from tqdm import tqdm
from sqlalchemy import create_engine
from multiprocessing import Pool
from tqdm.contrib.concurrent import process_map  # or thread_map

from .pipeline_transform_vrt_gdal import geocellid_from_file_name
from .pipeline_transform_sea_level import (
    basename_withoutext, 
    gpd_extract_features_by_location, 
    clip_vector_dataset, 
    extract_pixels_by_mask, 
    extract_coordinates_from_tile_name,
    collect_neighbouring_coastal_flood_files,
    merge_shapefiles, 
    polygonize_raster, 
    coastal_flooding_rasters_sum,
    gdal_select_pixels,
    has_empty_geometries,
    subtract_one_from_shapefile_attribute,
    get_coastlines)

import settings

config = settings.get_config()
database_url = settings.get_database_url()
engine = create_engine(database_url)
USGS_dem_regex = settings.get_config_params().USGS_dem_files_regex_match
ESA_dem_files_regex_match = settings.get_config_params().ESA_dem_files_regex_match
regex_match = ESA_dem_files_regex_match

# extract connected pixels

def process_region(region_label, labeled_image, mask_raster):
    region_mask = (labeled_image == region_label).astype(np.uint8)
    return region_label if np.any(region_mask * mask_raster) else None

def extract_intersecting_groups_2_multiprocess(main_raster, mask_raster):
    labeled_image = measure.label(main_raster)

    regions = measure.regionprops(labeled_image)
    
    # Using multiprocessing to process regions in parallel
    intersecting_groups = process_map(process_region, [(region.label, labeled_image, mask_raster) for region in regions])

    intersecting_groups = list(filter(None, intersecting_groups))
    result_raster = np.where(np.isin(labeled_image, intersecting_groups), main_raster, 0)

    return result_raster

def extract_intersecting_groups_3(main_raster, mask_raster):
    # label connected components
    labeled_image = measure.label(main_raster)

    # Initialize a list to store intersecting groups
    intersecting_groups = []

    # Iterate through labeled regions and check for intersection
    regions = measure.regionprops(labeled_image)

    for region in tqdm(regions):
        label_value = region.label
        region_mask = (labeled_image == label_value).astype(np.uint8)

        # Check if the region intersects with the rasterized line
        if np.any(region_mask * mask_raster):
            intersecting_groups.append(label_value)

    result_raster = np.where(np.isin(labeled_image, intersecting_groups), main_raster, 0)

    return result_raster

def extract_intersecting_groups_1(main_raster, mask_raster):
    """
    Extract pixel groups from the main raster based on their intersection with a given mask raster.

    Parameters:
    - main_raster (numpy.ndarray): The main raster containing labeled connected components.
    - mask_raster (numpy.ndarray): The rasterized line mask used for intersection.

    Returns:
    - numpy.ndarray: A new raster where pixel values are retained only for groups that intersect with the mask.
    
    Example:
    >>> main_raster = np.array([[1, 1, 0, 0],
                                [1, 0, 2, 0],
                                [0, 0, 2, 2],
                                [3, 3, 3, 0]])
    >>> mask_raster = np.array([[0, 1, 0, 0],
                                [0, 1, 1, 0],
                                [0, 0, 0, 0],
                                [0, 0, 1, 0]])
    >>> result = extract_intersecting_groups(main_raster, mask_raster)
    >>> print(result)
    array([[0, 0, 0, 0],
           [1, 0, 2, 0],
           [0, 0, 2, 2],
           [0, 0, 0, 0]])

    Note:
    - The function uses skimage.measure.label to identify connected components in the main raster.
    - It iterates through labeled regions and checks for intersection with the provided mask raster.
    - The resulting raster retains only pixel values for groups that intersect with the mask.
    """

    # label connected components
    labeled_image = measure.label(main_raster)

    # Initialize a list to store intersecting groups
    intersecting_groups = []

    # Iterate through labeled regions and check for intersection
    regions = measure.regionprops(labeled_image)

    for region in tqdm(regions):
        label_value = region.label
        region_mask = (labeled_image == label_value).astype(np.uint8)

        # Check if the region intersects with the rasterized line
        if np.any(region_mask * mask_raster):
            intersecting_groups.append(label_value)

    result_raster = np.where(np.isin(labeled_image, intersecting_groups), main_raster, 0)

    return result_raster

def get_intersection_between_raster_and_geometry(raster_path, output_tif_path, geometry_path):
    """
    Extract pixel groups from a raster based on their intersection with a vector geometry.

    Parameters:
    - raster_path (str): The file path to the input raster TIFF file.
    - output_tif_path (str): The file path to save the output GeoTIFF file.
    - geometry_path (str): The file path to the vector geometry (e.g., Shapefile, GeoJSON).

    Returns:
    None

    Example:
    >>> raster_path = "path/to/your/raster.tif"
    >>> output_tif_path = "path/to/your/output.tif"
    >>> geometry_path = "path/to/your/geometry.shp"
    >>> get_intersection_between_raster_and_geometry(raster_path, output_tif_path, geometry_path)

    Note:
    - The function reads the raster data and vector geometry using rasterio and geopandas, respectively.
    - It creates a mask based on the vector geometry using rasterio.features.geometry_mask.
    - The extract_intersecting_groups function is then used to obtain pixel groups that intersect with the mask.
    - The result is saved to a new GeoTIFF file with the specified output path.

    """
    with rasterio.open(raster_path) as src:
        # Read geometry from the vector file
        if "SELECT" in geometry_path:
            gdf = gpd.read_postgis(geometry_path, engine, geom_col="geom")
            gdf.to_file(raster_path.replace(".tif", "_queried_coastlines.shp"))
        # print(f"reference_gdf saved temporarily")
        else:
            gdf = gpd.read_file(geometry_path)
        # Create a mask based on the geometry
        mask = geometry_mask(gdf.geometry, out_shape=src.shape, transform=src.transform, invert=True)
        raster_data = src.read(1)

        # for debugging purposes, save the intermediate file
        intermedium_path = output_tif_path.replace(".tif", "_intermedium.tif")
        with rasterio.open(intermedium_path, 'w', driver='GTiff', 
                           width=src.width, height=src.height,
                           count=1, dtype=np.int8,
                           crs=src.crs, transform=src.transform) as dst:
            # Write the labeled data to the output raster
            dst.write(raster_data, 1)

        result = extract_intersecting_groups_1(raster_data, mask)

        # Create the output GeoTIFF file
        with rasterio.open(output_tif_path, 'w', driver='GTiff', 
                           width=src.width, height=src.height,
                           count=1, dtype=np.int8,
                           crs=src.crs, transform=src.transform) as dst:
            # Write the labeled data to the output raster
            dst.write(result, 1)



def coastal_flooding_pixel_prediction(input_file: str, sea_level = 10, coast_distance = 0, skip_if_exists = False, output_folder = None):
    """
    Main function that calculates potential rising sea levels from DEM file based on pixel selection principle.
    It loops over levels up to given maximum level in meters (Loop imitates flooding).

    Algorithm is run in steps:
        - 1. pixels are extracted for sea levels from 0 to 10 from dem.
        - 2. pixels are polygonized for topological analysis
        - 3. vectors intersections with coastline checked
        - 3.a Vector intersection with neighbour coastal flooding files, if there are some, check
        - 4. pixels for difference are extracted by vector mask
        - 5. difference between larger and smaller layer - pixel subtraction
        - 6. difference layers are merged - pixel addition
        - 7. pixels are polygonized for final product
        - 8. all columns in final product are subtracted for 1

    Note: all files except the final `output_file` are in stored in tmp folder during execution.

    Parameters:
        input_file(str): path to dem tile for which coastal flooding is calculated
        sea_level(int): max sea level wanted for calculation
    """
    sea_level = sea_level + 1  # otherwise 10 will not be selected because of range()
    original_file = input_file
    data_folder = config["USA_coastal_flooding_10m_dir"]

    output_file = f"{data_folder}/{basename_withoutext(original_file)}_coastal_flood.shp"

    # if skip_if_exists == True:
    #     if os.path.exists(output_file):
    #         logging.info(f"File {output_file} existed, so it was skipped")
    #         return
        
    # neighbouring_coastal_files = collect_neighbouring_coastal_flood_files(output_file, neighbours=coast_distance)
    
    lat_direction, lat, lon_direction, lon = extract_coordinates_from_tile_name(input_file, regex_match)
    lat = int(lat)
    lon = int(lon)

    with tempfile.TemporaryDirectory() as tmpdirname:
        # build_expanded_tif(input_file, coastal_dataset, tmpdirname)

        # # build neighbouring vector layer
        # neighbours_str = " ".join(neighbouring_coastal_files)
        # neighbours_shp = f"{tmpdirname}/{basename_withoutext(original_file)}_neighbours.shp"

        # # merge neighbouring shapefiles, -single option is for merging layers into a single shp
        # merge_neighbours_cmd = f"ogrmerge.py -o {neighbours_shp} {neighbours_str} -single -progress "
        # os.system(merge_neighbours_cmd)

        # # clip neighbours, only tiny bit of polygons is needed
        # neighbours_shp_clipped = neighbours_shp.replace(".shp", "_clipped.shp")
        neighbours_shp_clipped = None
        # try:
        #     if lat_direction == "N" and lon_direction == "W": # Americas
        #         clip_vector_dataset(neighbours_shp, neighbours_shp_clipped, -lon+1.001, lat-0.001, -lon - 0.001, lat + 1.001)
        #     elif lat_direction == "N" and lon_direction == "E": # Europe
        #         clip_vector_dataset(neighbours_shp, neighbours_shp_clipped, lon+1.001, lat-0.001, lon - 0.001, lat + 1.001)

        # except Exception:
        #     neighbours_shp_clipped = None

        # pre select
        # we change input file to the expanded_tif - or leave it with gdal_select pixels(input_file)
        preselected_file = (f"{tmpdirname}/{basename_withoutext(input_file)}_preselect.tif")
        gdal_select_pixels(input_file, preselected_file, sea_level, replace_pixels=False)

        sea_levels_simulation = (f"{tmpdirname}/{basename_withoutext(input_file)}_sea_levels_simulation.tif")

        for level in range(0, sea_level, 1):
            # selecting by height
            selected_sea_level_file = (f"{tmpdirname}/{basename_withoutext(input_file)}_sea_level_{level}.tif")
            gdal_select_pixels(preselected_file, selected_sea_level_file, level)

            # pixel connectivity
            polygonized_pixels_file = (f"{tmpdirname}/{basename_withoutext(input_file)}_polygonized_{level}.shp")
            # polygonize_raster(selected_sea_level_file, polygonized_pixels_file)

            # spatial relation with sea or spatial relation with flood
            intersected_features_file = f"{tmpdirname}/{basename_withoutext(input_file)}_sea_level_intersected_{level}.tif"
            
            if level == 0 and coast_distance == 0:
                geocellid = geocellid_from_file_name(input_file, regex_match)
                coastlines_query = get_coastlines(geocellid)
                
                # gpd_extract_features_by_location(polygonized_pixels_file, coastlines_query, intersected_features_file)

                get_intersection_between_raster_and_geometry(selected_sea_level_file, intersected_features_file, coastlines_query)

            else:
                try:
                    gpd_extract_features_by_location(polygonized_pixels_file, 
                                        reference_file, 
                                        intersected_features_file,
                                        condition = "rename")
                except UnboundLocalError: # if reference file is not yet initialised  (level is 0)
                    if neighbours_shp_clipped == None:
                        return
                    else:
                        gpd_extract_features_by_location(polygonized_pixels_file, 
                                                    neighbours_shp_clipped, 
                                                    intersected_features_file, 
                                                    condition = level, 
                                                    mode = "intersects", 
                                                    join = "inner")
                        neighbours_checked = True


            # intersection with neighbouring flooding files
            if neighbouring_coastal_files:
                try:
                    if neighbours_checked == True:
                        del neighbours_checked

                except UnboundLocalError: # if neighbours haven't yet been tested
                    intersected_neighbours_file = f"{tmpdirname}/{basename_withoutext(input_file)}_neighbours_intersected_{level}.shp"
                    gpd_extract_features_by_location(polygonized_pixels_file, 
                                                    neighbours_shp_clipped, 
                                                    intersected_neighbours_file, 
                                                    condition = level, 
                                                    mode = "intersects", 
                                                    join = "inner")
                    if has_empty_geometries(intersected_neighbours_file):
                        pass
                    else:
                        # merge two intersections
                        merge_shapefiles(intersected_features_file, intersected_neighbours_file, intersected_features_file, "geometry")

            # pixel connectivity and intersection overcome - clip the raster by mask
            intersected_tif = f"{tmpdirname}/{basename_withoutext(input_file)}_sea_level_intersected_{level}.tif"
            try:
                extract_pixels_by_mask(
                    selected_sea_level_file,
                    intersected_features_file,
                    intersected_tif,
                    mode="clip",
                    level=level,
                )
            except ValueError as e: # if layer is empty (no selection), this error is raised
                print(f"The extract_pixels_by_mask() for level: {level} of a {basename_withoutext(original_file)} couldn't be processed because of {e}")
                reference_file = intersected_features_file
                continue

            # this part calculates the difference
            if level > 0 and not has_empty_geometries(reference_file):
                # difference with rasters
                difference_tif = f"{tmpdirname}/{basename_withoutext(input_file)}_sea_level_difference_{level}.tif"
                try:
                    extract_pixels_by_mask(
                        intersected_tif,
                        reference_file,
                        difference_tif,
                        mode="difference",
                        level=level,
                    )
                except Exception as e:
                    print(f"The difference for level: {level} of a {basename_withoutext(original_file)} couldn't be processed because of {e}")
                    difference_tif = selected_sea_level_file

                # polygonize - difference shall serve as reference file in next iter
                difference_polygonized = f"{tmpdirname}/{basename_withoutext(input_file)}_difference_polygonized_{level}.shp"
                polygonize_raster(difference_tif, difference_polygonized)

                # reference file for difference extraction *after* the differences were extracted
                reference_file = intersected_features_file

                # edge cases - if level 0 was skipped there wont be `base_layer`
                try:
                    if base_layer == None:
                        base_layer = difference_tif

                except UnboundLocalError:
                    base_layer = difference_tif

                # add a tif layer to one tif file
                coastal_flooding_rasters_sum(base_layer, difference_tif, sea_levels_simulation, level)
                base_layer = sea_levels_simulation

            else:
                # in case that level is 0 (or other beginning level), existing file is used just for invert
                reference_file = intersected_features_file  # vector mask file
                difference_tif = f"{tmpdirname}/{basename_withoutext(input_file)}_sea_level_difference_{level}.tif"
                try:
                    extract_pixels_by_mask(
                        intersected_tif,
                        reference_file,
                        difference_tif,
                        mode="invert",
                        level=level,
                    )
                except:
                    difference_tif = selected_sea_level_file

                # this base layers is the first to be added to main file
                base_layer = difference_tif

        # copy final tif to data folder, if tile was empty raise warning and exit
        try:
            shutil.copyfile(sea_levels_simulation, output_file.replace(".shp", ".tif"))
        except FileNotFoundError as e:
            warnings.warn(f"File {sea_levels_simulation} is empty and therefore not complete")
            return

        # and convert the result
        polygonized_simulation = (f"{tmpdirname}/{basename_withoutext(input_file)}_coastal_flooding.shp")
        polygonize_raster(sea_levels_simulation, polygonized_simulation)
        
        # clip back to tile size
        clipped_dataset = (f"{tmpdirname}/{basename_withoutext(input_file)}_stimulation_clipped.shp")
        clip_vector_dataset(polygonized_simulation, clipped_dataset, lon, lat, lon + 1, lat + 1)

        # last subtraction
        subtract_one_from_shapefile_attribute(polygonized_simulation, "level", output_file)

    return print(output_file, " is done")


raster_path = "/tmp/tmpwhycykz2/Copernicus_DSM_30_N45_00_E014_00_DEM_sea_level_0.tif"
output_tif_path = "/tmp/tmpwhycykz2/Intersecting_pixels_for_level_0.tif"
vector_path = "/tmp/tmpwhycykz2/Copernicus_DSM_30_N45_00_E014_00_DEM_polygonized_0_queried_coastlines.shp"


# Example usage:
# get_intersection_between_raster_and_geometry(raster_path, output_tif_path, vector_path)
