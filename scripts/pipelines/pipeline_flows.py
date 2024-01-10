import yaml
from dotenv import load_dotenv
import os
import sys
import time
import multiprocessing

from .pipeline_load_localPG import import_to_local_db
from .pipeline_transform_vrt_gdal import geofilter_paths_list, gdal_build_vrt, absolute_file_paths, transform_raster, categorize_aspect, create_vrt_ovr_flow, split_list
from . import model_pipeline

sys.path.append("/home/nikola/4_north_america/scripts/")
import settings

# load configuration
config = settings.get_config()
conn_parameters = settings.get_conn_parameters()
schema = settings.get_schema()
database_name = conn_parameters["database"]

# TODO: DRY the load_countries_boundaries_flow, load_continents_boundaries_flow, load_ESA_dem_grid_flow
# TODO: Convert to OOP - convert functions into methods and clases, use polymorphism for similar functions

def download_tree_cover_density_flow():
    tree_cover_density_downloader = model_pipeline.GisDataDownloader(config["tree_cover_density_dir"])
    return tree_cover_density_downloader.download_tree_cover_density()


def load_countries_boundaries_flow():
    """
    load downloaded files for world boundaries into local pgdb
    """
    world_countries_boundaries_file = config["world_countries_boundaries_file"]
    table_name = config["world_countries_boundaries_table"]
    import_to_local_db(schema, table_name, world_countries_boundaries_file, database_name, "geopackage")
    return print(f"task {load_countries_boundaries_flow.__name__} done!")


def load_continents_boundaries_flow():
    """
    load downloaded files for boundaries of continents into local pgdb
    """
    world_continents_boundaries_file = config["world_continents_boundaries_file"]
    table_name = config["world_continents_boundaries_table"]
    import_to_local_db(schema, table_name, world_continents_boundaries_file, database_name, "geopackage")
    return print(f"task {load_continents_boundaries_flow.__name__} done!")


def load_ESA_dem_grid_flow():
    """
    load grid file for ESA global DEM tiles into local pgdb
    """
    esa_global_dem_grid_file = config["esa_global_dem_grid_file"]
    table_name = config["esa_global_dem_grid_table"]
    import_to_local_db(schema, table_name, esa_global_dem_grid_file, database_name, "geopackage")
    return print(f"task {load_ESA_dem_grid_flow.__name__} done!")


def load_aspect_flow():
    """
    load aspect file computed from esa global dem 90m into local pgb
    """
    load_file = config["NA_aspect_categorized"]
    table_name = config["NA_aspect_categorized_table"]
    aspect_loader = model_pipeline.DataLoader(table_name, load_file, database_name, "raster")
    return aspect_loader.load_local()


def load_slope_flow():
    """
    load slope file computed from esa global dem 90m into local pgb
    """
    load_file = config["NA_slope_rescaled"]
    table_name = config["NA_slope_rescaled"]
    slope_loader = model_pipeline.DataLoader(table_name, load_file, database_name, "raster")
    return slope_loader.load_local()


# transform flows: vrt/ovr (usgs10, esa_global_dem_90, esa_na_dem_30, esa_na_dem_90, esa_world_cover), slope, aspect, geomorphon

def usa_dem_vrt_ovr_flow():
    return create_vrt_ovr_flow("usgs_dem_dir", "usa_dem_vrt")

def transform_esa_usa_dem_vrt_flow():
    files_dir = config["esa_global_dem_90_dir"]
    output_vrt = config["esa_usa_dem_90_vrt"]
    esa_usa_dem_vrt_90_transformer = model_pipeline.DataTransformer(files_dir)
    return esa_usa_dem_vrt_90_transformer.build_vrt(output_vrt, extent="United States")

def transform_ESA_NA_dem_90_vrt_flow():
    """
    Create virtual raster (vrt file) for downloaded DEM 90m resolution dataset from ESA. 
    The dataset is filtered using `filter_paths_list` so that it covers the boundaries of North America/usa.
    """
    files_dir = config["esa_global_dem_90_dir"]
    output_vrt = config["esa_na_dem_90_vrt"]
    esa_NA_dem_vrt_90_transformer = model_pipeline.DataTransformer(files_dir)
    return esa_NA_dem_vrt_90_transformer.build_vrt(output_vrt, extent="North American")


def transform_ESA_NA_dem_30_vrt_flow():
    """
    Create virtual raster (vrt file) for downloaded DEM 30m resolution dataset from ESA. 
    The dataset was filtered using `filter_paths_list` in download step so that it covers the boundaries of North America/usa.
    """
    input_dir = config["esa_global_dem_30_dir"]
    output_vrt = config["esa_global_dem_30_vrt"]
    esa_NA_dem_vrt_30_transformer = model_pipeline.DataTransformer(input_dir)
    return esa_NA_dem_vrt_30_transformer.build_vrt(output_vrt)


def transform_ESA_global_dem_90_vrt_flow():
    """
    Create virtual raster (vrt file) for downloaded DEM 90m resolution dataset from ESA global dem (global boundary).  
    """
    input_dir = config["esa_global_dem_90_dir"]
    output_vrt = config["esa_global_dem_90_vrt"]
    esa_global_dem_vrt_90_transformer = model_pipeline.DataTransformer(input_dir)
    return esa_global_dem_vrt_90_transformer.build_vrt(output_vrt)


def transform_ESA_world_cover_vrt_flow():
    """
    create vrt for esa world cover tiles (global extent)
    """
    esa_world_cover_dir = config["esa_world_cover_dir"]
    esa_world_cover_vrt = config["esa_world_cover_vrt"]
    esa_world_cover_transformer = model_pipeline.DataTransformer(esa_world_cover_dir)
    return esa_world_cover_transformer.build_vrt(esa_world_cover_vrt)


# slope, aspect
def transform_slope_flow():
    """    
    Compute slope using the gdaldem command from shell
    """
    input_raster = config["esa_usa_dem_90_vrt"]
    output_raster = config["NA_slope"]
    transform_raster(input_raster, output_raster, transform_type="slope")
    return print(f"task transform_slope_flow done!")


def transform_aspect_flow():
    """    
    Compute aspect using the gdaldem command from shell. 
    Convert continuous raster values into 8 discrete categories afterwards using `categorize_slope` function.
    """
    input_raster = config["esa_na_dem_90_vrt"]
    aspect_path = config["NA_aspect"]
    transform_raster(input_raster, aspect_path, transform_type="aspect")
    output_category_path = config["NA_aspect_categorized"]
    categorize_aspect(aspect_path, output_category_path)
    return print(f"task transform_aspect_flow done!")


def transform_geomorphon_flow(process_partition = "whole or chunky"):
    """
    Transform raster tiles using geomorphon algorithm from grass gis.
    """
    input_raster = config["esa_usa_dem_90_vrt"]
    geomorphon_file = config["NA_geomorphon"]
    geomorphon_transformer = model_pipeline.DataTransformer(input_raster)
    return geomorphon_transformer.geomorphon(geomorphon_file)


def transform_geomorphon_chunky():
    os.makedirs(config["NA_geomorphon_dir"], exist_ok=True)
    files_list = geofilter_paths_list(config["esa_global_dem_90_dir"])
    files_list = split_list(files_list, 3)
    part_1 = files_list[0]
    part_2 = files_list[1]
    part_3 = files_list[2]
    geomorphon_transformer = model_pipeline.DataTransformer(files_list)
    return geomorphon_transformer.geomorphon_chunks()

    # files_list = "/mnt/volume-nbg1-1/shared/nikola/ESA_global_dem_90m/Copernicus_DSM_30_N45_00_E015_00_DEM.tif"

def geomorphon_process_files(file_list):
    """
    Process a list of files using the DataTransformer. 
    Transform mode is chunky - dataset is transformed chunk by chunk, not whole.

    Args:
        file_list (list): List of file paths to be processed.

    Returns:
        list: Result of processing the file list using the DataTransformer.
    """
    geomorphon_transformer = model_pipeline.DataTransformer(file_list)
    return geomorphon_transformer.geomorphon_chunks()


def transform_geomorphon_multiprocess_flow():
    """
    Transform dem data with geomorphon algorithm using parallel processing.

    This function reads a list of file paths from the ESA global DEM directory,
    splits it into three parts, and processes each part in parallel using
    multiprocessing. The processing involves creating a DataTransformer instance
    for each part and calling the geomorphon_chunks method to transform the data.

    Returns:
        list: List containing the results of processing each part of the file list.
    """
    os.makedirs(config["NA_geomorphon_dir"], exist_ok=True)
    
    files_list = geofilter_paths_list(config["esa_global_dem_90_dir"])

    # now, if a process has stopped at any moment, you will want to set this to true so that the proces doesn't have to be repeated again
    skip_created_files = True
    if skip_created_files == True:
        already_created_files = absolute_file_paths(config["NA_geomorphon_dir"], ".tif")
        processed_files_stripped = []
        for file in already_created_files:
            stripped_file = os.path.basename(file).replace('_geomorphon.tif', '.tif')
            stripped_file = os.path.join(config["esa_global_dem_90_dir"], stripped_file)
            processed_files_stripped.append(stripped_file)

        files_not_processed = [file for file in files_list if file not in processed_files_stripped]
        files_list = files_not_processed

    files_list = split_list(files_list, 3)
    # Create a multiprocessing pool with the number of processes - one per list
    pool = multiprocessing.Pool(processes=3)
    
    # Use the pool.map to apply the process_files function to each part of the list
    result = pool.map(geomorphon_process_files, files_list)
    
    # Close the pool to free up resources
    pool.close()
    pool.join()
    
    return result


def transform_tree_cover_density_vrt_flow():
    """
    Create virtual raster (vrt file) for downloaded tree cover density tiles of 30m resolution. 
    """
    input_dir = config["tree_cover_density_dir"]
    output_vrt = config["tree_cover_density_vrt"]
    TCD_transformer = model_pipeline.DataTransformer(input_dir)
    return TCD_transformer.build_vrt(output_vrt)

def transform_geomorphon_vrt_flow():
    """
    Create virtual raster (vrt file) for transformed geomorphon tiles (90m resolution). 
    """
    input_dir = config["NA_geomorphon_dir"]
    output_vrt = config["NA_geomorphon_vrt"]
    geomorphon_transformer = model_pipeline.DataTransformer(input_dir)
    return geomorphon_transformer.build_vrt(output_vrt)


def coastal_flooding_flow():
    """
    Create predicted areas where coastal flooding could occur using DEM tiles.
    """
    os.makedirs(config["sea_level_rise_tmp"], exist_ok=True)
    files_list = geofilter_paths_list(config["esa_global_dem_30_dir"])     
    
    coastal_flooding_transformer = model_pipeline.DataTransformer(files_list)
    return coastal_flooding_transformer.coastal_flooding()    
