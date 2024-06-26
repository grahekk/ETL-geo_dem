import yaml
from dotenv import load_dotenv
import os
import sys
import time
import multiprocessing
import tempfile

from .pipeline_load_localPG import import_to_local_db
from .pipeline_transform_vrt_gdal import geofilter_paths_list, gdal_build_vrt, absolute_file_paths, transform_raster, categorize_aspect, create_vrt_ovr_flow, split_list, geocell_regex_match
from . import model_pipeline
from .pipeline_transform_sea_level import sea_level_precheck, create_file_name_from_geocellid, coastal_flooding_pixel_prediction

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
    esa_NA_dem_vrt_90_transformer = model_pipeline.DataTransformer(files_dir, output_vrt)
    return esa_NA_dem_vrt_90_transformer.build_vrt(extent="North America")


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
    esa_world_cover_transformer = model_pipeline.DataTransformer(esa_world_cover_dir, esa_world_cover_vrt)
    return esa_world_cover_transformer.build_vrt()


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


def eu_resample_1km_flow():
    """
    Transformation flow for geomorphon for europe for resolution of 1km.
    1. resamples the 90m dem to 1km resolution.
    2. calculates geomorphon
    """

    # files_list = geofilter_paths_list(config["esa_global_dem_90_dir"], by="extent", extent = (-25, 35, 35, 72))

    EU_resample_1km_tif = config["EU_resample_1km_tif"]
    EU_resample_1km_vrt = config["EU_resample_1km_vrt"]

    # with tempfile.TemporaryDirectory() as tmpdirname:
    #     resample_transformer = model_pipeline.DataTransformer(files_list, tmpdirname)
    #     resample_transformer.resample_tif()
        
    #     for fname in os.listdir(tmpdirname):
    #         if "clipped" in fname:
    #             os.remove(os.path.join(tmpdirname, fname))
    #             print(f"{fname} removed!")


    #     files_list = absolute_file_paths(tmpdirname, extension=".tif")
    #     vrt_transformer = model_pipeline.DataTransformer(tmpdirname, EU_resample_1km_vrt)

    #     vrt_transformer.build_vrt(EU_resample_1km_vrt)
    #     vrt_transformer.output = EU_resample_1km_tif
    #     vrt_transformer.merge_vrt()

    EU_slope_1km_tif = config["EU_slope_1km_tif"]
    # slope_transformer = model_pipeline.DataTransformer(EU_resample_1km_tif, EU_slope_1km_tif)
    # slope_transformer.slope()
    
    EU_geomorphon_1km_tif = config["EU_geomorphon_1km_tif"]
    geomorphon_transformer = model_pipeline.DataTransformer(EU_resample_1km_tif, EU_geomorphon_1km_tif)
    geomorphon_transformer.geomorphon(search=17, skip=3, flat=1, dist=1)
    print(f"geomorphon done {EU_geomorphon_1km_tif}")
    

def transform_tree_cover_density_vrt_flow():
    """
    Create virtual raster (vrt file) for downloaded tree cover density tiles of 30m resolution. 
    """
    input_dir = config["tree_cover_density_dir"]
    output_vrt = config["tree_cover_density_vrt"]
    TCD_transformer = model_pipeline.DataTransformer(input_dir, output_vrt)
    TCD_transformer.build_vrt()
    

def transform_geomorphon_vrt_flow():
    """
    Create virtual raster (vrt file) for transformed geomorphon tiles (90m resolution). 
    """
    input_dir = config["NA_geomorphon_dir"]
    output_vrt = config["NA_geomorphon_vrt"]
    geomorphon_transformer = model_pipeline.DataTransformer(input_dir)
    return geomorphon_transformer.build_vrt(output_vrt)

def transform_coastal_flooding_vrt():
    """
    take dem tiles with potential flood and create vrt for them
    """
    tiles_path = config["coastal_flooding_tiles"]
    with open(tiles_path, "r") as file: #prevoiusly saved file from coastal_flooding precheck
        lines = file.read()
        tiles = lines.split(", ")

    output_vrt = config["coastal_flooding_vrt"]
    gdal_build_vrt(tiles, output_vrt)
    return

def coastal_flooding_90m_pixels_flow():
    """
    Create predicted areas where coastal flooding could occur using DEM 90m tiles.
    """
    coastal_flooding_transformer = model_pipeline.DataTransformer(config["esa_global_dem_90_dir"], config["data_folder"])
    coastal_flooding_transformer.coastal_flooding_tiles()
    return

def coastal_flooding_90m_pixels_multiprocessing_flow():
    """
    The same function as `coastal_flooding_90m_pixels_flow`, but with multiprocessing. 
    Create predicted areas where coastal flooding could occur using DEM 90m tiles.
    """
    coastal_flooding_transformer = model_pipeline.DataTransformer(config["esa_global_dem_90_dir"], config["data_folder"])
    coastal_flooding_transformer.coastal_flooding_tiles_multiprocessing()
    return

def coastal_flooding_90m_Europe_flow():
    """
    Create predicted areas where coastal flooding could occur using DEM 90m tiles.
    """
    coastal_flooding_transformer = model_pipeline.DataTransformer(config["esa_global_dem_90_dir"], config["coastal_flooding_EU_90m"])
    for i in range(0,2):
        coastal_flooding_transformer.coastal_flooding_tiles()
    return

def coastal_flooding_30m_NA_flow():
    """
    Create predicted areas where coastal flooding could occur using DEM 30m tiles for North America.
    """
    coastal_flooding_transformer = model_pipeline.DataTransformer(config["esa_na_dem_30_dir"], config["NA_coastal_flooding_30m_dir"])
    for i in range(0,2):
        coastal_flooding_transformer.coastal_flooding_tiles()
    return

def coastal_flooding_30m_NA_multiprocessing_flow():
    """
    Create predicted areas where coastal flooding could occur using DEM 30m tiles for North America.
    """
    coastal_flooding_transformer = model_pipeline.DataTransformer(config["esa_na_dem_30_dir"], config["NA_coastal_flooding_30m_dir"])
    for i in range(0,2):
        coastal_flooding_transformer.coastal_flooding_tiles_multiprocessing()
    return

def coastal_flooding_10m_USA_flow():
    """
    Create predicted areas where coastal flooding could occur using DEM 10m tiles for USA from USGS.
    """
    coastal_flooding_transformer = model_pipeline.DataTransformer(config["usgs_dem_dir"], config["USA_coastal_flooding_10m_dir"])
    for i in range(0,2):
        coastal_flooding_transformer.coastal_flooding_tiles()
    return

def coastal_flooding_10m_USA_multiprocessing_flow():
    """
    Create predicted areas where coastal flooding could occur using DEM 30m tiles for North America.
    """
    coastal_flooding_transformer = model_pipeline.DataTransformer(config["usgs_dem_dir"], config["USA_coastal_flooding_10m_dir"])
    for i in range(0,2):
        coastal_flooding_transformer.coastal_flooding_tiles_multiprocessing()
    return

def coastal_flooding_10m_Europe_flow():
    """
    Create predicted areas where coastal flooding could occur using DEM 10m tiles for Europe from ESA.
    """
    coastal_flooding_transformer = model_pipeline.DataTransformer(config["esa_eu_dem_10_dir"], config["EU_coastal_flooding_10m_dir"])
    for i in range(0,2):
        coastal_flooding_transformer.coastal_flooding_tiles()
    return

def tree_cover_filter_flow():
    """
    flow for filtering forests out of world cover dataset
    """
    tree_cover_dir = config["NA_tree_cover_10m_dir"]
    os.makedirs(tree_cover_dir, exist_ok=True)
    filter_dataset = config["esa_world_cover_dir"]

    filter_transformer = model_pipeline.DataTransformer(filter_dataset, tree_cover_dir)
    return filter_transformer.filter_tif(feature = "Forest")

def shrubland_filter_flow():
    """
    flow for filtering shrubland out of world cover dataset
    """
    shrubland_dir = config["NA_shrubs_10m_dir"]
    os.makedirs(shrubland_dir, exist_ok=True)
    filter_dataset = config["esa_world_cover_dir"]

    filter_transformer = model_pipeline.DataTransformer(filter_dataset, shrubland_dir)
    return filter_transformer.filter_tif(feature = "Shrubland")

def grassland_filter_flow():
    """
    flow for filtering shrubland out of world cover dataset
    """
    grassland_dir = config["NA_grassland_10m_dir"]
    os.makedirs(grassland_dir, exist_ok=True)
    filter_dataset = config["esa_world_cover_dir"]

    filter_transformer = model_pipeline.DataTransformer(filter_dataset, grassland_dir)
    return filter_transformer.filter_tif(feature = "Grassland")

def hawaii_repair_flow():
    """
    flow for adding geomorphon, slope and aspect data for Hawaii
    """
    hawaii_geocells = config["NA_hawai_geocells"].split(", ")

    # first vrt
    # then aspect, slope, geomorphon
    # then new geomorphon vrt
    # then new geomorphon out of geomorphon vrt
    # then merge aspect, merge slope
    esa_global_dem_90_dir = config["esa_global_dem_90_dir"]
    files_list = absolute_file_paths(esa_global_dem_90_dir, extension='.tif')
    hawaii_files = geocell_regex_match(files=files_list, regex_match = r'_(S|N)(\d+)_00_(W|E)(\d+)', geocell_filter_list=hawaii_geocells)
    print(hawaii_files)

    hawai_vrt = config["hawaii_vrt"]
    hawaii = model_pipeline.DataTransformer(hawaii_files, hawai_vrt)
    hawaii.build_vrt(extent="files_list")
    print("vrt done")
    hawaii.output = config["hawaii_aspect"]
    hawaii.aspect()
    print("aspect done")
    
    hawaii.data = hawai_vrt
    hawaii.output = config["hawaii_slope"]
    hawaii.slope()
    print("Slope done")

    hawaii.data = hawai_vrt
    hawaii.output = config["hawaii_geomorphon"]
    hawaii.geomorphon()
    print("Geomorphon done")

def hawaii_30m_coastal_flooding_flow():
    """
    flow for using only hawaii dem tiles for filling up the coastal flooding dataset
    """
    NA_hawaii_geocells = config["NA_hawai_geocells"].split(", ")

    file_paths = []
    for geocell in NA_hawaii_geocells:
        file_path = create_file_name_from_geocellid(geocell, config["esa_na_dem_30_dir"])
        file_paths.append(file_path)
        
    for i in file_paths:
        coastal_flooding_pixel_prediction(i, output_folder=config["NA_coastal_flooding_30m_dir"])


def load_to_mnt_flow():
    data_loader = model_pipeline.DataLoader(config["EU_coastal_flooding_90m_table"], config["EU_coastal_flooding_90m_shp"], "shapefile")
    data_loader.load_local()

def load_to_aws_flow():
    data_loader = model_pipeline.DataLoader(config["NA_aspect_table"], config["NA_aspect"], "raster")
    data_loader.dump_local_table()
    # data_loader.aws_load()

def load_coastal_flooding_30m_to_aws_flow():
    data_loader = model_pipeline.DataLoader("na_coastal_flooding_30m_4326", config["NA_coastal_flooding_30m_dir"], "shapefile")
    data_loader.aws_iterative_shp_load()

def load_coastal_flooding_10m_EU_to_mnt_flow():
    data_loader = model_pipeline.DataLoader("eu_coastal_flooding_10m_4326", config["EU_coastal_flooding_10m_dir"], "shapefile")
    data_loader.local_iterative_shp_load()

def load_coastal_flooding_10m_EU_to_aws_flow():
    data_loader = model_pipeline.DataLoader("eu_coastal_flooding_10m_4326", config["EU_coastal_flooding_10m_dir"], "shapefile")
    data_loader.aws_iterative_shp_load()

def load_coastal_flooding_10m_USA_to_aws_flow():
    data_loader = model_pipeline.DataLoader("usa_coastal_flooding_10m_4329", config["USA_coastal_flooding_10m_dir"], "shapefile")
    data_loader.aws_iterative_shp_load()