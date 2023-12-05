import yaml
from dotenv import load_dotenv
import os
import sys
import time

from .pipeline_load_localPG import import_to_local_db
from .pipeline_transform_vrt_gdal import geofilter_paths_list, gdal_build_vrt, absolute_file_paths, transform_raster, categorize_aspect, create_vrt_ovr_flow
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


def transform_geomorphon_flow():
    """
    Transform raster tiles using geomorphon algorithm from grass gis.
    """
    # input_raster = config["esa_usa_dem_90_vrt"]
    input_raster = "/mnt/volume-nbg1-1/shared/nikola/ESA_global_dem_90m/Copernicus_DSM_30_N45_00_E015_00_DEM.tif"
    # geomorphon_file = config["NA_geomorphon"]
    geomorphon_file = config["NA_geomorphon"]+"_45_15.tif"
    geomorphon_transformer = model_pipeline.DataTransformer(input_raster)
    return geomorphon_transformer.geomorphon(geomorphon_file)


def transform_geomorphon_chunky():
    os.makedirs(config["NA_geomorphon_dir"], exist_ok=True)
    files_list = geofilter_paths_list(config["esa_global_dem_90_dir"])
    geomorphon_transformer = model_pipeline.DataTransformer(files_list)
    return geomorphon_transformer.geomorphon_class()


def transform_tree_cover_density_vrt_flow():
    """
    Create virtual raster (vrt file) for downloaded tree cover density tiles of 30m resolution. 
    """
    input_dir = config["tree_cover_density_dir"]
    output_vrt = config["tree_cover_density_vrt"]
    TCD_transformer = model_pipeline.DataTransformer(input_dir)
    return TCD_transformer.build_vrt(output_vrt)