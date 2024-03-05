import psycopg2
from psycopg2 import sql
import sys
import os
from tqdm import tqdm
import subprocess
import numpy as np
import time
import warnings
import fiona
import rasterio
from rasterio.features import geometry_mask
from skimage import measure

import rasterio.mask
import scipy
import geopandas as gpd
import tempfile
import re
import pandas as pd
import math
import multiprocessing
import concurrent.futures
from shapely.ops import unary_union
from sqlalchemy import create_engine
import shutil
import logging

from .pipeline_transform_vrt_gdal import gdal_build_vrt, geocell_regex_match, geocellid_from_file_name, geofilter_paths_list
from .model_data import get_geocellid, get_product_name

# Define the database connection parameters
sys.path.append("/home/nikola/4_north_america/GeoDataPump/scripts")
import settings

config = settings.get_config()
Configuration = settings.get_config_params()
conn_parameters = settings.get_conn_parameters()
schema = settings.get_schema()

database_name = conn_parameters["database"]
password = conn_parameters["password"]
host = conn_parameters["host"]
username = conn_parameters["user"]

database_url = settings.get_database_url()
engine = create_engine(database_url)

coastlines = config["NA_coastlines"]
# coastlines = "/home/nikola/coastline_45_14_kvarner.shp"
# coastlines = kvarner_coastline
regex_match = Configuration.ESA_dem_files_regex_match
data_folder = config["data_folder"]

sys.path.append("/usr/lib/python3/dist-packages")
from osgeo import gdal, ogr

# warnings.filterwarnings("ignore", message="Normalized/laundered field name:", category=UserWarning)
# warnings.filterwarnings("ignore", message="Normalized/laundered field name:", category=UserWarning, module="geopandas")
# warnings.simplefilter(action="ignore", category=UserWarning)

def sea_level_with_database():
    # Establish a connection to the database
    with psycopg2.connect(**conn_parameters) as connection:
        with connection.cursor() as cursor:
            # Loop to create united_geometries tables from 10 to 1
            for sea_level in tqdm(range(1, 10, 1)):
                output_file = os.path.join(
                    config["sea_level_rise_dir"], f"sea_level_output_{sea_level}.gpkg"
                )
                ogr_table_args = f"ogr2ogr -f GPKG {output_file} PG:'dbname=mnt tables=sea_level_rise_{sea_level}'"
                process = subprocess.Popen(
                    stdin=[ogr_table_args], stdout=subprocess.PIPE, shell=True
                ).communicate()
                print(process)

    # Close the database connection
    cursor.close()
    connection.close()


def create_contours(input_dem: str, output_contour: str, option: str, interval=1):
    assert option in ("-p", "")

    # create contour polygons
    if option == "-p":
        gdal_contour_cmd = f'gdal_contour -amin "AMIN" -amax "AMAX" -b 1 -i 1.0 {option} -f "ESRI Shapefile" "{input_dem}" "{output_contour}"'
    # Create contour lines
    else:
        gdal_contour_cmd = f'gdal_contour -a height -b 1 -i {interval} {option} -f "ESRI Shapefile" "{input_dem}" "{output_contour}"'

    os.system(gdal_contour_cmd)
    return output_contour


def explode_geometries(input_file: str):
    """
    This function uses ogr2ogr utility to explode multygeometries into distinct features.

    Parameters:
        input_file (str): vector file which needs exploding
    Note:
        ogr2ogr takes two arguments - destination dataset and source dataset
    """
    exploded_file = input_file.replace(".shp", "_exploded.shp")
    explode_cmd = (
        f'ogr2ogr -f "ESRI Shapefile" {exploded_file} {input_file} -explodecollections'
    )
    os.system(explode_cmd)
    return exploded_file


def create_filtered_contours(input_dem: str, tmp_contour, elevation_range, option: str):
    """
    use gdaldem to create contours and later filter only the wanted levels
    """
    assert option in ("-p", "")

    create_contours(input_dem, tmp_contour, option)
    tmp_contour = explode_geometries(tmp_contour)

    # Open the contour shapefile
    contour_ds = ogr.Open(tmp_contour)
    contour_layer = contour_ds.GetLayer()

    # Create a new shapefile for filtered contours
    driver = ogr.GetDriverByName("ESRI Shapefile")
    filtered_contour_file = tmp_contour.replace(".shp", "_filtered.shp")
    filtered_contour_ds = driver.CreateDataSource(filtered_contour_file)

    if option == "":
        filtered_contour_layer = filtered_contour_ds.CreateLayer(
            "filtered_contours", geom_type=ogr.wkbLineString
        )
    elif option == "-p":
        filtered_contour_layer = filtered_contour_ds.CreateLayer(
            "filtered_contours", geom_type=ogr.wkbPolygon
        )

    # Create field to store elevation information
    field_defn = ogr.FieldDefn("height", ogr.OFTReal)
    filtered_contour_layer.CreateField(field_defn)

    # Iterate through contour lines and filter based on elevation range
    for feature in contour_layer:
        if option == "":
            elevation = feature.GetField("height")
        elif option == "-p":
            elevation = feature.GetField("AMIN")

        if elevation_range[0] <= elevation <= elevation_range[1]:
            geometry = feature.GetGeometryRef()

            # Create filtered contour feature
            filtered_contour_feature = ogr.Feature(
                filtered_contour_layer.GetLayerDefn()
            )
            filtered_contour_feature.SetGeometry(geometry)
            filtered_contour_feature.SetField("height", elevation)

            # Add the feature to the layer
            filtered_contour_layer.CreateFeature(filtered_contour_feature)

    # Close datasets
    contour_ds = None
    filtered_contour_ds = None

    return filtered_contour_file


def contour_sea_level(lat: int, lon: int, res: int, elevation: int):
    """
    Main function for calculating potential flooding on DEM tiles by contour principle.
    It takes lat, lon, resolution, and elevation arguments for constructing DEM tile name.

    Parameters:
    - lat (int): Latitude of the DEM tile.
    - lon (int): Longitude of the DEM tile.
    - res (int): Resolution of the DEM tile. Accepted values: 10, 30, or 90 meters.
    - elevation (int): Elevation for filtering contours.

    Returns:
    - str: Path to the generated shapefile containing flooded areas.

    TODO: flood selection in a loop

    Algorithm goes in steps:
    1. Tile for which contours are calculated is expanded with neighboring tiles.
    2. Pixel selection - advisable to select 2 to 10 times higher altitude than wanted.
    3. VRT creation.
    4. Contour creation with gdal_contour and filtering up to elevation_range.
    5. Flood selection in a loop - check for intersection with coastline and neighboring flooded tiles.

    """
    assert res in (10, 30, 90)

    start_time = time.time()
    # Input DEM file
    if res == 90:
        input_dem_file = f"Copernicus_DSM_30_N{lat}_00_E0{lon}_00_DEM.tif"
        input_dem_file = (
            f"/mnt/volume-nbg1-1/shared/nikola/ESA_global_dem_90m/{input_dem_file}"
        )
    if res == 10:
        input_dem_file = f"Copernicus_DSM_03_N{lat}_00_E0{lon}_00_DEM.tif"
        input_dem_file = f"/mnt/volume-nbg1-1/satellite/eu_dem/dem10m/{input_dem_file}"

    # Elevation ranges for each meter from 1 to 10
    elevation_range = (0, elevation)

    with tempfile.TemporaryDirectory() as tmpdirname:
        # Make wider tile so that contours are correct
        input_dem_file = expand_tile(input_dem_file, tmpdirname)

        # Select pixels so that contours are made faster (always select at least 30 more meters than wanted)
        selected_sea_level_file = (
            f"{basename_withoutext(input_dem_file)}_selected_level.tif"
        )
        selected_sea_level_file = os.path.join(tmpdirname, selected_sea_level_file)
        gdal_select_pixels(
            input_dem_file, selected_sea_level_file, sea_level=35, replace_pixels=False
        )

        # Output contour shapefile
        contours_file = f"contours_{lat}_{lon}.shp"
        tmp_contour = os.path.join(tmpdirname, contours_file)

        # Create contours
        contours_file = create_filtered_contours(
            selected_sea_level_file, tmp_contour, elevation_range, option="-p"
        )

        # Cut contours to tile extent
        clip_vector_dataset(
            contours_file, contours_file, lon, lat, lon + 1, lat + 1
        )

        # imitate flooding with loops
        contours_file = f"/home/nikola/contours_45_14.shp"
        dst_shapefile_path = f"{tmpdirname}/dst_contours_45_14.shp"
        coastline_flooding_shp = f"coastal_flooding_{lat}_{lon}.shp"
        # flooding_selection(contours_file, coastlines, coastline_flooding_shp)
        # the function that utilizes gdal and queryes: update_flood_attribute_iterative
        # update_flood_attribute_iterative(dst_shapefile_path, contours_file, coastlines)

    # select_contours_by_height_and_intersection(contours_file,  coastlines, dst_shapefile_path)

    end_time = round((time.time() - start_time), 2)
    print(f"File done {dst_shapefile_path} in time {end_time} seconds")
    return dst_shapefile_path

# geopandas function, deprecated
def gpd_select_contours_by_height_and_intersection(contours_path, coastlines_path, output_path):
    subprocess.run(['ogr2ogr', '-s_srs', 'EPSG:4326', '-t_srs', 'EPSG:4326','-f', 'ESRI Shapefile', output_path, contours_path])

    coastlines = gpd.read_file(coastlines_path)
    contours = gpd.read_file(contours_path)
    output_data = contours.copy()

    for level in range(0, 11):
        # Select non-flooded contours and for current sea level height
        selected_contours_none = contours[(contours['height'] <= level) & (contours['flood'].isnull())]

        # Perform spatial join to get only the selected contours that intersect with coastlines
        if level == 0:

            selected_contours_intersects = gpd.sjoin(selected_contours_none, coastlines, predicate='intersects', how='left')
            selected_contours_intersects["flood"] = level
            selected_contours_intersects = selected_contours_intersects[["flood", "height", "geometry"]]

        # if level is higher than 0 then check for flooded contours, recursion should be added
        else:
            contours_with_flood = output_data[~output_data['flood'].isnull()]

            selected_contours_intersects = gpd.sjoin(selected_contours_none, contours_with_flood, predicate='intersects', how='inner')
            selected_contours_intersects["flood"] = level
            selected_contours_intersects = selected_contours_intersects[["flood", "height", "geometry"]]

        # Reset index before updating to handle potential duplicate index values
        output_data = output_data.reset_index(drop=True)
        selected_contours_intersects = selected_contours_intersects.reset_index(drop=True)

        # Save the result to a new shapefile
        output_data.update(selected_contours_intersects)

    output_data.to_file(output_path, driver='ESRI Shapefile')
    return


def update_flood_attribute_iterative_recursive(dst_shapefile_path, src_shapefile_path, coastline_shapefile_path):
    """
    Update flood attributes iteratively and recursively.

    Parameters:
    - dst_shapefile_path (str): The path to the destination shapefile.
    - src_shapefile_path (str): The path to the source shapefile that contains contours
    - coastline_shapefile_path (str): The path to the coastline shapefile.

    Returns:
    None
    """

    dst_shapefile_path_main = dst_shapefile_path
    coastlines_table = basename_withoutext(coastline_shapefile_path)
    src_table = basename_withoutext(src_shapefile_path)

    # TODO: add spatial index and srs
    # create_spatial_index_cmd = f'ogrinfo -sql "CREATE SPATIAL INDEX ON {src_table}" {src_shapefile_path}'
    # os.system(create_spatial_index_cmd)
    # create_spatial_index_cmd = f'ogrinfo -sql "CREATE SPATIAL INDEX ON {coastlines_table}" {coastline_shapefile_path}'
    # os.system(create_spatial_index_cmd)


    # Create a temporary GeoPackage file
    # tmp_gpkg_cmd = f"ogrmerge.py -f GPKG -overwrite_ds -a_srs 'EPSG:4326' -o {tmp_gpkg_path} {src_shapefile_path} {coastline_shapefile_path}"
    tmp_gpkg_path = "tmp_merged.gpkg"
    if os.path.exists(tmp_gpkg_path):
        os.remove(tmp_gpkg_path)

    # Merge source and coastline shapefiles into the temporary GeoPackage
    tmp_gpkg_cmd = f'ogr2ogr -f GPKG -a_srs EPSG:4326 -overwrite {tmp_gpkg_path} {src_shapefile_path} -nln {src_table} -nlt POLYGON '
    os.system(tmp_gpkg_cmd)
    tmp_gpkg_cmd = f'ogr2ogr -f GPKG -a_srs EPSG:4326 -update {tmp_gpkg_path} {coastline_shapefile_path} -nln {coastlines_table} -nlt LINESTRING'
    os.system(tmp_gpkg_cmd)


    for level in tqdm(range(0, 11)):
        
        if level == 0:
            dst_shapefile_path = dst_shapefile_path.replace(".shp", f"_{level}.shp")

            # Select features where height is below the current level and intersect with coastlines
            ogr2ogr_intersect = (
                f'ogr2ogr -f "ESRI Shapefile" -a_srs "EPSG:4326" {dst_shapefile_path} {tmp_gpkg_path} -dialect SQLite -sql '
                f'"SELECT DISTINCT contours.geom, contours.flood '
                f'FROM \'{src_table}\' as contours '
                f'INNER JOIN \'{coastlines_table}\' as coastlines '
                # f'ON (ST_intersects(contours.geom, coastlines.geom) OR ST_touches(contours.geom, coastlines.geom) OR ST_within(contours.geom, coastlines.geom)) '
                f'ON ST_relate(contours.geom, coastlines.geom) '
                f'WHERE contours.height <= {level}"'
            )
            os.system(ogr2ogr_intersect)

            # Set flooding level
            dst_shapefile_table = basename_withoutext(dst_shapefile_path)
            ogr_update = (
            f'ogrinfo {dst_shapefile_path} -dialect SQLite -sql '
            f'"UPDATE \'{dst_shapefile_table}\' '
            f'SET flood = {level}" '
            )
            os.system(ogr_update)

            # Append flooding shapefile - update GeoPackage
            tmp_gpkg_cmd = f"ogr2ogr -f GPKG -overwrite {tmp_gpkg_path} {dst_shapefile_path} -nln {dst_shapefile_table} -nlt POLYGON"
            os.system(tmp_gpkg_cmd)

            # Update flooding attributes in source table
            ogr_update = (
            f'ogrinfo {tmp_gpkg_path} -dialect SQLite -sql '
            f'"UPDATE \'{src_table}\' as src_table '
            f'SET flood = {level} '
            f'FROM \'{dst_shapefile_table}\' as flood_table '
            f'WHERE ST_equals(src_table.geom, flood_table.geom) "'
            )
            os.system(ogr_update)

        else:
            # Recursive part
            for recursive_level in range(level, -1, -1):

                dst_shapefile_path = dst_shapefile_path_main.replace(f".shp", f"_{level}_{recursive_level}.shp")

                # Export intersecting features to shapefile
                # intersecting feature - last shapefile that was created (level before)
                ogr2ogr_intersect = (
                    f'ogr2ogr -f "ESRI Shapefile" -a_srs "EPSG:4326" {dst_shapefile_path} {tmp_gpkg_path} -dialect SQLite -sql '
                    f'"SELECT DISTINCT contours.* '
                    f'FROM \'{src_table}\' as contours '
                    f'INNER JOIN (SELECT * FROM \'{src_table}\' WHERE flood IS NOT NULL) as flood '
                    f'ON ST_Intersects(contours.geom, flood.geom) OR ST_touches(contours.geom, flood.geom) OR ST_within(contours.geom, flood.geom) '
                    f'WHERE contours.height <= {level} AND contours.flood is NULL"'
                )
                os.system(ogr2ogr_intersect)

                # update flooding level shp
                # dst_shapefile_table = basename_withoutext(dst_shapefile_path)
                # ogr_update = (
                # f'ogrinfo {dst_shapefile_path} -dialect SQLite -sql '
                # f'"UPDATE \'{dst_shapefile_table}\' '
                # f'SET flood = {level}" '
                # )
                # os.system(ogr_update)

                # Append flooding shapefile - update GeoPackage
                tmp_gpkg_cmd = f"ogr2ogr -f GPKG -overwrite {tmp_gpkg_path} {dst_shapefile_path} -nln {dst_shapefile_table} -nlt POLYGON"
                os.system(tmp_gpkg_cmd)

                # Update flooding attributes in source table by flooding shapefile
                ogr_update = (
                f'ogrinfo {tmp_gpkg_path} -dialect SQLite -sql '
                f'"UPDATE \'{src_table}\' as src_table '
                f'SET flood = {level} '
                f'FROM \'{dst_shapefile_table}\' as flood_table '
                f'WHERE ST_equals(src_table.geom, flood_table.geom) "'
                )
                os.system(ogr_update)
            
        # Get out the list of feature IDs, just informative
        feature_ids = get_feature_ids_from_shapefile(dst_shapefile_path)
        print(feature_ids)
        
    return


def update_flood_attribute_iterative_v2_recursive_query(dst_shapefile_path, src_shapefile_path, coastline_shapefile_path):
    """
    Update flood attributes iteratively using a recursive SQL query.

    Parameters:
    - dst_shapefile_path (str): The path to the destination shapefile.
    - src_shapefile_path (str): The path to the source shapefile.
    - coastline_shapefile_path (str): The path to the coastline shapefile.

    Returns:
    None
    """

    dst_shapefile_path_main = dst_shapefile_path
    coastlines_table = basename_withoutext(coastline_shapefile_path)
    src_table = basename_withoutext(src_shapefile_path)

    # TODO: add spatial index and srs
    # create_spatial_index_cmd = f'ogrinfo -sql "CREATE SPATIAL INDEX ON {src_table}" {src_shapefile_path}'
    # os.system(create_spatial_index_cmd)
    # create_spatial_index_cmd = f'ogrinfo -sql "CREATE SPATIAL INDEX ON {coastlines_table}" {coastline_shapefile_path}'
    # os.system(create_spatial_index_cmd)


    # create temporary gpkg file
    # tmp_gpkg_cmd = f"ogrmerge.py -f GPKG -overwrite_ds -a_srs 'EPSG:4326' -o {tmp_gpkg_path} {src_shapefile_path} {coastline_shapefile_path}"
    tmp_gpkg_path = "tmp_merged.gpkg"
    if os.path.exists(tmp_gpkg_path):
        os.remove(tmp_gpkg_path)

    # hardcoded
    tmp_gpkg_cmd = f'ogr2ogr -f GPKG -a_srs EPSG:4326 -overwrite {tmp_gpkg_path} {src_shapefile_path} -nln {src_table} -nlt POLYGON '
    os.system(tmp_gpkg_cmd)
    tmp_gpkg_cmd = f'ogr2ogr -f GPKG -a_srs EPSG:4326 -update {tmp_gpkg_path} {coastline_shapefile_path} -nln {coastlines_table} -nlt LINESTRING'
    os.system(tmp_gpkg_cmd)


    for level in tqdm(range(0, 11)):
        
        if level == 0:
            dst_shapefile_path = dst_shapefile_path.replace(".shp", f"_{level}.shp")

            ogr2ogr_intersect = (
                f'ogr2ogr -f "ESRI Shapefile" -a_srs "EPSG:4326" {dst_shapefile_path} {tmp_gpkg_path} -dialect SQLite -sql '
                f'"SELECT DISTINCT contours.geom, contours.flood '
                f'FROM \'{src_table}\' as contours '
                f'INNER JOIN \'{coastlines_table}\' as coastlines '
                # f'ON (ST_intersects(contours.geom, coastlines.geom) OR ST_touches(contours.geom, coastlines.geom) OR ST_within(contours.geom, coastlines.geom)) '
                f'ON ST_relate(contours.geom, coastlines.geom) '
                f'WHERE contours.height <= {level}"'
            )
            os.system(ogr2ogr_intersect)

             # set flooding level
            dst_shapefile_table = basename_withoutext(dst_shapefile_path)
            ogr_update = (
            f'ogrinfo {dst_shapefile_path} -dialect SQLite -sql '
            f'"UPDATE \'{dst_shapefile_table}\' '
            f'SET flood = {level}" '
            )
            os.system(ogr_update)

            # update gpkg file - append flooding shapefile
            tmp_gpkg_cmd = f"ogr2ogr -f GPKG -overwrite {tmp_gpkg_path} {dst_shapefile_path} -nln {dst_shapefile_table} -nlt POLYGON"
            os.system(tmp_gpkg_cmd)

            # update flooding atributes in src table
            ogr_update = (
            f'ogrinfo {tmp_gpkg_path} -dialect SQLite -sql '
            f'"UPDATE \'{src_table}\' as src_table '
            f'SET flood = {level} '
            f'FROM \'{dst_shapefile_table}\' as flood_table '
            f'WHERE ST_equals(src_table.geom, flood_table.geom) "'
            )
            os.system(ogr_update)

        else:
            # Recursive part
            dst_shapefile_path = dst_shapefile_path_main.replace(".shp", f"_{level}_recursive.shp")

            # Export intersecting features with flood to shapefile using recursive query
            ogr2ogr_intersect = (
                f'ogr2ogr -f "ESRI Shapefile" -a_srs "EPSG:4326" {dst_shapefile_path} {tmp_gpkg_path} -dialect SQLite -sql '
                f'"WITH RECURSIVE connected_features AS ('
                f'  SELECT * FROM {src_table} WHERE flood = {level-1}'
                f'  UNION ALL '
                f'  SELECT t.* FROM {src_table} t '
                f'  INNER JOIN connected_features cf ON ST_Relate(t.geom, cf.geom) '
                f'      AND t.height < {level} AND t.flood IS NULL'
                f') '
                f'SELECT DISTINCT t.* FROM {src_table} t '
                f'LEFT JOIN connected_features cf ON ST_Intersects(t.geom, cf.geom) OR ST_touches(t.geom, cf.geom) OR ST_within(t.geom, cf.geom) '
                f'WHERE t.height <= {level} AND t.flood IS NULL AND cf.flood IS NOT NULL"'
            )
            
            os.system(ogr2ogr_intersect)

            # update flooding level shp
            dst_shapefile_table = basename_withoutext(dst_shapefile_path)
            ogr_update = (
            f'ogrinfo {dst_shapefile_path} -dialect SQLite -sql '
            f'"UPDATE \'{dst_shapefile_table}\' '
            f'SET flood = {level}" '
            )
            os.system(ogr_update)

            # append flooding shapefile - update gpkg file
            tmp_gpkg_cmd = f"ogr2ogr -f GPKG -overwrite {tmp_gpkg_path} {dst_shapefile_path} -nln {dst_shapefile_table} -nlt POLYGON"
            os.system(tmp_gpkg_cmd)

            # TODO: update without ogrinfo
            # update flooding atributes in src table by flooding shapefile

            ogr_update = (
            f'ogrinfo {tmp_gpkg_path} -dialect SQLite -sql '
            f'"UPDATE \'{src_table}\' as src_table '
            f'SET flood = {level} '
            f'FROM \'{dst_shapefile_table}\' as flood_table '
            f'WHERE ST_equals(src_table.geom, flood_table.geom) "'
            )
            os.system(ogr_update)


def get_feature_ids(shapefile_path):
    command = f'ogrinfo -al -q {shapefile_path}'
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, error = process.communicate()

    if process.returncode == 0:
        lines = output.decode('utf-8').split('\n')
        feature_ids = [line.split()[-1] for line in lines if 'Feature ID:' in line]
        return feature_ids
    else:
        print(f"Error executing ogrinfo: {error.decode('utf-8')}")
        return []

def get_feature_ids_from_shapefile(shapefile_path):
    """
    Get a list of feature IDs from a shapefile.

    Parameters:
    - shapefile_path (str): Path to the shapefile.

    Returns:
    - list: List of feature IDs.
    """
    # Open the shapefile
    dataSource = ogr.Open(shapefile_path)

    if dataSource is None:
        print("Error opening the shapefile.")
        return None

    # Get the first layer (assuming only one layer is present)
    layer = dataSource.GetLayer()

    # Get feature IDs and return the list
    feature_ids_list = [feature.GetFID() for feature in layer]

    # Close the shapefile
    dataSource = None
    return feature_ids_list


def clip_vector_dataset(
    input_vector_path, output_vector_path, min_x, min_y, max_x, max_y
):
    """
    Clip a vector dataset by a given extent.

    Args:
        input_vector_path (str): Path to the input vector dataset.
        output_vector_path (str): Path to the output clipped vector dataset.
        min_x (float): Minimum X-coordinate of the clipping extent.
        min_y (float): Minimum Y-coordinate of the clipping extent.
        max_x (float): Maximum X-coordinate of the clipping extent.
        max_y (float): Maximum Y-coordinate of the clipping extent.

    Returns:
        None
    """
    # Open the input vector dataset
    input_ds = ogr.Open(input_vector_path)
    if input_ds is None:
        raise Exception(f"Failed to open input vector dataset at {input_vector_path}")

    # Create the specified extent
    extent_ring = ogr.Geometry(ogr.wkbLinearRing)
    extent_ring.AddPoint(min_x, min_y)
    extent_ring.AddPoint(max_x, min_y)
    extent_ring.AddPoint(max_x, max_y)
    extent_ring.AddPoint(min_x, max_y)
    extent_ring.AddPoint(min_x, min_y)

    extent_polygon = ogr.Geometry(ogr.wkbPolygon)
    extent_polygon.AddGeometry(extent_ring)

    # Open the output vector dataset and create a new layer
    output_driver = ogr.GetDriverByName("ESRI Shapefile")
    if output_driver is None:
        raise Exception("ESRI Shapefile driver not available")

    output_ds = output_driver.CreateDataSource(output_vector_path)
    if output_ds is None:
        raise Exception(f"Exception: Failed to create output vector dataset at {output_vector_path}")

    # get spatial ref, and set output layer
    input_layer = input_ds.GetLayer()
    input_srs = input_layer.GetSpatialRef()
    output_layer = output_ds.CreateLayer("clipped", geom_type=ogr.wkbUnknown, srs=input_srs)

    if "contours" in input_vector_path:
        # Add "height" field to the output layer
        height_field = ogr.FieldDefn("height", ogr.OFTReal)
        output_layer.CreateField(height_field)

        # Add "flood" field to the output layer with integer type
        flood_field = ogr.FieldDefn("flood", ogr.OFTInteger)
        output_layer.CreateField(flood_field)
    
    else:
        # Add "flood" field to the output layer with integer type
        flood_field = ogr.FieldDefn("level", ogr.OFTInteger)
        output_layer.CreateField(flood_field)


    for feature in input_layer:
        input_geom = feature.GetGeometryRef()

        # Use Intersection method to clip the geometry
        clipped_geom = input_geom.Intersection(extent_polygon)

        # If the intersection is not empty, create a new feature
        if not clipped_geom.IsEmpty():
            new_feature = ogr.Feature(output_layer.GetLayerDefn())
            new_feature.SetGeometry(clipped_geom)

            # Retrieve and set the "height" or "level" attribute value
            if "contours" in input_vector_path:
                height_value = feature.GetField("height")
                new_feature.SetField("height", height_value)
                new_feature.SetField("flood", None)
            else:
                level_value = feature.GetField("level")
                new_feature.SetField("level", level_value)

            output_layer.CreateFeature(new_feature)

    # Clean up and close datasets
    input_ds = None
    output_ds = None


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


def extract_coordinates_from_tile_name(
    file_name: str, regex_match=r"_(S|N)(\d+)_00_(W|E)(\d+)"
):
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


def flooding_selection(
    polygon_shapefile_path, line_shapefile_path, output_shapefile_path
):
    """
    Function iterates over contours layer and selects interesects features - imitation of coastal flooding.

    Parameters:
        polygon_shapefile_path (str): path to contours file
        line_shapefile_path (str): path to line/coastline file
        output_shapefile_path (str): path to output file
    """
    # Read the shapefiles into GeoDataFrames
    polygons = gpd.read_file(polygon_shapefile_path)
    lines = gpd.read_file(line_shapefile_path)

    # Check for intersection
    for level in tqdm(range(0, 11, 1)):
        if level == 0:
            intersecting_features = gpd.overlay(polygons, lines, how="intersection")
        else:
            intersecting_features = gpd.overlay(
                polygons, reference_features, how="intersection"
            )

        # Select features where "AMIN" attribute is equal to current level
        selected_features = intersecting_features[
            intersecting_features["AMIN"] == {level}
        ]

        # concatenating the selected features to file
        if level > 0:
            selected_features = gpd.GeoDataFrame(
                pd.concat([reference_features, selected_features], ignore_index=True)
            )
        else:
            pass

        reference_features = selected_features

        # Save the selected features to a new shapefile
        selected_features.to_file(output_shapefile_path, driver="ESRI Shapefile")

    return selected_features


def gdal_create_buffer(input_file, output_buffer, buffer_dist):
    inputds = ogr.Open(input_file)
    input_layer = inputds.GetLayer()

    shpdriver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(output_buffer):
        shpdriver.DeleteDataSource(output_buffer)
    output_buffer_ds = shpdriver.CreateDataSource(output_buffer)
    buffer_layer = output_buffer_ds.CreateLayer(output_buffer, geom_type=ogr.wkbPolygon)
    featureDefn = buffer_layer.GetLayerDefn()

    for feature in input_layer:
        ingeom = feature.GetGeometryRef()
        geomBuffer = ingeom.Buffer(buffer_dist)

        outFeature = ogr.Feature(featureDefn)
        outFeature.SetGeometry(geomBuffer)
        buffer_layer.CreateFeature(outFeature)
        outFeature = None


def gdal_select_pixels(input_dem_file, output_extracted:str, sea_level=10, replace_pixels=True):
    """
    this function calls gdal calculator utility (gdal_calc.py) to extract wanted pixels.
    """
    if replace_pixels == True:
        calc = f"((A>=-999)*(A<={sea_level}))*{1} + (A>{sea_level})*-9999"

    if replace_pixels == False:
        calc = f"((A>=-999)*(A<={sea_level}))*A + (A>{sea_level})*-9999"

    gdal_calc_cmd = f'gdal_calc.py -A {input_dem_file} --outfile={output_extracted} --calc="{calc}" --NoDataValue=-9999 --co "COMPRESS=DEFLATE" --co "NBITS=32" --overwrite'
    os.system(gdal_calc_cmd)
    return output_extracted


def polygonize_raster(input_raster, output_polygons):
    """
    Function utilizes gdal_polygonize.py to convert given raster into vector poylgons.
    It is advisable to use 8 connectedness.
    """
    gdal_calc_cmd = f'gdal_polygonize.py -8 {input_raster} -overwrite -f "ESRI Shapefile" {output_polygons} {output_polygons} "level"'
    os.system(gdal_calc_cmd)
    return output_polygons


def get_coastline_from_pg(output_shapefile, extent="North America"):
    """
    A function that uses pgsql2shp to export existing coastlines from world table for North America into shp.
    """
    # right now hardcoded for north america
    if extent == "North America":
        pgsql_cmd = f'pgsql2shp -f {output_shapefile} -h {host} -u {username} -P {password} -g geom {database_name} "SELECT * FROM {schema}.world_coastlines WHERE ST_Intersects(geom, ST_MakeEnvelope(-180, 5, -20, 90, 4326));"'
    elif extent == "kvarner":
        pgsql_cmd = f'pgsql2shp -f {output_shapefile} -h {host} -u {username} -P {password} -g geom {database_name} "SELECT * FROM {schema}.world_coastlines WHERE ST_Intersects(geom, ST_MakeEnvelope(14, 45, 15, 46, 4326));"'
    os.system(pgsql_cmd)
    return print("Coastlines gotten!")


def calculate_polygon_difference(smaller_polygon, bigger_polygon, output_difference):
    """
    Function that calculates geospatial difference between two layers - bigger and smaller polygon.
    Based on geopandas, works best with shp files.

    Parameters:

    """
    # Read the input polygon layers
    gdf_smaller = gpd.read_file(smaller_polygon)
    gdf_bigger = gpd.read_file(bigger_polygon)

    # Calculate the difference between the two polygon layers
    with warnings.catch_warnings():
        warnings.simplefilter(action="ignore", category=UserWarning)
        difference_gdf = gdf_bigger.difference(gdf_smaller, align=True)

        # Save the result to a new shapefile
        difference_gdf.to_file(output_difference)
    return output_difference


def gdal_polygon_difference(method_ds, input_shp, output_ds):
    """
    Calculate difference between bigger and smaller polygon shapefile.
    """
    command = f"ogr_layer_algebra.py SymDifference -input_ds {input_shp} -method_ds {method_ds} -output_ds {output_ds} -nlt POLYGON"
    os.system(command)


def gdal_extract_features_by_location(input_file, intersection_file, output_file):
    """
    A function that selects overlapping or intersecting features from given layer that refer to intersection/reference layer.

    Parameters:
        input_file(str): input file from which to select features
        intersection_file(str): reference file by which features are selected
        output_file(str): path and name of the file to be saved
    """
    # try ST_Relate
    cmd_test_intersection = f"""
    ogr2ogr -f "ESRI Shapefile" -overwrite -progress {output_file} {input_file} -dialect sqlite -sql "SELECT polygon.geometry FROM {basename_withoutext(input_file)} polygon, '{intersection_file}'.{basename_withoutext(intersection_file)} coastline WHERE ST_Relate(polygon.geometry, coastline.geometry)" 
    """
    os.system(cmd_test_intersection)
    print(f"Features selected in the {input_file}!")
    return output_file

def merge_shapefiles(shapefile1_path, shapefile2_path, output_path, identifier_field):
    """
    Merge two shapefiles, removing duplicate features based on a specified identifier field.

    Parameters:
    - shapefile1_path (str): Path to the first shapefile.
    - shapefile2_path (str): Path to the second shapefile.
    - output_path (str): Path to save the merged shapefile.
    - identifier_field (str): Field used to identify and remove duplicate features.

    """
    # Read shapefiles into GeoDataFrames
    gdf1 = gpd.read_file(shapefile1_path)
    gdf2 = gpd.read_file(shapefile2_path)

    # Merge GeoDataFrames
    merged_gdf = gpd.GeoDataFrame(pd.concat([gdf1, gdf2], ignore_index=True), crs=gdf1.crs)

    # Drop duplicate features based on the identifier field
    merged_gdf.drop_duplicates(subset=identifier_field, keep='first', inplace=True)

    # Save the merged GeoDataFrame to a new shapefile
    merged_gdf.to_file(output_path)


def gpd_extract_features_by_location(input_shapefile:str, reference_shapefile:str, output_shapefile:str, condition = None, mode = "intersects", join = 'inner'):
    """
    function that uses geopandas to select and save features in input file by location of reference shapefile.

    Parameters:
        input_shapefile(str): input file from which to select features
        reference_shapefile(str): reference file by which features are selected
        output_shapefile(str): path and name of the file to be saved
    """
    # Read the input and reference shapefiles
    input_gdf = gpd.read_file(input_shapefile)
    if "SELECT" in reference_shapefile:
        reference_gdf = gpd.read_postgis(reference_shapefile, engine, geom_col="geom")
        reference_gdf.to_file(input_shapefile.replace(".shp", "_queried_coastlines.shp"))
        # print(f"reference_gdf saved temporarily")
    else:
        reference_gdf = gpd.read_file(reference_shapefile)
        # subset condition is optional (should be)
        if isinstance(condition, int):
            subset_condition = reference_gdf['level'] <= condition
            reference_gdf = reference_gdf[subset_condition]

    try:
        joined_gdf = gpd.sjoin(input_gdf, reference_gdf, how=join, predicate=mode)
    except:
        input_gdf = input_gdf[['level_left', 'geometry']]
        input_gdf = input_gdf.rename(columns={'level_left': 'level'})
        try:
            reference_gdf = reference_gdf[['level_left', 'geometry']]
            reference_gdf = reference_gdf.rename(columns={'level_left': 'level'})
        except:
            pass

        joined_gdf = gpd.sjoin(input_gdf, reference_gdf, how=join, predicate=mode)

    joined_gdf = joined_gdf.drop_duplicates(["geometry"])

    with warnings.catch_warnings():
        warnings.simplefilter(action="ignore")
        # Save the selected features to a new shapefile
        try:
            joined_gdf.to_file(output_shapefile)
        except:
            duplicated_columns = joined_gdf.columns[joined_gdf.columns.duplicated()]
            joined_gdf = joined_gdf.drop(columns=duplicated_columns)
            joined_gdf.to_file(output_shapefile)


def extract_intersecting_features(target_shapefile, other_shapefiles):
    """
    Extracts intersecting features between a target shapefile and a list of other shapefiles.

    Parameters:
    - target_shapefile (str): Path to the target shapefile.
    - other_shapefiles (list): List of paths to other shapefiles for intersection testing.

    Returns:
    - gpd.GeoDataFrame: GeoDataFrame containing the intersected geometries.
    """
    # Read the target shapefile
    target_gdf = gpd.read_file(target_shapefile)

    # Create an empty GeoDataFrame to store the result
    result_gdf = gpd.GeoDataFrame()

    # Iterate over other shapefiles
    for other_shapefile in other_shapefiles:
        # Read the other shapefile
        other_gdf = gpd.read_file(other_shapefile)

        # Check if there is any intersection
        intersection = gpd.overlay(target_gdf, other_gdf, how="intersection")

        if not intersection.empty:
            # Append the intersected geometries to the result GeoDataFrame
            result_gdf = result_gdf.append(intersection, ignore_index=True)


def label_connected_pixels(raster_file, line_shapefile, output_raster, label):
    """
    Labels pixels connected to a line in a raster using rasterio.

    Parameters:
    - raster_file (str): Path to the raster file.
    - line_shapefile (str): Path to the line shapefile.
    - output_raster (str): Path to the output labeled raster.
    - label: Label value for the connected pixels.

    Returns:
    - None

    TODO: test this function to see if it is faster than intersection testing with vectors
    """
    # Read the line shapefile
    line_gdf = gpd.read_file(line_shapefile)

    # Read the raster file
    with rasterio.open(raster_file) as src:
        # Create a mask for the pixels connected to the line
        mask = geometry_mask(
            line_gdf.geometry, out_shape=src.shape, transform=src.transform, invert=True
        )
        raster_data = src.read(1, masked=True)

        # neighbour part
        # Create a 3x3 convolution kernel
        kernel = np.ones((3, 3), dtype=bool)
        # Apply the kernel to the mask - problem with multidimensional array - scipy
        neighbor_mask = scipy.ndimage.convolve(mask, kernel)
        # neighbor_mask = np.convolve(mask, kernel, mode='same')

        # Update metadata for the output raster
        profile = src.profile
        profile.update(dtype=rasterio.int32, count=1)

        result = np.where(neighbor_mask, raster_data, profile["nodata"])

    # save file
    with rasterio.open(output_raster, "w", **profile) as dst:
        dst.write(result, 1)


def extract_pixels_by_mask(input_raster: str, mask_input: str, output_raster: str, mode: str, level: int):
    """
    Extracts pixels from an input raster based on operations with a vector mask.

    Parameters:
    - input_raster (str): Path to the input raster.
    - mask_input (str): Path to the vector mask.
    - output_raster (str): Path to the output raster.
    - mode (str): Operation mode - 'clip', 'difference', or 'invert'.
    - level (int): Threshold level for the operation (used in 'difference' mode).

    Returns:
    - None
    """
    assert mode in ("clip", "difference", "invert")

    vector_layer = gpd.read_file(mask_input)

    with rasterio.open(input_raster) as src:
        # Read the raster data
        raster_data = src.read(1, masked=True)
        profile = src.profile
        profile.update(dtype=rasterio.int32, count=1)

        mask = geometry_mask(
            vector_layer.geometry,
            out_shape=src.shape,
            transform=src.transform,
            invert=False,
        )

        if mode == "clip":
            features_values = np.logical_and(raster_data, mask)

        elif mode == "difference":
            features_values = np.where(mask, 1 - raster_data, raster_data)

            if level > 0:
                profile.update(nodata=0)

        elif mode == "invert":
            features_values = np.logical_and(raster_data, mask)
            features_values = np.invert(features_values)
            if level > 0:
                features_values = np.where(features_values, level+1, 0)

    with rasterio.open(output_raster, "w", **profile) as dst:
        dst.write(features_values, 1)


def extract_pixels_by_mask_concurrent_futures(input_raster: str, mask_input: str, output_raster: str, mode: str, level: int):
    """
    Extracts pixels from an input raster based on operations with a vector mask.

    Parameters:
    - input_raster (str): Path to the input raster.
    - mask_input (str): Path to the vector mask.
    - output_raster (str): Path to the output raster.
    - mode (str): Operation mode - 'clip', 'difference', or 'invert'.
    - level (int): Threshold level for the operation (used in 'difference' mode).

    Returns:
    - None
    """
    assert mode in ("clip", "difference", "invert")
    vector_layer = gpd.read_file(mask_input)

    with rasterio.open(input_raster) as src:
        profile = src.profile
        profile.update(dtype=rasterio.int32, count=1)

        raster_data = src.read(1, masked=True)

        def process_chunk(chunk):
            mask = geometry_mask(
                vector_layer.geometry,
                out_shape=chunk.shape,
                transform=src.transform,
                invert=False,
            )

            if mode == "clip":
                features_values_chunk = np.logical_and(chunk, mask)

            elif mode == "difference":
                features_values_chunk = np.where(mask, 1 - chunk, chunk)

                if level > 0:
                    profile.update(nodata=0)

            elif mode == "invert":
                features_values_chunk = np.logical_and(chunk, mask)
                features_values_chunk = np.invert(features_values_chunk)

            return features_values_chunk

        # Define chunk size for parallel processing
        chunk_size = 512
        chunks = [raster_data[i:i + chunk_size, j:j + chunk_size] for i in range(0, raster_data.shape[0], chunk_size) for j in range(0, raster_data.shape[1], chunk_size)]

        with concurrent.futures.ThreadPoolExecutor() as executor:
            features_values_chunks = list(executor.map(process_chunk, chunks))

        # Find the maximum width of the chunks
        max_width = max(chunk.shape[1] for chunk in features_values_chunks)

        # Pad the smaller chunks to match the maximum width
        features_values_chunks_padded = [np.pad(chunk, ((0, 0), (0, max_width - chunk.shape[1])), mode='constant') for chunk in features_values_chunks]

        features_values = np.concatenate(features_values_chunks_padded, axis=0)

    with rasterio.open(output_raster, "w", **profile) as dst:
        dst.write(features_values.astype(np.int32), 1)


def process_chunk(chunk, vector_layer, mode, level, src_transform, profile):
    """
    Process a chunk of raster data based on operations with a vector mask.

    Parameters:
    - chunk (numpy.ndarray): A chunk of raster data.
    - vector_layer (geopandas.geodataframe.GeoDataFrame): GeoDataFrame representing the vector mask.
    - mode (str): Operation mode - 'clip', 'difference', or 'invert'.
    - level (int): Threshold level for the operation (used in 'difference' mode).
    - src_transform (rasterio.transform.Affine): Affine transform of the source raster.

    Returns:
    - numpy.ndarray: Processed chunk of raster data.
    """

    mask = geometry_mask(
        vector_layer.geometry,
        out_shape=chunk.shape,
        transform=src_transform,
        invert=False,
    )

    if mode == "clip":
        features_values_chunk = np.logical_and(chunk, mask)
    elif mode == "difference":
        features_values_chunk = np.where(mask, 1 - chunk, chunk)
        if level > 0:
            profile.update(nodata=0)
    elif mode == "invert":
        features_values_chunk = np.logical_and(chunk, mask)
        features_values_chunk = np.invert(features_values_chunk)

    return features_values_chunk

def extract_pixels_by_mask_chunks(input_raster: str, mask_input: str, output_raster: str, mode: str, level: int):
    """
    Extracts pixels from an input raster based on operations with a vector mask.
    Uses multiprocessing to divide dataset into chunks and utilize function `process_chunk`.

    Parameters:
    - input_raster (str): Path to the input raster.
    - mask_input (str): Path to the vector mask.
    - output_raster (str): Path to the output raster.
    - mode (str): Operation mode - 'clip', 'difference', or 'invert'.
    - level (int): Threshold level for the operation (used in 'difference' mode).

    Returns:
    - None
    """
    assert mode in ("clip", "difference", "invert")

    vector_layer = gpd.read_file(mask_input)

    with rasterio.open(input_raster) as src:
        profile = src.profile
        profile.update(dtype=rasterio.int32, count=1)

        raster_data = src.read(1, masked=True)
        src_transform = src.transform

        # chunk_size = 512
        chunk_size = 1024
        chunks = [
            raster_data[i:i + chunk_size, j:j + chunk_size]
            for i in range(0, raster_data.shape[0], chunk_size)
            for j in range(0, raster_data.shape[1], chunk_size)
        ]

        with multiprocessing.Pool() as pool:
            results = pool.starmap(
                process_chunk,
                [(chunk, vector_layer, mode, level, src_transform, profile) for chunk in chunks]
            )

        max_width = max(result.shape[1] for result in results)
        features_values_chunks_padded = [
            np.pad(result, ((0, 0), (0, max_width - result.shape[1])), mode='constant')
            for result in results
        ]

        features_values = np.concatenate(features_values_chunks_padded, axis=0)

    with rasterio.open(output_raster, "w", **profile) as dst:
        dst.write(features_values.astype(np.int32), 1)

def process_chunk_2(chunk, mask, mode, level):
    """
    Process a chunk of raster data based on operations with a vector mask.

    Parameters:
    - chunk (numpy.ndarray): A chunk of raster data.
    - mask (numpy.ndarray): Vector mask.
    - mode (str): Operation mode - 'clip', 'difference', or 'invert'.
    - level (int): Threshold level for the operation (used in 'difference' mode).

    Returns:
    - numpy.ndarray: Processed chunk of raster data.
    """

    if mode == "clip":
        features_values_chunk = np.logical_and(chunk, mask)
    elif mode == "difference":
        features_values_chunk = np.where(mask, 1 - chunk, chunk)
        if level > 0:
            features_values_chunk = np.where(mask, level, features_values_chunk)
    elif mode == "invert":
        features_values_chunk = np.logical_and(chunk, mask)
        features_values_chunk = np.invert(features_values_chunk)

    return features_values_chunk

def extract_pixels_by_mask_chunks_2(input_raster, mask_input, output_raster, mode, level):
    """
    Extracts pixels from an input raster based on operations with a vector mask.
    Uses multiprocessing to divide dataset into chunks and utilize function `process_chunk`.

    Parameters:
    - input_raster (str): Path to the input raster.
    - mask_input (str): Path to the vector mask.
    - output_raster (str): Path to the output raster.
    - mode (str): Operation mode - 'clip', 'difference', or 'invert'.
    - level (int): Threshold level for the operation (used in 'difference' mode).

    Returns:
    - None
    """
    assert mode in ("clip", "difference", "invert")

    with rasterio.open(input_raster) as src:
        profile = src.profile
        profile.update(dtype=rasterio.int32, count=1)

        raster_data = src.read(1, masked=True)
        src_transform = src.transform

        # Read the vector mask using gdal/ogr
        mask_ds = gdal.OpenEx(mask_input, gdal.OF_VECTOR)
        mask_layer = mask_ds.GetLayer()

        # Create a temporary in-memory raster for the mask
        mem_driver = gdal.GetDriverByName('MEM')
        mem_raster = mem_driver.Create('', raster_data.shape[1], raster_data.shape[0], 1, gdal.GDT_Byte)

        # Set the geotransform and projection of the temporary raster
        mem_raster.SetGeoTransform(src_transform.to_gdal())
        mem_raster.SetProjection(src.crs.to_wkt())

        # Rasterize the vector mask into the temporary raster
        gdal.RasterizeLayer(mem_raster, [1], mask_layer, burn_values=[1])

        # Read the mask from the temporary raster
        mask = mem_raster.GetRasterBand(1).ReadAsArray()

        mask = mask[:raster_data.shape[0], :raster_data.shape[1]]

        # chunk_size = 512
        chunk_size = 1201
        chunks = [
            raster_data[i:i + chunk_size, j:j + chunk_size]
            for i in range(0, raster_data.shape[0], chunk_size)
            for j in range(0, raster_data.shape[1], chunk_size)
        ]

        with multiprocessing.Pool() as pool:
            results = pool.starmap(
                process_chunk_2,
                [(chunk, mask, mode, level) for chunk in chunks]
            )

        max_width = max(result.shape[1] for result in results)
        features_values_chunks_padded = [
            np.pad(result, ((0, 0), (0, max_width - result.shape[1])), mode='constant')
            for result in results
        ]

        features_values = np.concatenate(features_values_chunks_padded, axis=0)

    with rasterio.open(output_raster, "w", **profile) as dst:
        dst.write(features_values.astype(np.int32), 1)


def coastal_flooding_rasters_sum(input_raster, raster_to_add, raster_sum, level: int):
    """
    Adds different layers of flooding to a GeoTIFF file.

    Parameters:
    - input_raster (str): Path to the input GeoTIFF file.
    - raster_to_add (str): Path to the GeoTIFF file containing additional flooding layers.
    - raster_sum (str): Path to the output GeoTIFF file where the sum of flooding layers will be stored.
    - level (int): An integer representing the level or weight to be assigned to the additional flooding layer.
                  If level > 0, each pixel in the additional layer will be assigned a value of level + 1.

    Returns:
    - None
    """
    with rasterio.open(input_raster) as src:
        input_raster_data = src.read(1)

    with rasterio.open(raster_to_add) as src:
        profile = src.profile
        add_raster_data = src.read(1)

    if level > 0:
        level = level + 1

    add_raster_data = np.where(add_raster_data == 1, level, profile["nodata"])
    sum_values = np.add(input_raster_data, add_raster_data)

    with rasterio.open(raster_sum, "w", **profile) as dst:
        dst.write(sum_values, 1)


def coastal_flooding_raster_subtraction(
    flooding_raster: str, output_raster: str, subtract=1
):
    """
    Subtracts a specified value from all pixels in a coastal flooding layer raster.

    Parameters:
    - flooding_raster (str): Path to the input GeoTIFF file containing the coastal flooding layer.
    - output_raster (str): Path to the output GeoTIFF file where the subtracted raster will be stored.
    - subtract (int, optional): The value to subtract from each pixel in the flooding layer.
                               Defaults to 1.

    Returns:
    - None
    """
    with rasterio.open(flooding_raster) as src:
        flood_data = src.read(1)
        profile = src.profile

    flood_data = np.where(flood_data, flood_data - subtract, profile["nodata"])

    with rasterio.open(output_raster, "w", **profile) as dst:
        dst.write(flood_data, 1)


def subtract_one_from_shapefile_attribute(shapefile_path, attribute_column_name, output_shapefile):
    """
    Subtracts 1 from a specified attribute/column in a shapefile.

    Parameters:
    - shapefile_path (str): Path to the input shapefile.
    - attribute_column_name (str): Name of the attribute/column to be subtracted by 1.
    - output_shapefile (str): Path to the output shapefile where the modified data will be stored.

    Returns:
    - None
    """
    # Read the shapefile
    gdf = gpd.read_file(shapefile_path)

    # Subtract 1 from the specified attribute column
    gdf[attribute_column_name] = gdf[attribute_column_name] - 1

    # Save the modified shapefile
    gdf.to_file(output_shapefile)


def basename_withoutext(input_file):
    """
    return basename of the file without extension
    """
    return os.path.splitext(os.path.basename(input_file))[0]


def give_tmp_file_name(input_file: str, level: int, name: str):
    """
    helper function for giving names to temporary files created in the loop function
    """
    return f"{basename_withoutext(input_file)}_{name}_{level}.shp"


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
    # if neighbours == 0: 
    #     return [central_file]
    
    if include_central == True:
        neighbour_files = [central_file]
    else:
        neighbour_files = []

    file_name = os.path.basename(central_file)
    tile_folder = os.path.dirname(central_file)
    # find file coords
    # match something like "N23_00_W123"
    geocell_id_match = re.search(r"_(S|N)(\d+)_00_(W|E)(\d+)", file_name)
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
                        tile_name = f"Copernicus_DSM_30_{lat_direction}{j:02d}_00_{k}{i:03d}_00_DEM_{coastal_flood_sufix}.shp"
                        tile_path = os.path.join(tile_folder, tile_name)
                        if os.path.exists(tile_path) and tile_path != central_file:
                            neighbour_files.append(tile_path)

                else:
                    tile_name = f"Copernicus_DSM_30_{lat_direction}{j:02d}_00_{lon_direction}{i:03d}_00_DEM_{coastal_flood_sufix}.shp"
                    tile_path = os.path.join(tile_folder, tile_name)
                    if os.path.exists(tile_path) and tile_path != central_file:
                        neighbour_files.append(tile_path)

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


def coastal_flooding_pixel_prediction(input_file: str, sea_level = 10, coast_distance = 0, skip_if_exists = False):
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
    data_folder = config["coastal_flooding_EU_90m"]

    output_file = f"{data_folder}/{basename_withoutext(original_file)}_coastal_flood.shp"

    if skip_if_exists == True:
        if os.path.exists(output_file):
            logging.info(f"File {output_file} existed, so it was skipped")
            return
        
    neighbouring_coastal_files = collect_neighbouring_coastal_flood_files(output_file, neighbours=coast_distance)
    
    lat_direction, lat, lon_direction, lon = extract_coordinates_from_tile_name(input_file)
    lat = int(lat)
    lon = int(lon)

    with tempfile.TemporaryDirectory() as tmpdirname:
        # build_expanded_tif(input_file, coastal_dataset, tmpdirname)

        # build neighbouring vector layer
        neighbours_str = " ".join(neighbouring_coastal_files)
        neighbours_shp = f"{tmpdirname}/{basename_withoutext(original_file)}_neighbours.shp"

        # merge neighbouring shapefiles, -single option is for merging layers into a single shp
        merge_neighbours_cmd = f"ogrmerge.py -o {neighbours_shp} {neighbours_str} -single -progress "
        os.system(merge_neighbours_cmd)

        # clip neighbours, only tiny bit of polygons is needed
        neighbours_shp_clipped = neighbours_shp.replace(".shp", "_clipped.shp")
        try:
            if lat_direction == "N" and lon_direction == "W": # Americas
                clip_vector_dataset(neighbours_shp, neighbours_shp_clipped, -lon+1.001, lat-0.001, -lon - 0.001, lat + 1.001)
            elif lat_direction == "N" and lon_direction == "E": # Europe
                clip_vector_dataset(neighbours_shp, neighbours_shp_clipped, lon+1.001, lat-0.001, lon - 0.001, lat + 1.001)

        except Exception:
            neighbours_shp_clipped = None

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
            polygonize_raster(selected_sea_level_file, polygonized_pixels_file)

            # spatial relation with sea or spatial relation with flood
            intersected_features_file = f"{tmpdirname}/{basename_withoutext(input_file)}_sea_level_intersected_{level}.shp"
            
            if level == 0 and coast_distance == 0:
                geocellid = geocellid_from_file_name(input_file)
                coastlines_query = f"""
                SELECT DISTINCT wc.geom as geom FROM osm.world_coastlines wc
                INNER JOIN osm.esa_global_dem_grid egdg
                ON ST_Intersects(wc.geom, egdg.geometry)
                WHERE egdg.geocellid = '{geocellid}'
                """
                gpd_extract_features_by_location(polygonized_pixels_file, coastlines_query, intersected_features_file)

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

def altitude_filter_files_list(dem_files: str, threshold_altitude=10, direction = 'descending'):
    """
    filter all dem tiles with altitude lower than specified treshold

    Parameters:
        dem_files(str): a list of dem files
        threshold_altitude(int): a number by which files are filtered
    """
    tiles_with_low_altitude = []

    for file_path in tqdm(dem_files):
        dataset = gdal.Open(file_path)

        if dataset is not None:
            # Get altitude values as a NumPy array
            altitude_array = dataset.ReadAsArray()

            # Check if any pixel has altitude less than the threshold
            if direction == "descending":
                if (altitude_array < threshold_altitude).any():
                    tiles_with_low_altitude.append(file_path)
            elif direction == "ascending":
                if (altitude_array > threshold_altitude).any():
                    tiles_with_low_altitude.append(file_path)

            # Close the dataset
            dataset = None

    return tiles_with_low_altitude



def has_empty_geometries(shapefile_path):
    """
    Check if a shapefile contains empty geometries.

    Parameters:
    - shapefile_path (str): The path to the shapefile.

    Returns:
    - bool: True if the shapefile contains empty geometries, False otherwise.

    """
    try:
        # Read the shapefile into a GeoDataFrame
        gdf = gpd.read_file(shapefile_path)

        # Check for empty geometries
        empty_geometries = gdf[gdf.is_empty]

        return gdf.is_empty[0]
    except Exception as e:
        return True
    

def create_file_name_from_geocellid(geocell_id, folder_path):
    regex_match = r'([NS])(\d+)([WE])(\d+)'
    geocell_id_match = re.match(regex_match, geocell_id)
    
    if geocell_id_match:
        lat_direction = geocell_id_match.group(1)
        lat_value = geocell_id_match.group(2)
        lon_direction = geocell_id_match.group(3)
        lon_value = geocell_id_match.group(4)
        if folder_path == config["esa_global_dem_90_dir"]:
            file_name = f"{folder_path}/Copernicus_DSM_30_{lat_direction}{lat_value}_00_{lon_direction}{lon_value}_00_DEM.tif"

        elif folder_path == config["esa_global_dem_30_dir"]:
            file_name = f"{folder_path}/Copernicus_DSM_10_{lat_direction}{lat_value}_00_{lon_direction}{lon_value}_00_DEM.tif"

        return file_name
    else:
        return None
    

def get_grid_and_coastline_gdf(continent = 'Europe'):
    """
    The intersection of coastline and esa_dem_grid file is queried and put into gdf object.
    """
    world_continents_boundaries = config["world_continents_boundaries_table"]
    
    # only coastline tiles
    grid_coastlines_query = f"""
    SELECT DISTINCT geocellid, grid.geometry
    FROM osm.esa_global_dem_grid as grid 
    JOIN osm.world_continents_boundaries wcb
    ON ST_Intersects(grid.geometry, wcb.geometry)
    JOIN osm.world_water_bodies_esri wwb
    ON ST_Intersects(grid.geometry, wwb.geom)
    WHERE wcb.continent = '{continent}' AND wwb.type = 'Ocean or Sea';
    """
    # coastline tiles and tiles close to those tiles (distances 0 and 1)
    grid_coastlines_query = f""" 
    SELECT DISTINCT  ON (geocellid) geocellid, geometry, distance
    FROM (
	SELECT DISTINCT geocellid, grid.geometry, 0 as distance
	FROM osm.esa_global_dem_grid as grid 
	JOIN osm.world_continents_boundaries wcb
	ON ST_Intersects(grid.geometry, wcb.geometry)
	JOIN osm.world_water_bodies_esri wwb
	ON ST_Intersects(grid.geometry, wwb.geom)
	WHERE wcb.continent = 'Europe' AND wwb.type = 'Ocean or Sea'
	UNION
	SELECT DISTINCT grid.geocellid, grid.geometry, 1 as distance
	FROM (
		SELECT DISTINCT grid.geometry
		FROM osm.esa_global_dem_grid as grid 
		JOIN osm.world_continents_boundaries wcb
		ON ST_Intersects(grid.geometry, wcb.geometry)
		JOIN osm.world_water_bodies_esri wwb
		ON ST_Intersects(grid.geometry, wwb.geom)
		WHERE wcb.continent = 'Europe' AND wwb.type = 'Ocean or Sea') as subquery,
	osm.esa_global_dem_grid as grid 
	WHERE ST_touches(grid.geometry, subquery.geometry)) as united
	;
    """

    coastline_grid_gdf = gpd.read_postgis(grid_coastlines_query, engine, geom_col = 'geometry')

    return coastline_grid_gdf


def sea_level_precheck(files_dir: str, continent = "Europe"):
    """
    this function iterates over grid of dem files and check if dem tiles are close to sea and if they contain pixels with latitude lower than 10m.

    Algorithm steps:
        1. Vector check: vector dem grid tiles proximity to coastline and geo extent (query + update proximity attribute)
        2. Raster check: altitude and extent - which raster dem tiles are lower than 10m and in given extent   
        3. Combine - create list by which the sea_level iterration will occur
    """
    low_altitude_files = config["low_altitude_dem_files"]

    overwrite = True
    if overwrite == True or os.path.exists(low_altitude_files) == False:
        files_list = geofilter_paths_list(files_dir, by = "continent", c_name=continent)
        low_dem_tiles = altitude_filter_files_list(files_list)
    #     file_paths_to_write = ", ".join(low_dem_tiles)

    #     with open(low_altitude_files, 'w') as file:
    #         file.write(file_paths_to_write) 

    # a part for processing (or skipping the processing) of low altitude files
    # elif os.path.exists(low_altitude_files):
    #     with open(low_altitude_files, 'r') as file:
    #         low_dem_tiles = file.read()
    #         low_dem_tiles = low_dem_tiles.split(', ')


    if continent == "North America":
        coastal_neighbours_df = gpd.read_file(config["na_coastal_flooding_neighbours"]) # tiles for NA
        # if tiles grid file was created by merging vector layers n stuff
        coastal_neighbours_df['distance'] = coastal_neighbours_df['layer'].str.extract(r'(\d+)', expand=False).astype(int)

    elif continent == "Europe":
        coastal_neighbours_df = get_grid_and_coastline_gdf(continent) #tiles for EU

    # now the sorting part
    folder_path = config["esa_global_dem_90_dir"]
    coastal_neighbours_df['file_path'] = coastal_neighbours_df['geocellid'].apply(create_file_name_from_geocellid, args=(folder_path,))

    # remove all the files that are not in low altutide filtered list, and sort by coastal distance
    coastal_neighbours_df = coastal_neighbours_df[coastal_neighbours_df['file_path'].isin(low_dem_tiles)]
    coastal_neighbours_df = coastal_neighbours_df.sort_values(by=["distance", "geocellid"])

    # dict is returned
    coastal_neighbours_dict = dict(zip(coastal_neighbours_df['file_path'], coastal_neighbours_df['distance']))

    return coastal_neighbours_dict
