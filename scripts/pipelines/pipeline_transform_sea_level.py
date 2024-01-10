import psycopg2
from psycopg2 import sql
import sys
import os
from tqdm import tqdm
import subprocess
import numpy as np
import time
from scipy.ndimage import label
import warnings
import rasterio
from rasterio.features import geometry_mask
import scipy
import geopandas as gpd
from skimage.morphology import binary_dilation
from skimage import measure
from skimage.io import imread, imsave
from skimage.color import label2rgb
import tempfile


# Define the database connection parameters
sys.path.append("/home/nikola/4_north_america/GeoDataPump/scripts")
import settings
config = settings.get_config()
conn_parameters = settings.get_conn_parameters()
schema = settings.get_schema()

database_name = conn_parameters["database"]
password = conn_parameters["password"]
host = conn_parameters["host"]
username = conn_parameters["user"]
coastlines = config["NA_coastlines"]

sys.path.append('/usr/lib/python3/dist-packages')
from osgeo import gdal, ogr


def sea_level_with_database():
    # Establish a connection to the database
    with psycopg2.connect(**conn_parameters) as connection:
        with connection.cursor() as cursor:

            # Loop to create united_geometries tables from 10 to 1
            for sea_level in tqdm(range(1, 10, 1)):
                output_file = os.path.join(config["sea_level_rise_dir"], f"sea_level_output_{sea_level}.gpkg")
                ogr_table_args = f"ogr2ogr -f GPKG {output_file} PG:'dbname=mnt tables=sea_level_rise_{sea_level}'"
                process = subprocess.Popen(stdin=[ogr_table_args], stdout=subprocess.PIPE, shell=True).communicate()
                print(process)

    # Close the database connection
    cursor.close()
    connection.close()

def sea_level_split():
    for sea_level in tqdm(range(1, 10, 1)):
        output_file = os.path.join(config["sea_level_rise_dir"], f"sea_level_output_{sea_level}.gpkg")
        ogr_table_args = f"ogr2ogr -f GPKG {output_file} PG:'dbname=mnt tables=sea_level_rise_{sea_level}'"
        process = subprocess.Popen(stdin=[ogr_table_args], stdout=subprocess.PIPE, shell=True).communicate()
        print(process)    
    # ogr2ogr -f GPKG sea_level_output_10.gpkg PG:"dbname='mnt' tables='sea_level_rise_10'"        
    # ogr2ogr -f GPKG valid_sea_level_output_10.gpkg sea_level_output_10.gpkg -dialect sqlite -sql "select ST_MakeValid(geom) as geom, * from sea_level_rise_10"
    # ogr2ogr -overwrite sea_level_difference_10.gpkg sea_level_output_10.gpkg -dialect SQLite -sql "select st_difference(a.geometry, b.geometry) as geometry from 'sea_level_rise_10' as a, 'sea_level_output_9.gpkg' as b"
    # ogr2ogr -overwrite sea_level_difference_10.gpkg sea_levels.gpkg -dialect SQLite -sql "SELECT ST_Difference(a.geom, b.geom) AS geometry FROM 'sea_level_output_10_sea_level_rise_10' AS a, 'sea_level_output_9_sea_level_rise_9' AS b" -progress
    # ogr_layer_algebra.py SymDifference -input_ds sea_level_output_10.gpkg -method_ds sea_level_output_9.gpkg -output_ds sea_level_difference_10.shp -nlt POLYGON
    # ogr_layer_algebra.py SymDifference -input_ds valid_sea_level_10.shp -method_ds sea_level_output_9.gpkg -output_ds sea_level_difference_10.shp -nlt POLYGON



def valid_geometries_shell():
    for i in tqdm(range(1, 9, 1)):
        input_file = f"/mnt/volume-nbg1-1/shared/nikola/sea_level_rise/sea_level_output_{i}.gpkg"
        output_file = f"/mnt/volume-nbg1-1/shared/nikola/sea_level_rise/valid_sea_level_{i}.shp"
        ogr_makevalid_cmd = f'ogr2ogr {output_file} {input_file} -makevalid -progress'
        os.system(ogr_makevalid_cmd)
        # create spatial index 
        # ogrinfo -sql "CREATE SPATIAL INDEX ON valid_sea_level_10" valid_sea_level_10.shp


def sea_level_gdal_raster():
    # Open the DEM raster
    dem_file = 'path/to/your/dem.tif'
    dem_dataset = gdal.Open(dem_file)

    # Get DEM array
    dem_array = dem_dataset.ReadAsArray()

    # Define the sea level range
    sea_level_range = (0, 10)

    # Create a binary mask for pixels within the sea level range
    sea_level_mask = np.logical_and(dem_array >= sea_level_range[0], dem_array <= sea_level_range[1])

    # Load the coastal zone raster or shapefile
    coastal_zone_file = 'path/to/your/coastal_zone.tif'
    coastal_zone_dataset = gdal.Open(coastal_zone_file)
    coastal_zone_array = coastal_zone_dataset.ReadAsArray()

    # Create a mask for pixels close to the sea
    coastal_mask = np.logical_and(sea_level_mask, coastal_zone_array == 1)

    # Save the result as a new raster
    output_file = 'path/to/your/output_mask.tif'

    # Create a new raster with the same geotransform and projection as the input DEM
    output_dataset = gdal.GetDriverByName('GTiff').Create(output_file, dem_dataset.RasterXSize, dem_dataset.RasterYSize, 1, gdal.GDT_Byte)

    # Set the geotransform and projection
    output_dataset.SetGeoTransform(dem_dataset.GetGeoTransform())
    output_dataset.SetProjection(dem_dataset.GetProjection())

    # Write the mask to the raster
    output_dataset.GetRasterBand(1).WriteArray(coastal_mask.astype(np.uint8))

    # Close the datasets
    dem_dataset = None
    output_dataset = None

def create_contours(input_dem, output_contour, interval):
    # Create contour lines
    gdal_contour_cmd = f'gdal_contour -a height -b 1 -i {interval} -p -f "ESRI Shapefile" "{input_dem}" "{output_contour}"'
    os.system(gdal_contour_cmd)
    return


def create_filtered_contours(input_dem, output_contour, elevation_range, option):
    assert option in ("-p", "")
    # Create contoures
    if option == "-p":
        gdal_contour_cmd = f'gdal_contour -amin AMIN -amax AMAX -b 1 -i 1.0 {option} -f "ESRI Shapefile" "{input_dem}" "{output_contour}"'
    else:
        gdal_contour_cmd = f'gdal_contour -a height -b 1 -i 1.0 {option} -f "ESRI Shapefile" "{input_dem}" "{output_contour}"'

    os.system(gdal_contour_cmd)

    # Open the contour shapefile
    contour_ds = ogr.Open(output_contour)
    contour_layer = contour_ds.GetLayer()

    # Create a new shapefile for filtered contours
    driver = ogr.GetDriverByName('ESRI Shapefile')
    filtered_contour_file = output_contour.replace('.shp', '_filtered.shp')
    filtered_contour_ds = driver.CreateDataSource(filtered_contour_file)
    filtered_contour_layer = filtered_contour_ds.CreateLayer('filtered_contours', geom_type=ogr.wkbLineString)

    # Create field to store elevation information
    field_defn = ogr.FieldDefn('height', ogr.OFTReal)
    filtered_contour_layer.CreateField(field_defn)

    # Iterate through contour lines and filter based on elevation range
    for feature in contour_layer:
        if option == "":
            elevation = feature.GetField('height')
        else:
            elevation = feature.GetField('AMIN')

        if elevation_range[0] <= elevation <= elevation_range[1]:
            geometry = feature.GetGeometryRef()

            # Create filtered contour feature
            filtered_contour_feature = ogr.Feature(filtered_contour_layer.GetLayerDefn())
            filtered_contour_feature.SetGeometry(geometry)
            filtered_contour_feature.SetField('height', elevation)

            # Add the feature to the layer
            filtered_contour_layer.CreateFeature(filtered_contour_feature)

    # Close datasets
    contour_ds = None
    filtered_contour_ds = None

    return filtered_contour_file


def contour_sea_level():
    # Input DEM file
    start_time = time.time()
    lat = 45
    lon = 14
    input_dem_file = '/mnt/volume-nbg1-1/shared/nikola/ESA_global_dem_90m/Copernicus_DSM_30_N45_00_E014_00_DEM.tif'

    # Output contour shapefile
    output_contour_file = '/home/nikola/4_north_america/GeoDataPump/data/contours_45_14.shp'

    # Elevation ranges for each meter from 1 to 10
    elevation_ranges = [(i, i+1) for i in range(1, 11)]
    elevation_range = (1, 10)
    # elevation_range = (10, 10)

    # Create contours
    # create_contours(input_dem_file, output_contour_file, 1)
    file = create_filtered_contours(input_dem_file, output_contour_file, elevation_range, option="")
    end_time = round((time.time()-start_time), 2)
    print(f"File done {file} in time {end_time} seconds")
    return file


def gdal_create_buffer(input_file, output_buffer, buffer_dist):
    inputds = ogr.Open(input_file)
    input_layer = inputds.GetLayer()

    shpdriver = ogr.GetDriverByName('ESRI Shapefile')
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


def gdal_select_pixels(input_dem_file, output_extracted, sea_level = 10):
    """
    this function calls gdal calculator utility (gdal_calc.py) to extract wanted pixels.
    """
    calc = f"((A>=-999)*(A<={sea_level}))*{1} + (A>{sea_level})*-9999"

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
        pgsql_cmd = f'pgsql2shp -f {output_shapefile} -h {host} -u {username} -P {password} -g geom {database_name} "SELECT * FROM {schema}.world_coastlines WHERE ST_Intersects(geom, ST_MakeEnvelope(-140, 0, -30, 80, 4326));"'
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
        warnings.simplefilter(action='ignore', category=UserWarning)
        difference_gdf = gdf_bigger.difference(gdf_smaller, align = True)

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


def gpd_extract_features_by_location(input_shapefile, reference_shapefile, output_shapefile):
    """
    function that uses geopandas to select and save features in input fileby location of reference shapefile.

    Parameters:
        input_shapefile(str): input file from which to select features
        reference_shapefile(str): reference file by which features are selected
        output_shapefile(str): path and name of the file to be saved
    """
    # Read the input and reference shapefiles
    input_gdf = gpd.read_file(input_shapefile)
    reference_gdf = gpd.read_file(reference_shapefile)

    # Do the spatial join to identify intersecting features
    joined_gdf = gpd.sjoin(input_gdf, reference_gdf, how='inner', predicate='intersects')

    with warnings.catch_warnings():
        warnings.simplefilter(action='ignore', category=UserWarning)
        # Save the selected features to a new shapefile
        joined_gdf.to_file(output_shapefile)


def label_connected_pixels(raster_file, line_shapefile, output_raster, label):
    """
    functon uses rasterio to label intersected and connected pixels to coastline
    TODO: test this function to see if it is faster than intersection testing with vectors
    """
    # Read the line shapefile
    line_gdf = gpd.read_file(line_shapefile)

    # Read the raster file
    with rasterio.open(raster_file) as src:
        # Create a mask for the pixels connected to the line
        mask = geometry_mask(line_gdf.geometry, out_shape=src.shape, transform=src.transform, invert=True)
        raster_data = src.read(1, masked=True)

        #neighbour part
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
    with rasterio.open(output_raster, 'w', **profile) as dst:
        dst.write(result, 1)


def extract_pixels_by_mask(input_raster:str, mask_input:str, output_raster:str, mode:str, level:int):
    """
    Function that takes input raster and does operation with vector mask
    """
    assert mode in ("clip", "difference", "invert")
    vector_layer = gpd.read_file(mask_input)

    with rasterio.open(input_raster) as src:
        # Read the raster data
        raster_data = src.read(1, masked=True)
        profile = src.profile
        profile.update(dtype=rasterio.int32, count=1)

        if mode == "clip":
            mask = geometry_mask(vector_layer.geometry, out_shape=src.shape, transform=src.transform, invert=False)
            features_values = np.logical_and(raster_data, mask)

        elif mode == "difference":
            mask = geometry_mask(vector_layer.geometry, out_shape=src.shape, transform=src.transform, invert=False)
            # features_values = np.where(mask, raster_data, profile['nodata'])
            features_values = np.where(mask, 1 - raster_data, raster_data)

            if level > 0:
                profile.update(nodata=0)

        elif mode == "invert":
            mask = geometry_mask(vector_layer.geometry, out_shape=src.shape, transform=src.transform, invert=False)
            features_values = np.logical_and(raster_data, mask)
            features_values = np.invert(features_values)

    with rasterio.open(output_raster, 'w', **profile) as dst:
        dst.write(features_values, 1)


def coastal_flooding_rasters_sum(input_raster, raster_to_add, raster_sum, level:int):
    """
    function that adds different layers of flooding to a tif file
    """
    with rasterio.open(input_raster) as src:
        input_raster_data = src.read(1)

    with rasterio.open(raster_to_add) as src:
        profile = src.profile
        add_raster_data = src.read(1)

    if level > 0:
        level = level +1

    add_raster_data = np.where(add_raster_data==1, level, profile['nodata'])
    sum_values = np.add(input_raster_data, add_raster_data)    

    with rasterio.open(raster_sum, 'w', **profile) as dst:
        dst.write(sum_values, 1)


def coastal_flooding_raster_subtraction(flooding_raster:str, output_raster:str, subtract = 1):
    """
    function takes the flooding layer raster and subtracts all pixels for 1
    """
    with rasterio.open(flooding_raster) as src:
        flood_data = src.read(1)
        profile = src.profile

    flood_data = np.where(flood_data, flood_data-subtract, profile["nodata"])

    with rasterio.open(output_raster, 'w', **profile) as dst:
        dst.write(flood_data, 1)


def subtract_one_from_shapefile_attribute(shapefile_path, attribute_column_name, output_shapefile):
    """
    Take an attribute/column in a shapefile and subtract 1 from it (sea level)
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


def give_tmp_file_name(input_file:str, level:int, name:str):
    """
    helper function for giving names to temporary files created in the loop function
    """
    return f"{basename_withoutext(input_file)}_{name}_{level}.shp"


def sea_levels_loop(input_file, sea_level: int):
    """
    function that extracts rising sea levels from DEM file in a loop. 
    It loops over levels for given and up to given maximum level in meters.

    1. pixels are extracted for sea levels from 0 to 10.
    2. pixels are polygonized
    3. vectors intersections with coastline checked
    4. difference is made by vector mask
    5. difference layers are merged
    """
    sea_level = sea_level + 1 # otherwise 10 will not be selected

    #TODO: add tmp folder and context manager
    with tempfile.TemporaryDirectory() as tmpdirname:
        sea_levels_simulation = f"{tmpdirname}/{basename_withoutext(input_file)}_sea_levels_simulation.tif"

        for level in tqdm(range(0, sea_level, 1)):
            # selecting by height
            selected_sea_level_file = f"{tmpdirname}/{basename_withoutext(input_file)}_sea_level_{level}.tif"
            selected_sea_level = gdal_select_pixels(input_file, selected_sea_level_file, level)

            # pixel connectivity overcome
            polygonized_pixels_file = f"{tmpdirname}/{basename_withoutext(input_file)}_sea_level_{level}.shp"
            polygonized_pixels = polygonize_raster(selected_sea_level, polygonized_pixels_file)

            # spatial relation with sea
            intersected_features_file = f"{tmpdirname}/{basename_withoutext(input_file)}_sea_level_intersected_{level}.shp"
            gpd_extract_features_by_location(polygonized_pixels_file, coastlines, intersected_features_file)

            # pixel connectivity and intersection overcome - clip the raster by mask
            intersected_tif = f"{tmpdirname}/{basename_withoutext(input_file)}_sea_level_intersected_{level}.tif"
            extract_pixels_by_mask(selected_sea_level_file, intersected_features_file, intersected_tif, mode="clip", level=level)
            
            # this part calculates the difference (and polygonizes the difference if needed)
            if level > 0:
                # difference with rasters
                difference_tif = f"{tmpdirname}/{basename_withoutext(input_file)}_sea_level_difference_{level}.tif"
                extract_pixels_by_mask(intersected_tif, reference_file, difference_tif, mode="difference", level=level)

                # polygonize - difference shall serve as reference file in next iter
                # difference_polygonized = f"{tmpdirname}/{basename_withoutext(input_file)}_difference_polygonized_{level}.shp"
                # polygonize_raster(difference_tif, difference_polygonized)

                # reference file for difference extraction *after* the differences were extracted
                reference_file = intersected_features_file

                # edge cases
                if base_layer==None:
                    base_layer = difference_tif
                
                # add a tif layer to one tif file
                coastal_flooding_rasters_sum(base_layer, difference_tif, sea_levels_simulation, level)
                base_layer = sea_levels_simulation

            else:
                # in case that level is 0, existing file is used just for invert
                reference_file = intersected_features_file # vector mask file
                difference_tif = f"{tmpdirname}/{basename_withoutext(input_file)}_sea_level_difference_{level}.tif"
                extract_pixels_by_mask(intersected_tif, reference_file, difference_tif, mode="invert", level = level)

                # this base layers is the forst to be added to main file
                base_layer = difference_tif

        # and convert the result
        polygonized_simulation = f"{tmpdirname}/{basename_withoutext(input_file)}_coastal_flooding.shp"
        polygonize_raster(sea_levels_simulation, polygonized_simulation)

        # last subtraction
        output_file = f"{basename_withoutext(input_file)}_coastal_flooding.shp"
        subtract_one_from_shapefile_attribute(polygonized_simulation, "level", output_file)

    return print(output_file, " is done")


def altitude_filter_files_list(dem_files:str, threshold_altitude=10):
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
            if (altitude_array < threshold_altitude).any():
                tiles_with_low_altitude.append(file_path)

            # Close the dataset
            dataset = None

    return tiles_with_low_altitude

def sea_level_precheck(files_list:str):
    """
    this function iterates over grid of dem files and check if dem tiles are close to sea and if they contain pixels with latitude lower than 10m.
    """
    low_dem_tiles = altitude_filter_files_list(files_list)
    return low_dem_tiles


def contours_polygons_v2():
    # gdal_contour -b 1 -i 1.0 -f "ESRI Shapefile" -fl 0 10 -p /mnt/volume-nbg1-1/shared/nikola/ESA_global_dem_30m/Copernicus_DSM_10_N38_00_W076_00/DEM/Copernicus_DSM_10_N38_00_W076_00_DEM.tif .contour_range_pol_38_76.shp
    pass    

    """
    General process of selecting areas under risk of coastal flooding up to 10m:

    1. Select dem tiles where min(value) <= 10: gdal_info
    2. create vrt for dem tiles: gdal_vrt
    2. Create contours - filter contours where value < 10: gdal_contours
    3. Buffer contours - create mask: gdal_create_buffer
    5. Select pixels where values < 10: gdal_calc.py
    4. Raster clip by mask: gdal_clip
    6. Raster clip by extent: gdal_clip
    
    """


if __name__ == "__main__":
    dem_eu_10_45_14 = "/mnt/volume-nbg1-1/satellite/eu_dem/dem10m/Copernicus_DSM_03_N45_00_E014_00_DEM.tif"
    lat, lon = 45, 14
    output_extracted = f"/home/nikola/4_north_america/GeoDataPump/data/selected_dem_pixels_10m_{lat}_{lon}.tif"
    # gdal_select_pixels(dem_eu_10_45_14, output_extracted)
    kvarner_coastline = "coastline_45_14_kvarner.shp"
    coastlines = config["NA_coastlines"]
    coastlines = kvarner_coastline
    # get_coastline_from_pg(kvarner_coastline, extent="kvarner")
    sea_levels_loop(dem_eu_10_45_14, sea_level=10)
