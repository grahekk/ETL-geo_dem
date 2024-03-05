from test_pixel_utils import get_pixel_value_from_geotiff, get_tif_pixel_value, land_cover_percentage, clip_by_extent_gdal, query_pixel_value_rasterio
from test_create_coordinates import read_coordinates, enrich_existing_shp
import yaml
import pandas as pd
import geopandas as gpd
import sys
import time
from statistics import mean
import multiprocessing
import random

sys.path.append("/home/nikola/4_north_america/scripts/")
import settings
config = settings.get_config()

testing_file_usa = config["usa_dem_vrt"]
testing_file_cad = config["canada_dem_vrt"]
testing_file_global = config["nasa_global_dem_vrt"]
esa_global_dem_90_vrt = config["esa_global_dem_90_vrt"]
na_aspect_categorized = config["NA_aspect_categorized"]
# test_coords_path = "/home/nikola/grid_clipped.shp"
test_coords_path = "/home/nikola/grid_clipped_enrich.shp"
# test_coords_path = "/home/nikola/test_coordinates.csv"
test_coords_path = "/home/nikola/continents_grid_clipped_enrich.shp"
test_coords_path = "/home/nikola/4_north_america/tests/NA_test_coordinates.csv"

test_coords = read_coordinates(test_coords_path)


def test_points(testing_file, test_coords, test_func = ""):
    """
    Test a list of coordinates for pixel values in a GeoTIFF file and measure execution time.

    This function takes a GeoTIFF file path and a list of coordinates (latitude, longitude, city, country) to
    test. It extracts pixel values from the GeoTIFF file at each coordinate, records execution times, and returns
    the results as a list of tuples.

    Parameters:
    - testing_file (str): The file path of the GeoTIFF to be tested.
    - test_coords (list of tuple): A list of tuples containing (latitude, longitude, city, country) for testing.

    Returns:
    - results (list of tuple): A list of tuples with the following format:
      (latitude, longitude, city, country, pixel_value, execution_time).
    """
    results = []
    fun_times = []
    total_count = len(test_coords)
    for count, (lat, lon, city, country) in enumerate(test_coords):

        # results.append(get_pixel_value_from_geotiff(testing_file, lat, lon))
        start_time = time.time()
        result = get_tif_pixel_value(lat, lon, testing_file)
        fun_time = round(time.time() - start_time, 2)
        results.append((lat, lon, city, country, result, fun_time))
        fun_times.append(fun_time)

        # pretty print out for task completeness
        if count % 10 == 0:
            mean_time = mean(fun_times)
            time_left = round(((total_count-count)*mean_time)/60,2)

        percentage = round(count / total_count * 100, 2)
        sys.stdout.write('\r')
        sys.stdout.write(f'task for {testing_file} completeness: {percentage:.2f}% time left: {time_left:.2f} mins')
        sys.stdout.flush()

    print(f"\nTesting done for {testing_file}!")
    return results


def process_coordinates(coord, testing_file, total_count, fun_times, mean_time, shared_count, num_processes, lock):
    lat, lon, city, country = coord
    start_time = time.time()
    result = get_tif_pixel_value(lat, lon, testing_file)
    fun_time = round(time.time() - start_time, 2)

    with lock:
        count = shared_count.value
        fun_times.append(fun_time)
        if count % 10 == 0:
            mean_time.value = mean(fun_times)
            time_left = round((((total_count-count)/num_processes)*mean_time.value)/60,2)
        shared_count.value += 1

    percentage = round(count / total_count * 100, 2)
    time_left = round(((total_count - count) * mean_time.value) / 60, 2)
    sys.stdout.write(f'\rTask for {testing_file} completeness: {percentage:.2f}% time left: {time_left:.2f} mins')
    sys.stdout.flush()
    return (lat, lon, city, country, result, fun_time)


def test_points_with_multiprocessing(testing_file, test_coords, test_func=get_tif_pixel_value):
    """
    Test a list of coordinates for pixel values in a GeoTIFF file and measure execution time.

    This function takes a GeoTIFF file path and a list of coordinates (latitude, longitude, city, country) to
    test. It extracts pixel values from the GeoTIFF file at each coordinate, records execution times, and returns
    the results as a list of tuples.

    Parameters:
    - testing_file (str): The file path of the GeoTIFF to be tested.
    - test_coords (list of tuple): A list of tuples containing (latitude, longitude, city, country) for testing.

    Returns:
    - results (list of tuple): A list of tuples with the following format:
      (latitude, longitude, city, country, pixel_value, execution_time).
    """
    results = []
    fun_times = []
    total_count = len(test_coords)

    with multiprocessing.Manager() as manager:
        lock = manager.Lock()
        shared_count = manager.Value('i', 0)
        mean_time = manager.Value('d', 0.00)
        fun_times = manager.list()

        # Number of processes to use for multiprocessing
        num_processes = multiprocessing.cpu_count()
        
        with multiprocessing.Pool(num_processes) as pool:
            result_list = pool.starmap(process_coordinates, [(coord, testing_file, total_count, fun_times, mean_time, shared_count, num_processes, lock) for coord in test_coords])

    results.extend(result_list)
    # fun_times.extend([item[-1] for item in result_list])

    print(f"\nTesting done for {testing_file}!")
    return results


def process_test_points(input_file, test_coords):
    results = test_points(input_file, test_coords)
    df = pd.DataFrame(results, columns=["lat", "lon", "city", "country", "value", "fun_time"])
    return df


def save_test_results_multiprocess():
    # Initialize an empty DataFrame
    df = pd.DataFrame(columns=["lat", "lon", "city", "country", "value", "fun_time"])

    # Define input files for test_points
    testing_files = [testing_file_usa, testing_file_cad, testing_file_global]
    # testing_files = [testing_file_usa, testing_file_cad]

    # Create a pool of worker processes
    with multiprocessing.Pool(processes=3) as pool:
        # Process test_points in parallel
        dfs = pool.starmap(process_test_points, [(file, test_coords) for file in testing_files])
    # Concatenate the DataFrames from multiprocessing
    df = pd.concat([df, *dfs], ignore_index=True)

    df.to_csv("test_results.csv", index=False)

    # convert to shp
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.lon, df.lat), crs="EPSG:4329")
    gdf.to_file("test_results.shp")


def save_test_results(filename, test_file, test_coords):
    """
    run one of methods over testing file and call testing coords to test the data.
    the results are saved into csv and shp so that they can be interpreted easire using qgis gui.

    Parameters:
    - filename (str): file where results are saved
    - test_file (str): The file path of the GeoTIFF to be tested.
    - test_coords (list of tuple): A list of tuples containing (latitude, longitude, city, country) for testing.
    """
    # results = test_points(test_file, test_coords)
    results = test_points_with_multiprocessing(test_file, test_coords)
    df = pd.DataFrame(results, columns=["lat", "lon", "city", "country", "value", "fun_time"])
    df.to_csv(f"{filename}_test_results.csv", index=False)
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.lon, df.lat), crs="EPSG:4329")
    gdf.to_file(f"{filename}_test_results.shp")
    return print(f"Testing results saved into {filename}_test_results.csv")


def test_land_cover():
    input_path = config["esa_world_cover_vrt"]
    #lat, lon = 45.31067674891032, 14.38130187962704
    lat, lon = 45.80078805790288, 15.995128030974094 # križanje
    lat, lon = 45.31329434266051, 14.474612677074116
    lat, lon = 45.31585277138846, 14.455254506203037 # more
    lat, lon = 45.88743928524575, 13.131780970967833 # italyy
    lat, lon = 40.75089442096933, -74.20510176323846 # nueva york
    lat, lon = 50.71860137019576, -90.77808839956987 # canada ontario
    result = land_cover_percentage(lat, lon, input_path)
    print(result)


def test_data_esa_dem_continents(esa_global_dem_90_vrt, singleprocess_test_points):
    enrich_existing_shp("/home/nikola/continents_grid_clipped.shp")
    test_coords_path = "/home/nikola/continents_grid_clipped_enrich.shp"
    test_coords = read_coordinates(test_coords_path)
    singleprocess_test_points("esa_world_dem", esa_global_dem_90_vrt, test_coords)


def test_data_aspect_categories():
    """
    Function flow that calls the `save_test_results` function.
    Test coords are taken from the path provided in the variable `test_coords` above.
    """
    test_coords = read_coordinates(test_coords_path)
    na_aspect = config["NA_aspect"]
    # na_aspect_categorized = config["NA_aspect_categorized"]
    save_test_results("na_aspect_categorized", na_aspect_categorized, test_coords)


def test_data_slope_rescaled():
    """
    Function flow that calls the `save_test_results` function.
    Test coords are taken from the path provided in the variable `test_coords` above.
    """
    test_coords = read_coordinates(test_coords_path)
    na_slope_rescaled = config["NA_slope_rescaled"]
    save_test_results("na_slope_rescaled", na_slope_rescaled, test_coords)


if __name__ == "__main__":
    # singleprocess_test_points("global", testing_file_global)
    # test_land_cover()
    # esa_dem_continents_flow(esa_global_dem_90_vrt, singleprocess_test_points)

    # # the test
    # test_coords = read_coordinates(test_coords_path)
    # save_test_results("na_aspect_categorized", na_aspect_categorized, test_coords)

    # # cutting test
    lon, lat = 35.133991822016135, -119.71258190229405 #san andreas fault
    # lon, lat = 38.69014,-92.26247
    # lon, lat = 62.8487,-128.9834 #canada
    lon, lat = 37.59747, -113.87625
    lon, lat = -129.2448788, 62.4509572
    lat, lon = 62.7105955, -128.7794029
    lat, lon = 33.606916,-115.266720
    # lat, lon = 33.611793,-115.262013
    # lat, lon = 33.603675,-115.272724
    # lat, lon = 33.675886,-115.353531
    # lat, lon = 33.673413,-115.368096
    # lat, lon = 33.672987,-115.365180
    # lat, lon = 33.683578,-115.367670
    # # slope_raster = "/mnt/volume-nbg1-1/shared/nikola/NA_slope.tif"
    # # slope_cut = "./slope_result.tif"
    
    raster_in = config["NA_slope_rescaled"]
    raster_cut = "./slope_result.tif"
    clip_by_extent_gdal(lat, lon, raster_cut, raster_in)

    # test_data_aspect_categories()
    # test_data_slope_rescaled()

    # value = query_pixel_value(aspect_raster, lat, lon) #rasterio novi
    # value = get_pixel_value_from_geotiff(aspect_raster, lat, lon) #rasterio početni
    # print(value, " rasterio")
    # value = get_tif_pixel_value(lat, lon, aspect_raster) #gdal
    # print(value, " gdal")