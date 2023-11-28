import os
import sys
import subprocess
import yaml
import time
import datetime
import psycopg2
import re
import numpy as np
from statistics import mean

# sys.path.append("/home/nikola/4_north_america/scripts/")

# load configuration
import settings
config = settings.get_config()
conn_parameters = settings.get_conn_parameters()

#import python path from qgis so that gdal from osgeo can be loaded
sys.path.append('/usr/lib/python3/dist-packages')
from osgeo import gdal 


def absolute_file_paths(directory:str, extension:str): # ensure absolute paths
    """
    Returns generator object - a list of absolute file paths for given directory.

    Args:
        directory (str): path to directory where files are.
        extension (str): extension of wanted files
    """
    for dirpath,_,filenames in os.walk(directory):
        for f in filenames:
            if f.endswith(extension):
                yield os.path.abspath(os.path.join(dirpath, f))


def split_list(input_list:list, split_count:int):
    """
    Split a list into smaller lists with an approximately equal number of items.

    Args:
        input_list (list): The list to be split.
        split_count (int): The number of smaller lists to create.

    Returns:
        list of lists: A list of smaller lists, each containing an approximately
        equal number of items from the input_list.
    """
    # Calculate the length of each smaller list
    sublist_length = len(input_list) // split_count

    # Check if split_count is greater than the length of the list
    if sublist_length == 0:
        return [input_list]

    # Use list comprehension to split the input_list into smaller lists
    result = [input_list[i:i + sublist_length] for i in range(0, len(input_list), sublist_length)]
    return result


def geocellid_by_continents(conn_parameters, continent = "North America"):
    """
    Filter spatial data from a PostgreSQL table by a specified continent or country.

    This function connects to a PostgreSQL database, executes a stored procedure to obtain
    geocell IDs associated with the specified continent or country, and returns the filtered results.

    Args:
        table_to_filter (str, optional): The name of the table containing spatial data in the PostgreSQL database.
            Defaults to "esa_grid".
        continent (str, optional): The continent or country by which spatial data will be filtered.
            Defaults to "North America".
    """
    # Connect to the PostgreSQL database and fetch the geocell IDs
    with psycopg2.connect(**conn_parameters) as conn:
        with conn.cursor() as cursor: 
            cursor.callproc('osm.get_geocellids_in_continent', [continent])
            result = cursor.fetchall()
    conn.close()
    geocell_ids = [row[0] for row in result] #convert list of tuples containing string to list of strings
    return geocell_ids

def geocellid_by_country(conn_parameters, country = "United States"):
    """
    Filter spatial data from a PostgreSQL table by a specified continent or country.

    This function connects to a PostgreSQL database, executes a stored procedure to obtain
    geocell IDs associated with the specified continent or country, and returns the filtered results.

    Args:
        table_to_filter (str, optional): The name of the table containing spatial data in the PostgreSQL database.
            Defaults to "esa_grid".
        country (str, optional): The country by which spatial data will be filtered.
            Defaults to "North America".
    """
    # Connect to the PostgreSQL database and fetch the geocell IDs
    with psycopg2.connect(**conn_parameters) as conn:
        with conn.cursor() as cursor: 
            cursor.callproc('osm.get_geocellids_in_country', [country])
            result = cursor.fetchall()
    conn.close()
    geocell_ids = [row[0] for row in result] #convert list of tuples containing string to list of strings
    return geocell_ids


def filter_by_geocellid(links, coc = "continent", continent = "North America", country = "United States"):
    """
    Pure function to filter a list of links based on geocell IDs extracted from their filenames and a filter list of geocell IDs.
    Function uses helper function `filter_by_continents()` that extreacts geocell IDs from pg db.

    Args:
        links (list): A list of URLs or file paths.
        coc (str): (coc - Continent or country) either "continent" or "country" - a switch whether to filter by continent or a country
        continent (str): A name of continent by which to filter extent of geocellid's

    Returns:
        list: A filtered list of links that match the geocell IDs specified in the filter list for "North America."
    """
    filtered_links = []
    if coc == "continent":     
        filter_list = geocellid_by_continents(conn_parameters, continent = continent)
    elif coc=="country":
        filter_list = geocellid_by_country(conn_parameters, country = country)
    
    for link in links:
        file_name = os.path.basename(link)
            # match something like "N23_00_W123"
        geocell_id_match = re.search(r'_(S|N)(\d+)_00_(W|E)(\d+)', file_name)
        if geocell_id_match:
            lat_direction = geocell_id_match.group(1)
            lat_value = geocell_id_match.group(2)
            lon_direction = geocell_id_match.group(3)
            lon_value = geocell_id_match.group(4)
                
                # Construct the geocell ID based on the extracted values
            geocell_id = f"{lat_direction}{lat_value}{lon_direction}{lon_value}"

            if geocell_id in filter_list:
                filtered_links.append(link)
    return filtered_links


def create_vrt(output_vrt:str, input_path:str, vrt_chunks=1, use_prefix = '.tif'):
    """
    Create a VRT (Virtual Raster) file from a list of input raster files.

    Args:
        output_vrt (str): The output VRT file path.
        input_path (str): Directory of input raster files to include in the VRT.
        vrt_chunks (int): number of VRT files/chunks to be created
        use_prefix (str): this is prefix of all the files to be used to create VRT
    """
    # os.makedirs(os.path.dirname(output_vrt), exist_ok = True)
    input_files = absolute_file_paths(input_path, use_prefix) # returns generator object
    files_list = sorted(list(input_files))
    
    if vrt_chunks <= 1: #the usual case when only one vrt is created
        input_files_str = " ".join(files_list)
        create_single_vrt(output_vrt, input_files_str)
    else:
        input_lists = split_list(files_list, vrt_chunks)
        for i, chunk in enumerate(input_lists):
            input_files_str = " ".join(chunk)
            chunk_output_vrt = f"{output_vrt}_{i}.vrt"
            create_single_vrt(chunk_output_vrt, input_files_str)
            print(f"VRT chunk '{chunk_output_vrt}' created successfully.")


def create_single_vrt(output_vrt:str, input_files_str:str):
    """    
    Build a single VRT using the gdalbuildvrt command in shell
    """
    command = f"gdalbuildvrt {output_vrt} {input_files_str} -separate"

    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error creating VRT: {e}")


def create_overviews(input_vrt:str, overview_levels, resampling_method="average"):
    """
    Create overview files for a VRT (Virtual Raster).

    Args:
        input_vrt (str): The input VRT file path.
        overview_levels (list): List of overview levels to generate.
        resampling_method (str): Resampling method (e.g., 'average', 'nearest', 'bilinear').
    """
    overview_levels_str = ",".join(map(str, overview_levels))

    # Generate overviews using the gdaladdo command
    command = f"gdaladdo -r {resampling_method} {input_vrt} {overview_levels_str}"
    
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Overviews created for '{input_vrt}' with levels: {overview_levels_str}")
    except subprocess.CalledProcessError as e:
        print(f"Error creating overviews: {e}")


def convert_float32_to_int16(input_dir:str, convert_output_type="Int16"):
    """
    Convert and rename TIFF files with `Float32` data type to `Int16` in the specified directory.

    This function scans the provided directory for TIFF files, checks if their data type is Float32,
    and if so, renames the original file by appending "_float32" to the filename and converts the
    data type to `Int16`, saving the result with the original filename (without "_float32").

    Parameters:
    - input_dir (str): The path to the directory containing TIFF files to process.
    - convert_output_type (str): Type of number format written in memory (default is `Int16`)

     Requirements:
    - GDAL command-line tools must be installed on your system for this function to work.
    """
    # List all TIFF files in the directory
    tif_files = [f for f in os.listdir(input_dir) if f.endswith(".tif")]

    for tif_file in tif_files:
        try:
            tif_path = os.path.join(input_dir, tif_file)
            
            # Open the TIFF file using GDAL
            cmd = ["gdalinfo", tif_path]
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            output, _ = process.communicate()

            # Check if the output indicates a Float32 data type
            if "Float32" in output:
                # Construct the output filename for the Int16 version
                output_file = os.path.join(input_dir, os.path.splitext(tif_file)[0] + f'_{convert_output_type}.tif')
                
                # # Construct the filename for the original Float32 version with "_float32" appended
                # float32_filename = os.path.join(input_dir, os.path.splitext(tif_file)[0] + f'_{convert_output_type}.tif')
                
                # # Rename the original file to include "_float32" in the filename
                # os.rename(tif_path, float32_filename)
                
                # Use gdal_translate through subprocess to convert the data type to Int16 and save it with new filename
                cmd = ["gdal_translate", "-ot", convert_output_type, tif_path, output_file]
                process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                process.communicate()

                print(f"Converted {tif_path} to {output_file}")

        except Exception as e:
            print(f"Error processing {tif_file}: {str(e)}")

    print("Conversion and renaming completed.")


def geofilter_paths_list(files_dir:str):
    """
    Filter and return a sorted list of file paths from a directory based on geocell IDs.

    Args:
        files_dir (str): The directory containing the files to filter.

    Returns:
        list: A sorted list of file paths filtered based on geocell IDs.
    """
    input_files = absolute_file_paths(files_dir, ".tif") # returns generator object
    filtered_files = filter_by_geocellid(input_files, coc = "country", country = "United States")
    files_list = sorted(list(filtered_files))
    return files_list


def transform_raster(input_raster:str, output_raster:str, transform_type="slope"):
    """    
    Compute raster transformation using the gdaldem command from shell
    Parameters:
        - input_raster (str): The input raster file path.
        - output_raster (str): The output raster file path.
        - transform_type (str): The type of transformation (e.g., "slope" or "aspect").
    """
    if transform_type=="slope":
        command = f"gdaldem slope {input_raster} {output_raster} -p -b 1 -s 1.0"
    elif transform_type == "aspect":
        command = f"gdaldem aspect {input_raster} {output_raster} -zero_for_flat -b 1"
    else:
        raise ValueError("Invalid transform_type. Use 'slope' or 'aspect'.")
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"{transform_type} done for {output_raster}!")
    except subprocess.CalledProcessError as e:
        print(f"Error transforming raster: {e}")


def categorize_aspect(input_aspect_path, output_category_path, chunk_size=(1000, 1000)):
    """
    Categorize aspect values into 8 categories based on predefined ranges.

    Parameters:
    - input_aspect_path (str): Path to the input aspect raster file.
    - output_category_path (str): Path to the output categorized raster file.
    - chunk_size (tuple): Size of the processing chunks. Default is (1000, 1000).

    Returns:
    - None
    """
    # Open the input aspect raster
    start_time = time.time()
    input_ds = gdal.Open(input_aspect_path)
    if input_ds is None:
        raise Exception(f"Failed to open the input aspect raster: {input_aspect_path}")

    # Get aspect band
    aspect_band = input_ds.GetRasterBand(1)
    
    # Get raster dimensions
    rows, cols = input_ds.RasterYSize, input_ds.RasterXSize

    # Define aspect categories
    categories = [
        (0.0001, 22.4999),
        (22.5, 67.4999),
        (67.5, 112.499),
        (112.5, 157.4999),
        (157.5, 202.499),
        (202.5, 247.499),
        (247.5, 292.499),
        (292.5, 337.499),
        (337.5, 359.999),
    ]

    # Create the output raster
    driver = gdal.GetDriverByName('GTiff')
    output_ds = driver.Create(output_category_path, cols, rows, 1, gdal.GDT_Byte)
    output_ds.SetProjection(input_ds.GetProjection())
    output_ds.SetGeoTransform(input_ds.GetGeoTransform())

    # Process the raster in chunks, so that the process isn't memory intensive
    count = 0
    fun_times = []
    # total_count = max(range(0, rows, chunk_size[0]))
    # total_count = (rows/chunk_size[0]) * (cols/chunk_size[1])
    total_count = rows/chunk_size[0]

    for i in range(0, rows, chunk_size[0]):
        start_time = time.time()
        for j in range(0, cols, chunk_size[1]):
            chunk_rows = min(chunk_size[0], rows - i)
            chunk_cols = min(chunk_size[1], cols - j)

            # Read aspect data as a numpy array for the current chunk
            aspect_data = aspect_band.ReadAsArray(j, i, chunk_cols, chunk_rows)

            # Create an empty array for the categorized aspect chunk
            categorized_aspect_chunk = np.zeros_like(aspect_data, dtype=np.uint8)

            # Categorize aspect values for the current chunk
            for k, (lower, upper) in enumerate(categories, start=1):
                if lower < upper:
                    mask = np.logical_and(aspect_data >= lower, aspect_data < upper)
                else:
                    mask = np.logical_or(aspect_data >= lower, aspect_data < upper)
                categorized_aspect_chunk[mask] = k
            
            # category 9 and 1 are halfs of one whole/same category, 1 (north). category 0 are null values
            categorized_aspect_chunk[categorized_aspect_chunk == 9] = 1
            # Write the categorized aspect chunk to the output raster
            output_band = output_ds.GetRasterBand(1)
            output_band.WriteArray(categorized_aspect_chunk, j, i)

        # print progress
        fun_time = time.time() - start_time 
        fun_times.append(fun_time)

        if count % 10 == 0:
            mean_time = mean(fun_times)
            time_left = round(((total_count-count)*mean_time)/60,2)

        count = count+1
        percentage = round(count / total_count * 100, 2)
        sys.stdout.write('\r')
        sys.stdout.write(f'task for {__name__} completeness: {percentage:.2f}% time left (d:hh:mm.ss): {datetime.timedelta(seconds=time_left)}')
        sys.stdout.flush()

    # Close datasets
    input_ds = None
    output_ds = None


def gdal_build_vrt(files_list:list, output_vrt:str, resample_method = "bilinear"):
    """
    Create a Virtual Raster (VRT) dataset from a list of raster files.

    Args:
        files_list (list): The list of paths to the input raster files
        output_vrt (str): The path of the output .vrt file
        resample_method (str): The resampling method to use when creating the VRT. Options: 'cubic', 'average', 'bilinear', 'nearest'.

    This function constructs a Virtual Raster (VRT) dataset using GDAL tools to mosaic multiple raster files
    into a single virtual dataset. It uses the specified 'raster_files' path to locate the input raster files
    and 'output_vrt' to define the output VRT filename. You can also specify the resampling method.

    The VRT options include adding an alpha band and displaying nodata values.

    example:
        raster_files = "esa_global_dem_90"
        resample_method = "cubic, average"
    """
    assert resample_method in ["cubic", "average", "bilinear", "nearest"], "Invalid resample_method"

    # input_files = absolute_file_paths(raster_files_dir, ".tif") # returns generator object
    # files_list = sorted(list(input_files))

    start_time = time.time()
    vrt_options = gdal.BuildVRTOptions(resampleAlg=resample_method, addAlpha=True, hideNodata=False)
    my_vrt = gdal.BuildVRT(output_vrt, files_list, options=vrt_options)
    my_vrt = None
    fun_time = round(time.time()-start_time, 2)
    return print(f"Gdal building vrt done for file: {output_vrt} in time: {datetime.timedelta(seconds=fun_time)}")


def create_vrt_ovr_flow(input_raster_directory:str, output_vrt:str, vrt_chunks=1, use_prefix='.tif'):
    """
        Create a Virtual Raster (VRT) file and its associated overviews from a dataset of rasters.

        This function performs the following tasks:
        1. Takes a directory containing input raster files specified in the 'config' dictionary.
        2. Creates a Virtual Raster (VRT) file that combines the input rasters into a single dataset.
        3. Generates overviews for the VRT file with specified overview levels.
        
        Args:
            input_raster_directory (str): The key in the 'config' dictionary representing the input raster directory.
            output_vrt (str): The key in the 'config' dictionary representing the output VRT file.
            vrt_chunks (int): Number of vrt files to be created from number of files in a directory
            use_prefix (str): Prefix on each file to be used for creation of VRT - example: "_float32.tif"
    """
    input_path = config[input_raster_directory]
    output_vrt = config[output_vrt]
    # overview_levels = [2, 4, 8, 16]
    overview_levels = [256]

    create_vrt(output_vrt, input_path, vrt_chunks, use_prefix)
    create_overviews(output_vrt, overview_levels)


# Create flows:
# deprecated flows
def canada_dem_vrt_ovr_flow():
    """
    conversion was needed since gdalbuildvrt cannot build vrt from different types of tif files.
    Some tiff files are `int16` while others are of `float32` type.
    """
    convert_float32_to_int16(config["canada_dem_dir"])
    return create_vrt_ovr_flow("canada_dem_dir", "canada_dem_vrt", use_prefix=".tif") # vrt from files downloaded


def canada_WMS_dem_vrt_ovr_flow():
    """
    Unlike `canada_dem_vrt_ovr_flow`, this flow downloads dem tiles from WMS. It acts like reserve
    this flow first creates chunks of vrt's and then creates one main vrt
    """
    create_vrt("canada_dem_WMS_dir", "canada_dem_WMS_dir", vrt_chunks=3, use_prefix=".tif")
    return create_vrt_ovr_flow("canada_dem_WMS_dir", "canada_dem_WMS_vrt", vrt_chunks=1, use_prefix=".vrt")


def NASA_global_dem_vrt_ovr_flow():
    """
    Slightly different flow than the two ones before.
    This flow creates chunks of .vrt files because there is too much files to be read all at once into `gdalbuildvrt`.
    Then these chunks are again iteraded over so that one main VRT is created. 
    There is no need to create .ovr file since the downloaded data from nasa already contains calculated ovr.
    Nasa's global dem `_num.tif` files are in byte format while `_dem.tif` files are int16.
    """
    # return create_vrt_ovr_flow("nasa_global_dem_dir", "nasa_global_dem_vrt", 30)
    
    input_dir = config["nasa_global_dem_dir"]
    output_vrt = config["nasa_global_dem_vrt"]
    input_vrts_dir = config['nasa_global_dem_vrt_dir']
    os.makedirs(input_vrts_dir, exist_ok=True)
    create_vrt(output_vrt, input_dir, 30, "dem.tif")
    create_vrt(output_vrt, input_vrts_dir, 1, ".vrt")
    return 

