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
import rasterio

# load configuration
import settings
config = settings.get_config()
conn_parameters = settings.get_conn_parameters()
schema = settings.get_schema()

from .model_data import get_geocellid
# from .pipeline_transform_sea_level import basename_withoutext

#import python path from qgis so that gdal from osgeo can be loaded
sys.path.append('/usr/lib/python3/dist-packages')
# from osgeo import gdal 


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
    
    Returns:
        geocell_ids (list): list of strings
    """
    # Connect to the PostgreSQL database and fetch the geocell IDs
    with psycopg2.connect(**conn_parameters) as conn:
        with conn.cursor() as cursor: 
            cursor.callproc(f'{schema}.get_geocellids_in_continent', [continent])
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
            cursor.callproc(f'{schema}.get_geocellids_in_country', [country])
            result = cursor.fetchall()
    conn.close()
    geocell_ids = [row[0] for row in result] #convert list of tuples containing string to list of strings
    return geocell_ids


def filter_by_geocellid(files, by:str, c_name = "United States", extent = (-25, 35, 35, 72), regex_match = r'_(S|N)(\d+)_00_(W|E)(\d+)'):
    """
    Filter a list of files based on geocell IDs extracted from their filenames and a filter list of geocell IDs.
    The function uses helper functions to extract geocell IDs from a PostgreSQL database (for continents)
    or to generate geocell IDs based on a specified country or extent.

    Args:
        files (list): A list of URLs or file paths.
        by (str): Specifies whether to filter by "continent," "country," or "extent."
        c_name (str, optional): Name of the continent, country, or extent depending on the 'by' parameter. Defaults to "United States."
        extent (tuple, optional): A tuple specifying the extent (latitude and longitude) for filtering by extent. Defaults to (-25, 35, 35, 72).
        regex_match (str, optional): Regular expression pattern to extract geocell coordinates from file names. Defaults to r'_(S|N)(\d+)_00_(W|E)(\d+)'.

    Returns:
        list: A filtered list of files that match the geocell IDs specified in the filter list for the specified continent, country, or extent.
    """

    if by == "continent":     
        geocell_filter_list = geocellid_by_continents(conn_parameters, continent = c_name)
    elif by == "country" or by == "extent":
        geocell_filter_list = get_geocellid(c_name = c_name, by = by, extent = extent)

    filtered_files = geocell_regex_match(files, geocell_filter_list, regex_match)
    return filtered_files


def geocell_regex_match(files:list, geocell_filter_list:list, regex_match):
    """
    Filter a list of files based on geocell coordinates extracted using a regular expression.

    Parameters:
    - files (list): List of file paths to be filtered.
    - geocell_filter_list (list): List of geocell coordinates to filter the files.
    - regex_match (str): Regular expression pattern to extract geocell coordinates from file names.

    Returns:
    - list: A filtered list of file paths that match the specified geocell coordinates.
    """
    filtered_files = []
    for file in files:
        file_name = os.path.basename(file)
        # match something like "N23_00_W123"
        # regex_match = r'_(S|N)(\d+)_00_(W|E)(\d+)'

        geocell_id_match = re.search(regex_match, file_name)
        if geocell_id_match:
            lat_direction = geocell_id_match.group(1)
            lat_value = geocell_id_match.group(2)
            lon_direction = geocell_id_match.group(3)
            lon_value = geocell_id_match.group(4)
                
            # Construct the geocell ID based on the extracted values
            geocell_id = f"{lat_direction}{lat_value}{lon_direction}{lon_value}"

            if geocell_id in geocell_filter_list:
                filtered_files.append(file)

    return filtered_files


def geocellid_from_file_name(file_path:str, regex_match = r'_(S|N)(\d+)_00_(W|E)(\d+)'):
    """
    Take file name and extract a geocellid from it, if it matches.

    Args:
        file_path(str): a path to file
    """
    file_name = os.path.basename(file_path)


    geocell_id_match = re.search(regex_match, file_name)
    if geocell_id_match:
        lat_direction = geocell_id_match.group(1).upper()
        lat_value = geocell_id_match.group(2)
        lon_direction = geocell_id_match.group(3).upper()
        lon_value = geocell_id_match.group(4)

        if "USGS" in file_name:
            lat_value = int(lat_value)-1
            
        # Construct the geocell ID based on the extracted values
        geocell_id = f"{lat_direction}{lat_value}{lon_direction}{lon_value}"
    
    return geocell_id

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


def extract_raster_features_rasterio(input_path:str, output_path:str, feature:int):
    with rasterio.open(input_path) as src:
        # Read the raster data as a NumPy array
        data = src.read(1)

        profile = src.profile
        # Apply the mask to set non-matching pixels to nodata_value
        filtered_data = np.where(data == feature, 1, profile["nodata"])

        # Get metadata from the source raster
        metadata = src.meta.copy()

    metadata.update(
        dtype=rasterio.uint8,  # Set the data type as needed
        count=1  # Number of bands
    )

    # Create a new raster TIF file and write the filtered data
    with rasterio.open(output_path, 'w', **metadata) as dst:
        dst.write(filtered_data, 1)


def extract_values_gdal(input_path:str, output_path:str, feature:int):
    """
    Extracts and keeps all values that are equal to `feature` from an input raster TIF file.

    Parameters:
    - input_path (str): The file path to the input raster TIF file.
    - output_path (str): The file path to the output raster TIF file containing extracted values.
    - feature (int): Value to be filtered out of raster .

    Raises:
    - Exception: If unable to open the input raster file or create the output raster file.

    Notes:
    - The function reads the input raster, extracts values equal to 10, and sets all other values to 0.
    - The output raster will have the same geotransform and projection as the input raster.

    Example:
    >>> input_raster_path = "input_raster.tif"
    >>> output_raster_path = "output_extracted_raster.tif"
    >>> extract_values(input_raster_path, output_raster_path)
    """
    # Open the input raster file
    input_ds = gdal.Open(input_path)

    if input_ds is None:
        raise Exception(f"Unable to open input raster file: {input_path}")

    # Get raster band
    band = input_ds.GetRasterBand(1)

    # Read the raster data
    data = band.ReadAsArray()
    # extract
    extracted_data = data.copy()
    extracted_data[data != feature] = 0  # Set values not equal to feature to 0
    output_ds = gdal.GetDriverByName('GTiff').Create(output_path, input_ds.RasterXSize, input_ds.RasterYSize, 1, band.DataType, options=['COMPRESS=DEFLATE', 'PREDICTOR=2', 'TILED=YES'])

    if output_ds is None:
        raise Exception(f"Unable to create output raster file: {output_path}")

    # Set 0 as nodata value
    output_band = output_ds.GetRasterBand(1)
    output_band.SetNoDataValue(0)

    # Set the geotransform and projection from the input raster
    output_ds.SetGeoTransform(input_ds.GetGeoTransform())
    output_ds.SetProjection(input_ds.GetProjection())

    # Write the extracted data to the new raster
    output_band = output_ds.GetRasterBand(1)
    output_band.WriteArray(extracted_data)

    # Close the datasets
    input_ds = None
    output_ds = None


def geofilter_paths_list(files_dir:str, by = "continent", c_name = "Europe", extent = (-25, 35, 35, 72)):
    """
    Filter and return a sorted list of file paths from a directory based on geocell IDs.

    Args:
        files_dir (str): The directory containing the files to filter.
        by (str): (coc - Continent or country) either "continent" or "country" - a switch whether to filter by continent or a country
        c_name (str): a name of continent or a country to be queried
        extent (tuple): tuple of coordiantes (min_x, min_y, max_x, max_y)

    Returns:
        list: A sorted list of file paths filtered based on geocell IDs.
    """
    assert by in ("country", "continent", "extent")
    assert len(extent) == 4

    input_files = absolute_file_paths(files_dir, ".tif") # returns generator object
    filtered_files = filter_by_geocellid(input_files, by = by, c_name = c_name, extent=extent)
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


def rescale_raster(input_raster:str, output_path:str, band=1, new_min=0, new_max=90, nodata=None):
    """
    Rescales pixel values in a raster layer and writes the result to a new raster file.

    Parameters:
        input_raster (str): Path to the input raster file.
        output_path (str): Path to the output raster file.
        band (int, optional): The band index to rescale. Default is 1.
        new_min (float, optional): The minimum value for rescaling. Default is 0.
        new_max (float, optional): The maximum value for rescaling. Default is 90.
        nodata (float, optional): The nodata value to be set in the output raster. Default is None.

    Returns:
        None: The function writes the rescaled raster to the specified output path.
    """
    
    # Open the raster dataset, read band as np array
    raster_dataset = gdal.Open(input_raster, gdal.GA_ReadOnly)
    band_data = raster_dataset.GetRasterBand(band)
    data_array = band_data.ReadAsArray()

    # Apply the rescaling
    scaled_data = np.clip((data_array - np.nanmin(data_array)) / (np.nanmax(data_array) - np.nanmin(data_array)) * (new_max - new_min) + new_min, new_min, new_max)

    # Set nodata values
    if nodata is not None:
        scaled_data[np.isnan(data_array)] = nodata

    # Create a new raster dataset for the scaled data
    new_dataset = gdal.GetDriverByName('GTiff').Create(output_path, raster_dataset.RasterXSize, raster_dataset.RasterYSize, 1, band_data.DataType)
    new_dataset.SetGeoTransform(raster_dataset.GetGeoTransform())
    new_dataset.SetProjection(raster_dataset.GetProjection())

    # Write the scaled data to the new raster dataset
    new_dataset.GetRasterBand(1).WriteArray(scaled_data)

    # Set nodata values in the new raster dataset
    if nodata is not None:
        new_dataset.GetRasterBand(1).SetNoDataValue(nodata)

    # Close the datasets
    raster_dataset = None
    new_dataset = None


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
    Make sure to provide correct path to `osgeo` to import `gdal`.

    The VRT options include adding an alpha band and displaying nodata values.

    example:
        raster_files = "esa_global_dem_90"
        resample_method = "cubic, average"
    """
    assert resample_method in ["cubic", "average", "bilinear", "nearest"], "Invalid resample_method"

    vrt_options = gdal.BuildVRTOptions(resampleAlg=resample_method, addAlpha=False, hideNodata=True)
    my_vrt = gdal.BuildVRT(output_vrt, files_list, options=vrt_options)
    my_vrt = None
    # return print(f"Gdal building vrt done for file: {output_vrt} in time: {datetime.timedelta(seconds=fun_time)}")


def gdal_clip(input_tif:str, output_tif:str, projWin = [-75.3, 5.5, -73.5, 3.7]):
    """
    Call GDAL translate utility to clip out tile by specified coords.

    Parameters:
        input_tif(str): a file to be clipped
        output_tif(str): a file where clip is saved
        projwin(str): dimensions for clipping the file
    """
    assert input_tif.endswith(".tif")

    ds = gdal.Open(input_tif)
    ds = gdal.Translate(output_tif, ds, projWin)
    # close the dataset
    ds = None

def clip_tile_by_dimensions(input_tile_path, reference_tile_path, output_tile_dir):
    """
    Clips a raster using the dimensions of another raster tile.

    Parameters:
    - input_tile_path (str): Path to the input raster tile to be clipped.
    - reference_tile_path (str): Path to the reference raster tile whose dimensions will be used for clipping.
    - output_tile_path (str): Path to save the clipped raster tile.
    """
    # Open the input and reference tiles
    input_ds = gdal.Open(input_tile_path)
    reference_ds = gdal.Open(reference_tile_path)

    if input_ds is None or reference_ds is None:
        print("Error: Unable to open raster tiles.")
        return

    # Get the bounding box of the reference tile
    reference_ulx, xres, xskew, reference_uly, yskew, yres = reference_ds.GetGeoTransform()
    reference_lrx = reference_ulx + (reference_ds.RasterXSize * xres)
    reference_lry = reference_uly + (reference_ds.RasterYSize * yres)

    # Set projwin for gdal.Translate using the dimensions of the reference tile
    projwin = [reference_ulx, reference_uly, reference_lrx, reference_lry]

    output_file_path = make_file_path(input_tile_path, output_tile_dir, sufix = ".tif")

    # Clip the input tile using the dimensions of the reference tile
    gdal.Translate(output_file_path, input_tile_path, projWin=projwin)

    # Close the datasets
    input_ds = None
    reference_ds = None
    return output_file_path


def make_file_path(input_file, output_file_dir, sufix):
    """
    Short pure function that takes file path and ouput dir and makes a path.

    Parameters:
        input_file(str): a file which will be copied
        output_file_dir(str): dir where new file 
        sufix(str)
    """
    output_file = f'{os.path.splitext(os.path.basename(input_file))[0]}{sufix}'
    output_file_path = os.path.join(output_file_dir, output_file)
    return output_file_path


def create_neighbour_vrt(file, output_folder):
    """
    Function takes a dem file, finds neighbouring tiles and then creates virtual raster.
    Virtual raster contains given tile and 8 surounding tiles (9 tiles total)

    Parameters:
        file(str): central file for creating vrt
        tmp_folder(str): temporal folder for storing this vrt
    
    Returns:
        vrt containing 9 tiles
    """
    os.makedirs(output_folder, exist_ok=True)
    file_name = os.path.basename(file)
    tile_folder = os.path.dirname(file)
    # find file coords
    # match something like "N23_00_W123"
    geocell_id_match = re.search(r'_(S|N)(\d+)_00_(W|E)(\d+)', file_name)
    if geocell_id_match:
        lat_direction = geocell_id_match.group(1)
        lat_value = int(geocell_id_match.group(2))
        lon_direction = geocell_id_match.group(3)
        lon_value = int(geocell_id_match.group(4))

        input_tiles = []
        # find neighbour files
        for i in range(lon_value - 1, lon_value + 2):
            for j in range(lat_value - 1, lat_value + 2):
                # Construct the file name for each surrounding tile
                tile_name = f'Copernicus_DSM_30_{lat_direction}{j:02d}_00_{lon_direction}{i:03d}_00_DEM.tif'
                tile_path = os.path.join(tile_folder, tile_name)
                if os.path.exists(tile_path):
                    input_tiles.append(tile_path)

        # create vrt with gdal
        vrt_path = os.path.join(output_folder, f'{file_name}.vrt')
        gdal_build_vrt(input_tiles, vrt_path)
        return vrt_path


def clip_layers_by_attribute(input_file, output_file, attribute_field='sea_level'):
    """
    Clip vector layers in the input file based on the 'sea_level' attribute and store the results in the output file.

    Parameters:
        input_file (str): The path to the input file containing vector layers.
        output_file (str): The path to the output file where clipped layers will be stored.
        attribute_field (str): The name of the attribute field used for categorization (default is 'sea_level').

    Returns:
        None

    Notes:
        This function reads vector layers from the input file, calculates the union
        of features based on the specified attribute, and then clips each feature
        with the feature having the next 'sea_level' value. The clipped features are
        stored in the output file.

        The function assumes that the input file contains multiple vector layers with
        the specified attribute for categorization.

    """
    # Open the input layer
    # driver = ogr.GetDriverByName("GPKG")
    driver = ogr.GetDriverByName("ESRI Shapefile")
    input_ds = ogr.Open(input_file, 0)  # 0 for read-only access

    # Create the output layer
    output_ds = driver.CreateDataSource(output_file)
    output_layer = output_ds.CreateLayer("clipped_layers", geom_type=ogr.wkbPolygon)

    # Loop through each pair of consecutive features and calculate the difference
    for i in range(1, input_ds.GetLayerCount()):
        layer1 = input_ds.GetLayerByIndex(i - 1)
        layer2 = input_ds.GetLayerByIndex(i)

        # Iterate through unique 'sea_level' values
        sea_levels = set(feature.GetField(attribute_field) for feature in layer1)
        for sea_level in sea_levels:
            # Get features with the current 'sea_level' in both layers
            features_layer1 = [feature for feature in layer1 if feature.GetField(attribute_field) == sea_level]
            features_layer2 = [feature for feature in layer2 if feature.GetField(attribute_field) == sea_level]

            # Create a union geometry for features in the first layer
            union_geom = ogr.Geometry(ogr.wkbPolygon)
            for feature in features_layer1:
                union_geom = union_geom.Union(feature.GetGeometryRef())

            # Calculate the difference for each feature in the second layer
            for feature2 in features_layer2:
                difference_geom = feature2.GetGeometryRef().Difference(union_geom)

                # Create a new feature in the output layer with the difference geometry
                output_feature = ogr.Feature(output_layer.GetLayerDefn())
                output_feature.SetGeometry(difference_geom)

                # Set the 'sea_level' attribute for the output feature
                output_feature.SetField(attribute_field, sea_level)

                # Add the feature to the output layer
                output_layer.CreateFeature(output_feature)

    # Clean up
    input_ds = None
    output_ds = None

def gdal_resample(input_tif:str, output_dir:str, resolution = 0.01, crs = 4326):
    """
    function for resampling/reprojecting dem file using gdal warp.
    """
    output_tif = os.path.join(output_dir, basename_withoutext(input_tif)+"_resampled.tif")
    warp_cmd = f"gdalwarp -overwrite -tr {resolution} {resolution} -r bilinear -of GTiff {input_tif} {output_tif}"
    os.system(warp_cmd)
    return 

# flows

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

