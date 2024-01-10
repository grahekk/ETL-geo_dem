import requests
import zipfile
import tarfile
import os
import sys
import shutil
import time
import datetime
from statistics import mean

import settings
config = settings.get_config()


def print_progress(subject_function:callable, function_args:list, count:int, total_count:int, fun_times):
    """
    Execute a subject function while printing progress information.

    This function executes a subject function with specified arguments and
    tracks the progress of the function's execution. It prints progress
    information including completion percentage and estimated time remaining.

    Args:
    - subject_function (callable): The function to execute and track.
    - function_args (list): List of arguments to pass to the subject function.
    - count (int): Current count or index of the task being tracked.
    - total_count (int): Total number of tasks to complete.
    - fun_times (list): A list to store execution times of the subject function.

    Returns:
    - result: The result of executing the subject function.
    - fun_times (list): Updated list of execution times.
    """
    start_time = time.time()
    result = subject_function(*function_args)
    fun_time = time.time() - start_time 
    # results.append((lat, lon, city, country, result, fun_time))
    fun_times.append(fun_time)

    # pretty print out for task completeness
    if count % 10 == 0:
        mean_time = mean(fun_times)
        global time_left
        time_left = round(((total_count-count)*mean_time)/60,2)

    percentage = round(count / total_count * 100, 2)
    sys.stdout.write('\r')
    sys.stdout.write(f'task for {subject_function.__name__} completeness: {percentage:.2f}% time left (d:hh:mm.ss): {datetime.timedelta(seconds=time_left)}')
    sys.stdout.flush()
    return result, fun_times


def download_data(url:str, output_dir, print_output = False, skip_existing_files = True):
    """
    Downloads the data using request library from given url into given filename. 
    
    Args:
        url (str): The URL of the file to be downloaded.
        filename (str): The local path where the file will be saved.
    """
    filename = os.path.join(output_dir, os.path.basename(url))
    
    #skip if file exists
    if skip_existing_files==True:
        if os.path.exists(filename):
            print(f"{filename} exists, skipping download")
            return filename

    try:
        # Send an HTTP GET request to the URL
        response = requests.get(url)
        # Check if the request was successful (status code 200)
        if response.status_code == 200:
            # Open the file in binary write mode and write the content of the response
            with open(filename, "wb") as file:
                file.write(response.content)
                
            if print_output == True:
                print(f"File '{filename}' downloaded successfully.")
        else:
            print(f"Failed to download file. Status code: {response.status_code}")
    except Exception as e:
        print(f"An error occurred: {str(e)}")
    return filename


def unzip_file(zip_filename:str, output_dir:str):
    """
    Unzips a single file from a given path into given directory. 
    If file ends with ".tar" it utilizes `tarfile` module.
    
    Args:
        zip_filename (str): The path of the file to be unzipped.
        output_dir (str): The local directory where the file will be saved.
    """
    os.makedirs(output_dir, exist_ok=True)
    if zip_filename.endswith('.tar'):
        try:
            with tarfile.open(zip_filename) as file_ref:
                file_ref.extractall(output_dir)
            # print(f"File '{zip_filename}' successfully unzipped to '{output_dir}'.")
        except Exception as e:
            print(f"An error occurred while unpacking(un-taring) '{zip_filename}': {str(e)}")
    else:
        try:
            with zipfile.ZipFile(zip_filename, 'r') as zip_ref:
                zip_ref.extractall(output_dir)
            # print(f"File '{zip_filename}' successfully unzipped to '{output_dir}'.")
        except Exception as e:
            print(f"An error occurred while unzipping '{zip_filename}': {str(e)}")


def delete_non_tif_files(directory:str, extension:str):
    """
    Deletes all files in the given directory that do not end with .tif or some other specified extension.

    Args:
        directory (str): Path to the directory.
        extension (str): Wanted extension to be left in a directory
    """
    for dirpath, _, filenames in os.walk(directory):
        for filename in filenames:
            if not filename.endswith(extension):
                file_path = os.path.join(dirpath, filename)
                os.remove(file_path)
                print(f"Deleted: {file_path}")


def move_files_from_subfolder_to_parent(subfolder_path):
    """
    Moves all files from a subfolder to the directory above the subfolder.
    Or move files from subdirectories in directory into directory.

    Args:
        subfolder_path (str): Path to the subfolder whose files need to be moved.

    Returns:
        None
    """
    # Get the path of the parent directory
    parent_directory = os.path.dirname(os.path.dirname(subfolder_path))

    # List all files in the subfolder
    files_to_move = [f for f in os.listdir(subfolder_path) if os.path.isfile(os.path.join(subfolder_path, f))]

    # Move the files to the parent directory
    for file_to_move in files_to_move:
        source_path = os.path.join(subfolder_path, file_to_move)
        destination_path = os.path.join(parent_directory, file_to_move)
        shutil.move(source_path, destination_path)


def delete_non_matching_files_and_subdirectories(directory:str):
    """
    Recursively deletes files and subdirectories in the specified directory and its
    subdirectories that do not end with ".tar" or "DEM.tif".

    Parameters:
    - directory (str): The path to the directory from wh ich non-matching files and
      subdirectories should be deleted
    """
    for root, dirs, files in os.walk(directory, topdown=False):
        for file in files:
            file_path = os.path.join(root, file)
            if not (file.endswith(".tar") or file.endswith("DEM.tif")):
                os.remove(file_path)
                print(f"Deleted file: {file_path}")

        for dir_name in dirs:
            dir_path = os.path.join(root, dir_name)
            try:
                os.rmdir(dir_path)
                print(f"Deleted directory: {dir_path}")
            except OSError:
                pass

def move_tif_files_to_parent_directory(root_dir:str):
    """
    Take all the tif files in subdirectories into one directory above.

    Parameters:
        root_dir(str): directory where subdirectories with tif files are
    """
    # Iterate through the root directory and its subdirectories
    files_to_move = [f for f in os.listdir(root_dir) if os.path.isfile(os.path.join(root_dir, f))]
    total_count = len([file.endswith(".tif") for file in files_to_move])

    for root, _, files in os.walk(root_dir):
        count = 0
        fun_times = []
        for file in files:
            if file.endswith(".tif"):
                # Construct the source file path
                source_path = os.path.join(root, file)

                # Construct the destination file path in the parent directory
                destination_path = os.path.join(root_dir, file)

                # Check if the file is not already in the parent directory
                if source_path != destination_path:
                    # Move the file to the parent directory
                    # shutil.move(source_path, destination_path)
                    # print(f"Moved '{file}' to '{root_dir}'")
                    fun_args = [source_path, destination_path]
                    print_progress(shutil.move, fun_args, count, total_count, fun_times)
                count = count+1


def save_geojson_from_url(url:str, file_dir:str, file_path:str):
    """
    Download GeoJSON data from a URL and save it to a file.

    Args:
        url (str): The URL of the GeoJSON data.
        file_path (str): The file path where the GeoJSON data will be saved.
    """
    response = requests.get(url)
    os.makedirs(file_dir, exist_ok=True)

    if response.status_code == 200:
        # Get the GeoJSON data from the response
        geojson_data = response.text

        # Write the GeoJSON data to the file
        with open(file_path, "w") as file:
            file.write(geojson_data)

        print(f"GeoJSON data saved to {file_path}")
    else:
        print(f"Failed to retrieve data. HTTP status code: {response.status_code}")


def download_flow(data_name:str):
    """
    Download and process data for a specified data_name.

    This function downloads data from the URL specified in the provided configuration
    dictionary, stores it in the directory specified, and then unzips the downloaded
    file. It also prints a message indicating the completion of the download process.

    Args:
        data_name (str): A unique identifier for the data you want to download.
        config (dict): A dictionary containing configuration parameters, including the
            URL and directory information needed for downloading and processing the data.

    Returns:
        str: The filename of the downloaded and unzipped data.
    """
    url = config[f"{data_name}_url"]
    dir = config[f"{data_name}_dir"]
    os.makedirs(dir, exist_ok=True)
    filename = download_data(url, dir)
    unzip_file(filename, dir)
    print(f"{__name__} for {data_name} done")
    return filename


# flows here

def download_us_soil_flow():
    """
    Download and process data for US soil.

    Calls the `download_flow` function to download and process data for US soil
    based on the configuration specified in `config`.
    """
    return download_flow("US_soil")

def download_cec_lc_flow():
    """
    Download and process data from CEC for Land Cover in North America.

    Calls the `download_flow` function to download and process data for CEC Land Cover
    in North America based on the configuration specified in `config`.
    """
    return download_flow("cec_lc_north_america")

def download_canada_soil_flow():
    """
    Download and process data for Canadian soil.

    Calls the `download_flow` function to download and process data for Canadian soil
    based on the configuration specified in `config`.
    """
    return download_flow("canada_soil")

def download_world_countries_boundaries_flow():
    """
    Download and process data for boundaries of all the countries in the world.

    Calls the `download_flow` function to download and process data for Canadian soil
    based on the configuration specified in `config`.
    """
    return download_flow("world_countries_boundaries")

def download_world_continent_boundaries_flow():
    """
    Download geojson data from ESRI's api. The file holds info of continent boundaries.
    """
    url = config["world_continents_boundaries_url"]
    file_dir = config["world_continents_boundaries_dir"]
    file_path = config["world_continents_boundaries_file"]
    return save_geojson_from_url(url, file_dir, file_path)

def download_ESA_global_dem_grid_flow():
    """
    download grid file from given url.
    """
    return download_flow("esa_global_dem_grid")