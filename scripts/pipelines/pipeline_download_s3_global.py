import requests
import os
import sys
import re
from bs4 import BeautifulSoup
import geopandas as gpd
from shapely.geometry import shape
from botocore.exceptions import NoCredentialsError
import subprocess
import xml.etree.ElementTree as ET
import lxml
import time
from statistics import mean
import datetime
from dotenv import load_dotenv
import psycopg2
import settings

from .pipeline_download_utils_soils import unzip_file, download_data, print_progress, move_tif_files_to_parent_directory, delete_non_matching_files_and_subdirectories
from .pipeline_transform_vrt_gdal import absolute_file_paths, filter_by_geocellid

# load configuration from settings.py
config = settings.get_config()
schema = settings.get_schema()
conn_parameters = settings.get_conn_parameters()
data_folder_path = settings.get_data_path()

# TODO: print progress wrapper

def list_s3_bucket_files(bucket_url:str):
    """
    Use AWS CLI and list files in an AWS S3 bucket under the specified URL.

    Args:
        bucket_url (str): The S3 bucket URL including the path to the source files.
    """
    try:
        files_list_command = f'aws s3 ls --no-sign-request {bucket_url}'
        process = subprocess.Popen(files_list_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        files_list, err = process.communicate()

        if process.returncode == 0:
            file_lines = files_list.split('\n')
            files = [line.split()[-1] for line in file_lines if line.strip()]
            return files
        else:
            print("Error:", err)
            return []
        
    except Exception as e:
        print("An error occurred:", str(e))
        return []


def download_s3_bucket_contents(bucket_url:str, output_directory:str):
    """
    Download files from an AWS S3 bucket to a local directory with progress percentage.
    Serves mostly for tiles found on https://registry.opendata.aws/esa-worldcover-vito/
    This function iterates over the list of file names and downloads each file individually.

    Args:
        bucket_url (str): The S3 bucket URL including the path to the source files.
        output_directory (str): The local directory where files will be downloaded.

    Note:
    - Make sure you have the AWS CLI installed and configured on your system.
    """
    files_list = list_s3_bucket_files(bucket_url)
    files_downloaded = 0
    total_files = len(files_list)

    for file in files_list:
        try:
            download_command = f'aws s3 cp --no-sign-request {bucket_url}{file} {output_directory}/{file}'
            print(f"Downloading {file} \n with command: {download_command}")
            process = subprocess.Popen(download_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

            progress_percentage = (files_downloaded / total_files) * 100
            print(f"Downloaded: {file}, Overall Progress: {progress_percentage:.2f}%")

        except Exception as e:
            print(f"An error occurred while downloading {file}: {str(e)}")

    return print("All files downloaded")


def extract_links_from_xml(url, url_tag="url", save_output = True, output_file = "urls_to_download.txt"):
    """ 
    Get links from a list of URLs obtained from a given URL of xml document using bs4.

    Args:
        url (str): The URL that contains the list of download links.
    """
    # Send a GET request to fetch the links
    response = requests.get(url)

    if response.status_code == 200:
        # Extract and split the response content into individual links
        soup = BeautifulSoup(response.content, 'xml')
        links = soup.find_all(url_tag)
        links = [link.text for link in links]

        if save_output == True:
            # links have to be string, not list, to be written into .txt file
            links_str = "\n".join(links)
            with open(output_file, 'w') as file:
                file.write(links_str)
        
        return links
    
    else:
        print(f"Failed to retrieve links. Status code: {response.status_code}")


def save_urls_to_file(urls, file_path):
    with open(file_path, 'w') as file:
        for url in urls:
            file.write(f"{url}\n")

def read_urls_from_file(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file.readlines()]


def download_tree_cover_density(urls_file = "TCD_urls_global_30m.txt", extent = "North America"):
    """
    Function iterates over range for North America, creates urls that hold Tree cover density (TCD) tiles to be downloaded and then downloads tiles.

    Parameters:
        urls_file(str): .txt file where urls are saved into
    """
    tree_cover_density_dir = config["tree_cover_density_dir"]
    tile_urls = []
    
    # Generate tile URLs within the bounding box of North America
    if extent == "North America":
        for lat in range(10, 90, 10):
            for lon in range(40, 180, 10):
                if lon < 100:
                    lon = f"0{lon}"
                tile_url = f"https://storage.googleapis.com/earthenginepartners-hansen/GFC2015/Hansen_GFC2015_treecover2000_{lat}N_{lon}W.tif"
                tile_urls.append(tile_url)
    
    else:
        raise ValueError("Invalid extent. Use 'North America'.")

    # Save the tile URLs to a text file
    urls_file = os.path.join(data_folder_path, urls_file)
    save_urls_to_file(tile_urls, urls_file)
    print(f"Urls for tree cover density tiles saved into {urls_file}")

    os.makedirs(tree_cover_density_dir, exist_ok=True)

    # Download tiles from the saved URLs
    total_count = len(tile_urls)
    for count, url in enumerate(tile_urls):
        # print_progress()
        file = download_data(url, tree_cover_density_dir, print_output=True)
        print(f"Downloaded, {count}/{total_count}: {file}")

    
    # total_count = len(tile_urls)
    # fun_times = []
    # for count, link in enumerate(tile_urls):
    #     download_args = [link, tree_cover_density_dir]
    #     result, fun_times = print_progress(download_data, download_args, count, total_count, fun_times)
    
    return print(f"Downloading task done and saved into {tree_cover_density_dir}")

    
# flows here
def download_ESA_worldcover_flow():
    """
    Download and unzip ESA World Cover data from an S3 bucket.

    This function downloads ESA World Cover data files from a specified S3 bucket
    and unzips them to an output directory. It then prints a completion message.
    """
    bucket_url = config["esa_world_cover_url"]
    output_directory = config["esa_world_cover_dir"]
    os.makedirs(output_directory, exist_ok=True)
    download_s3_bucket_contents(bucket_url, output_directory)
    file_paths = absolute_file_paths(output_directory, '.zip')
    for i in file_paths:
        unzip_file(i, output_directory)
    print("All files downloaded and unzipped")
    # TODO: remove the zip files
    return


def download_ESA_global_dem(resolution = "90", files_filter = True, filter_list = []):
    """
    Download ESA Global Digital Elevation Model data for a specified resolution, either 90m or 30m.

    This function downloads data from ESA's Global Digital Elevation Model (DEM)
    dataset for a specified resolution. It retrieves URLs from an XML file and
    downloads the data files to a specified directory. URLs are saved to .txt file.
    Source of data: https://spacedata.copernicus.eu/collections/copernicus-digital-elevation-model

    Args:
    - resolution (str): The resolution of the DEM data (default is "90").
    - files_filter (bool): To apply or not to apply filtering on files to be downloaded.
    - filter_list (list): Filter which is a list to be iterated over, geospatial extent by which files are downloaded.
    """
    esa_global_dem_xml = config[f"esa_global_dem_{resolution}_xml"]
    esa_global_dem_dir = config[f"esa_global_dem_{resolution}_dir"]
    esa_global_dem_urls = config[f"esa_global_dem_{resolution}_urls"]

    # creating directory where the downloaded data will be saved
    os.makedirs(esa_global_dem_dir, exist_ok=True)

    # if txt file that where urls are written already exists, then it's not created again
    if not os.path.exists(esa_global_dem_urls):
        links = extract_links_from_xml(esa_global_dem_xml, "nativeDemUrl", output_file = esa_global_dem_urls)
    else:
        with open(esa_global_dem_urls, 'r') as file:
            links = file.readlines()
            links  = [line.strip() for line in links]

    # filter files to be downloaded
    if files_filter == True:
        filtered_links = filter_by_geocellid(links)

        # finally, change the links variable into filtered_links
        links = filtered_links

    # data downloading part
    total_count = len(links)
    fun_times = []
    for count, link in enumerate(links):
        download_args = [link, esa_global_dem_dir]
        result, fun_times = print_progress(download_data, download_args, count, total_count, fun_times)

    return print("Downloading task is done")


def download_ESA_global_dem_90_flow():
    download_ESA_global_dem("90", files_filter = False)
    esa_global_dem_90_dir = config[f"esa_global_dem_90_dir"]
    # unzip
    file_paths = absolute_file_paths(esa_global_dem_90_dir, '.tar')
    total_count = len(file_paths)
    fun_times = []
    for count, i in file_paths:
        fun_args = [i, esa_global_dem_90_dir]
        print_progress(unzip_file, fun_args, count, total_count, fun_times)

    move_tif_files_to_parent_directory(esa_global_dem_90_dir)
    delete_non_matching_files_and_subdirectories(esa_global_dem_90_dir)
    return print("task for ESA dem 90m downloading and unzipping done!")


def download_ESA_global_dem_30_flow():
    download_ESA_global_dem("30", files_filter = True)
    esa_global_dem_30_dir = config[f"esa_global_dem_30_dir"]
    # unzip
    file_paths = absolute_file_paths(esa_global_dem_30_dir, '.tar')
    file_paths = list(file_paths)
    total_count = len(file_paths)
    fun_times = []
    for count, i in enumerate(file_paths):
        fun_args = [i, esa_global_dem_30_dir]
        print_progress(unzip_file, fun_args, count, total_count, fun_times)

    delete_non_matching_files_and_subdirectories(esa_global_dem_30_dir)
    return print("task for ESA dem 30m downloading and unzipping done!")
