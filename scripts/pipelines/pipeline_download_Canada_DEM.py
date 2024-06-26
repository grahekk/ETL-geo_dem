import os
import requests
from bs4 import BeautifulSoup
from urllib.parse import urljoin

from .pipeline_download_utils_soils import unzip_file, delete_non_tif_files, move_files_from_subfolder_to_parent
from .pipeline_transform_vrt_gdal import absolute_file_paths

import settings
config = settings.get_config()



def download_file(file_url, destination_directory, count, total_count):
    """
    Download a single file from a given URL.
    
    Args:
        file_url (str): The URL of the file to be downloaded.
        destination_directory (str): The local directory where the file will be saved.
        count (int): Current file count.
        total_count (int): Total number of files to download.
    """
    filename = os.path.basename(file_url)
    local_file_path = os.path.join(destination_directory, filename)

    file_response = requests.get(file_url)
    if file_response.status_code == 200:
        # Save the file to the destination directory
        with open(local_file_path, "wb") as f:
            f.write(file_response.content)
        print(f"Downloaded {local_file_path} ({count}/{total_count}) - {count/total_count:.2%} complete")
    else:
        print(f"Failed to download: {file_url}")


def download_files_from_directory_recursive(directory_url, destination_directory):
    """
    Recursively download files from a directory URL and its subdirectories
    to a destination directory.
    
    Args:
        directory_url (str): The URL of the directory containing the files to be downloaded.
        destination_directory (str): The local directory where the files will be saved.
    """
    # Create the destination directory if it doesn't exist
    os.makedirs(destination_directory, exist_ok=True)

    # Send an HTTP GET request to the directory URL
    response = requests.get(directory_url)

    # Check if the request was successful
    if response.status_code == 200:
        # Parse the HTML content of the page using BeautifulSoup
        soup = BeautifulSoup(response.content, "html.parser")

        # Find all links (href attributes) in the page
        links = soup.find_all("a")
        links_digit = []
        for link in links:
            subdirectory_url = link.get("href")
            # Check if the link starts with a number (exclude others)
            if subdirectory_url and subdirectory_url[0].isdigit():
                full_url = urljoin(directory_url, subdirectory_url)
                links_digit.append(full_url)

        # finding links that contain tif files
        download_content = []
        for link in links_digit:
            response = requests.get(link)
            soup = BeautifulSoup(response.content, "html.parser")
            page_content = soup.find_all("a")

            for d in page_content:
                file_url = d.get("href")
                if 'tif' in file_url:
                    full_file_url = urljoin(link, file_url)
                    download_content.append(full_file_url)

        # Iterate through the links and download files 
        for count, file in enumerate(download_content):
            filename = os.path.basename(file)
            local_file_path = os.path.join(destination_directory, filename)

            file_response = requests.get(file)
            if file_response.status_code == 200:
                # Save the file to the destination directory
                with open(local_file_path, "wb") as f:
                    f.write(file_response.content)
                # print(f"Downloaded: {local_file_path}")
                print(f"Downloaded {local_file_path} ({count}/{len(download_content)}) - {count/len(download_content):.2%} complete")
            else:
                print(f"Failed to download: {full_file_url}")

        # with Pool(processes=12) as pool:  # Adjust the number of processes as needed
        #     for count, file_url in enumerate(download_content, start=1):
        #         pool.apply_async(download_file, args=(file_url, destination_directory, count, len(download_content)))

    else:
        print(f"Failed to access: {directory_url}")

# all the functions packed into flow
def download_canadian_dem_flow():
    canada_directory_url = config['canada_dem_url']
    canada_destination_directory_zip = config['canada_dem_dir_zip']
    canada_destination_directory_unzip = config['canada_dem_dir_unzip']

    # first download the Canadian dem tiles
    download_files_from_directory_recursive(canada_directory_url, canada_destination_directory_zip)

    # and then unzip these files
    file_paths = absolute_file_paths(canada_destination_directory_zip, '.zip')
    for i in file_paths:
        unzip_file(i, canada_destination_directory_unzip)

    # then delete these tif files and put them to parent directory
    delete_non_tif_files(canada_destination_directory_unzip, '.tif')
    move_files_from_subfolder_to_parent(canada_destination_directory_unzip)
