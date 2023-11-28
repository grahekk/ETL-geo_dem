import os
import requests


def download_files_from_urls(urls:str, output_folder:str):
    """
    Downloads files from the given URLs and saves them to the specified output folder.

    Args:
        urls (str or list of str): The URL or list of URLs to download files from.
        output_folder (str): The directory where downloaded files will be saved.

    Notes:
        - If `output_folder` doesn't exist, it will be created.
        - Progress is printed out after downloading each file, showing the completion status.
    """
    # Create the output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    total_urls = len(urls)

    # Iterate through the URLs and download the files
    for i, url in enumerate(urls, start=1):

        filename = os.path.basename(url)
        local_filepath = os.path.join(output_folder, filename)

        # Download the file
        response = requests.get(url)

        # Check if the request was successful (HTTP status code 200)
        if response.status_code == 200:
            with open(local_filepath, 'wb') as f:
                f.write(response.content)
            print(f"Downloaded {filename} ({i}/{total_urls}) - {i/total_urls:.2%} complete")
        else:
            print(f"Failed to download {filename} - HTTP status code {response.status_code}")

    print("Download completed.")


def get_url_from_file(data_urls:str):
    """
    Reads URLs from a text file and returns them as a list of strings.

    Args:
        data_urls (str): The path to the text file containing the URLs.

    Returns:
        list of str: A list of URLs read from the file.

    Notes:
        - The provided text file should contain one URL per line.
        - The function removes any trailing newline characters from the URLs.
    """
    with open(data_urls, 'r') as urls_file:
        urls = urls_file.readlines()
    urls = [line[:-1] for line in urls if line[-1]=='\n']
    return urls


# Usage:
USGS_DEM_FOLDER = "/mnt/volume-nbg1-1/shared/nikola/USGS_DEM_10m"
USGS_DEM_URLS = "/home/nikola/4_north_america/USGS_dem_urls.txt"
NASA_GLOBAL_DEM_FOLDER = "/mnt/volume-nbg1-1/shared/nikola/NASA_DEM"
NASA_GLOBAL_DEM_URLS = "/home/nikola/4_north_america/data/NASA_global_dem_urls.txt"

output_folder = NASA_GLOBAL_DEM_FOLDER
data_urls = NASA_GLOBAL_DEM_URLS

def download_us_dem_flow():
    urls = get_url_from_file(data_urls)
    download_files_from_urls(urls, output_folder)