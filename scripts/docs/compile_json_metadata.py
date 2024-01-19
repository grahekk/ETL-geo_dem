import os
import sys
import json
import subprocess
import yaml
from pathlib import Path
import datetime
import re
import ast

from scripts import settings

# load configuration
config = settings.get_config()
metadata_description_path = "/home/nikola/4_north_america/GeoDataPump/scripts/docs/metadata_description.yaml"
with open(metadata_description_path, 'r') as file:
    metadata_description = yaml.safe_load(file)


def update_file_size(file_path:str):
    """
    Update the file_size attribute with the human-readable file size.
    """
    # Get the file size in bytes
    file_size_bytes = os.path.getsize(file_path)

    # Convert file size to a human-readable format (e.g., KB, MB, GB)
    suffixes = ['B', 'KB', 'MB', 'GB', 'TB']
    index = 0
    while file_size_bytes > 1024 and index < len(suffixes) - 1:
        file_size_bytes /= 1024.0
        index += 1

    return "{:.2f} {}".format(file_size_bytes, suffixes[index]) 

def get_directory_size(directory_path):
    """
    Get the total size of a directory and return it in a human-readable format.
    """
    total_size_bytes = 0

    for dirpath, dirnames, filenames in os.walk(directory_path):
        for filename in filenames:
            file_path = os.path.join(dirpath, filename)
            total_size_bytes += os.path.getsize(file_path)

    # Convert total size to a human-readable format (e.g., KB, MB, GB)
    suffixes = ['B', 'KB', 'MB', 'GB', 'TB']
    index = 0
    while total_size_bytes > 1024 and index < len(suffixes) - 1:
        total_size_bytes /= 1024.0
        index += 1

    return "{:.2f} {}".format(total_size_bytes, suffixes[index])

def get_resolution(file_path):
    try:
        # Use gdalinfo to get raster resolution
        gdal_info_output = subprocess.check_output(["gdalinfo", file_path], text=True)
        for line in gdal_info_output.split("\n"):
            if "Pixel Size" in line:
                resolution_parts = line.split(" ")
                pixel_size = resolution_parts[3]
                resolution_tuple = ast.literal_eval(pixel_size)
                first_num = round(resolution_tuple[0]/0.00001, 3)
                second_num = round(resolution_tuple[1]/0.00001, 3)
                resolution = f"Pixel size in degrees is {pixel_size} and in meters it is ({first_num, second_num})"
                return resolution
    except subprocess.CalledProcessError:
        print(f"Error: Unable to get resolution for {file_path}")
    return None, None

def get_gdal_crs(file_path):
    try:
        # Use gdalinfo to get raster resolution
        gdal_info_output = subprocess.check_output(["gdalinfo", file_path], text=True)
        for line in gdal_info_output.split("\n"):
            if 'ID["EPSG"' in line:
                # get epsg numbers
                crs = re.findall(r'\d+', line.strip())
                # make epsg string
                crs = f'EPSG: {crs[0]}'
                return crs
    except subprocess.CalledProcessError:
        print(f"Error: Unable to get resolution for {file_path}")
    return None, None

def ogr_get_coordinate_system(file_path):
    try:
        # Use ogrinfo to get coordinate system
        ogr_info_output = subprocess.check_output(["ogrinfo", file_path, "-al", "-so"], text=True)
        for line in ogr_info_output.split("\n"):
            if "Layer SRS WKT" in line:
                return line.strip().replace("Layer SRS WKT:", "")
    except subprocess.CalledProcessError:
        print(f"Error: Unable to get coordinate system for {file_path}")
    return None


class FileMetadata:
    def __init__(self, file_path):
        """
        Initialize a FileMetadata object.

        Parameters:
            file_path (str): The path to the file.

        Attributes:
            file_path (str): The path to the file.
            date_created (datetime): The date and time when the file was created.
            date_modified (datetime): The date and time when the file was last modified.
            resolution (tuple): The resolution of the file (e.g., (width, height)).
            coordinate_system (str): The coordinate system of the file.
            file_size (str): The human-readable file size (e.g., "2.45 MB").
            description (str): A description of the file.
        """
        self.file_path = file_path
        self.file_size = None
        self.date_created = None
        self.date_modified = None
        self.resolution = None
        self.coordinate_system = None
        # from config file
        self.description = None
        self.time_to_make = None
        self.download_url = None
        self.creation_method = None
        self.script_path = None
        self.version = "1.0"
        self.data_attributes = None
        self.image_path = None

    def update_metadata(self):
        """
        Update metadata attributes such as date_created, date_modified, resolution, and coordinate_system.
        """
        # Get the date created (in seconds since the epoch)
        created_timestamp = os.path.getctime(self.file_path)
        self.date_created = datetime.datetime.fromtimestamp(created_timestamp)

        # Get the date modified (in seconds since the epoch)
        modified_timestamp = os.path.getmtime(self.file_path)
        self.date_modified = datetime.datetime.fromtimestamp(modified_timestamp)
        
        # if a file is vrt, it is pointing to a dataset, which is why wee need size of a directory
        if self.file_path.endswith('.vrt'):
            self.file_size = get_directory_size(self.file_path.replace(".vrt",""))
        else:
            self.file_size = update_file_size(self.file_path)

        self.resolution = get_resolution(self.file_path)

        if self.file_path.endswith(".tif") or self.file_path.endswith(".vrt"):
            self.coordinate_system = get_gdal_crs(self.file_path) 
        else:
            self.coordinate_system = ogr_get_coordinate_system(self.file_path)  


def create_metadata_list(file_paths, descriptions, script_paths, creation_methods, times_to_make, download_urls, versions, data_attributes, image_paths):
    """
    Create a list of FileMetadata objects based on provided file_paths and descriptions.

    Parameters:
        file_paths (list): List of file paths.
        descriptions (list): List of file descriptions.

    Returns:
        list: List of FileMetadata objects.
    """
    metadata_list = []

    for file_path, description, script_path, creation_method, time_to_make, download_url, version, data_attribute, image_path in zip(file_paths, descriptions, script_paths, creation_methods, times_to_make, download_urls, versions, data_attributes, image_paths):
        metadata_obj = FileMetadata(file_path)
        metadata_obj.update_metadata()
        metadata_obj.description = description
        metadata_obj.time_to_make = time_to_make
        metadata_obj.download_url = download_url
        metadata_obj.creation_method = creation_method
        metadata_obj.script_path = script_path
        metadata_obj.version = version
        metadata_obj.data_attributes = data_attribute
        metadata_obj.image_path = image_path
        # append the object to metadata list
        metadata_list.append(metadata_obj)

    return metadata_list

def export_metadata_to_json(metadata_list, json_file_path):
    """
    Export metadata attributes of FileMetadata objects to a JSON file.

    Parameters:
        metadata_list (list): List of FileMetadata objects.
        json_file_path (str): Path to the JSON file.
    """
    serialized_data = []

    for metadata_obj in metadata_list:
        serialized_data.append({
            "file_path": metadata_obj.file_path,
            "date_created": metadata_obj.date_created.strftime("%Y-%m-%d %H:%M:%S"),
            "date_modified": metadata_obj.date_modified.strftime("%Y-%m-%d %H:%M:%S"),
            "resolution": metadata_obj.resolution,
            "coordinate_system": metadata_obj.coordinate_system,
            "file_size": metadata_obj.file_size,
            "version": metadata_obj.version,
            "download_url": metadata_obj.download_url,
            "creation_method": metadata_obj.creation_method,
            "script_path": metadata_obj.script_path,
            "time_to_make": metadata_obj.time_to_make,
            "attributes": metadata_obj.data_attributes,
            "description": metadata_obj.description,
            "image_path": metadata_obj.image_path,
        })

    with open(json_file_path, "w") as json_file:
        json.dump(serialized_data, json_file, indent=2)


def main():
    file_paths = [config["esa_world_cover_vrt"], 
                  config["NA_slope"],
                  config["NA_aspect"],
                  config["NA_geomorphon"],
                  config["NA_tree_cover_density_tif"],
                  config["usgs_dem_vrt"],
                  config["esa_global_dem_30_vrt"],
                  config["esa_global_dem_90_vrt"]]

    download_urls = [config["esa_world_cover_url"], 
                    "Transformed from esa global dem 90m",
                    "Transformed from esa global dem 90m",
                    "Transformed from esa global dem 90m",
                    config["NA_tree_cover_density_url"],
                    config["usgs_dem_url"],
                    config["esa_global_dem_30_xml"],
                    config["esa_global_dem_90_xml"]]
    
    keys = ["esa_world_cover", "NA_slope", "NA_aspect", "NA_geomorphon", "NA_tree_cover_density", "USGS_dem", "ESA_global_dem_30m", "ESA_global_dem_90m"]

    descriptions = [metadata_description[key]["description_hrv"] for key in keys]
    script_paths = [metadata_description[key]["script_path"] for key in keys]
    creation_methods = [metadata_description[key]["creation_method"] for key in keys]
    times_to_make = [metadata_description[key]["time_to_make"] for key in keys]
    versions = [metadata_description[key]["version"] for key in keys]
    data_attributes = [metadata_description[key]["attributes"] for key in keys]
    image_paths = [metadata_description[key]["image_path"] for key in keys]

    metadata_list = create_metadata_list(file_paths, descriptions, script_paths, creation_methods, times_to_make, download_urls, versions, data_attributes, image_paths)
    export_metadata_to_json(metadata_list, config["metadata_list"])

if __name__ == "__main__":
    main()
