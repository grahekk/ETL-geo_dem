import os
from owslib.wms import WebMapService
import fiona
import rasterio
from multiprocessing import Pool, Manager
import re
import yaml
import settings

config = settings.get_configuration()

def rename_files_in_directory(directory):
    # Change to the specified directory
    os.chdir(directory)

    # Define the pattern to match the filenames
    pattern = re.compile(r'cdem\.color-shaded-relief_tile_(-?\d+)_(\d+)\.tiff')

    # Loop through the files in the directory
    for filename in os.listdir(directory):
        if os.path.isfile(filename):
            match = pattern.match(filename)
            if match:
                # Extract the parts you want from the filename using regular expressions
                new_name = f'cdem_-{match.group(1)}_{match.group(2)}.tiff'
                # Rename the file
                os.rename(filename, new_name)


def tile_contains_only_nodata(tile_path):
    with rasterio.open(tile_path, 'r') as src:
        data = src.read()
        nodata_value = src.nodata

        # Check if all values in the tile are equal to the nodata value
        return (data == nodata_value).all()


def delete_empty_nodata_tiles(tile_directory):
    for root, _, files in os.walk(tile_directory):
        for filename in files:
            tile_path = os.path.join(root, filename)
            if tile_contains_only_nodata(tile_path):
                os.remove(tile_path)
                print(f"Deleted empty nodata tile: {tile_path}")


def create_WMS_tile(wms_url:str, layer_name:str, tile_bbox, output_directory, crs, total_tiles, counter, lock, timeout=2500):
    """
    Download a Digital Elevation Model (DEM) tile from a Web Map Service (WMS) and save it as a GeoTIFF file.

    Parameters:
        wms_url (str): The URL of the WMS service.
        layer_name (str): The name of the DEM layer on the WMS service.
        tile_bbox (tuple): A tuple representing the bounding box (extent) of the tile as (min_lon, min_lat, max_lon, max_lat).
        output_directory (str): The directory where the downloaded tile will be saved.
        crs (str): The Coordinate Reference System (CRS) to use for the tile.
        total_tiles (int): The total number of tiles to be downloaded.
        counter (multiprocessing.Value): A multiprocessing shared counter for tracking downloaded tiles.
        lock (multiprocessing.Lock): A multiprocessing lock for synchronizing counter updates.
        timeout (int, optional): The timeout in milliseconds for the WMS request (default is 2500).

    Note:
        timeout is set high because a lot of tiles need some time for response
    """
    tile_name = f"cdem_{tile_bbox[0]}_{tile_bbox[1]}.tif"
    output_file = os.path.join(output_directory, tile_name)
    try:
        wms = WebMapService(wms_url, version="1.3.0", timeout=timeout)
        # Create a request for the WMS service for the tile
        response = wms.getmap(layers=[layer_name],
                            styles=["default"],
                            srs=crs,
                            bbox=tile_bbox,
                            size=(4000, 4000),  
                            format="image/tiff",
                            transparent=True,
                            version="1.3.0")

        # Save the downloaded GeoTIFF data to a file
        with open(output_file, "wb") as f:
            f.write(response.read())
        print(f"Tile {tile_name} downloaded and saved.")

        # Update progress
        with lock:
            counter.value += 1
            tiles_downloaded = counter.value
            progress = (tiles_downloaded / total_tiles) * 100
            print(f"Downloaded ({tiles_downloaded}/{total_tiles}) Progress: {progress:.2f}%")

    except Exception as e:
            print(f"Error downloading tile {tile_name}: {str(e)}")


def download_WMS_tiles(wms_url:str, layer_name:str, bounding_box, output_directory:str, crs="EPSG:4269", tile_size=1):
    """
    Download tiles from a Web Map Service (WMS) and save them to a directory.

    Parameters:
        wms_url (str): The URL of the WMS service.
        layer_name (str): The name of the layer to download from the WMS service.
        bounding_box (shapely.geometry.box): The bounding box (extent) specifying the area of interest.
        output_directory (str): The directory where downloaded tiles will be saved.
        crs (str, optional): The Coordinate Reference System (CRS) to use (default is "EPSG:4269").
        tile_size (int, optional): The size of each tile in degrees (default is 1 degree).

    Description:
        This function connects to a WMS service and downloads tiles of a specified layer within a given
        bounding box. It splits the area into a grid of tiles and downloads them concurrently using a
        multiprocessing pool. The tiles are saved to the specified output directory.

    Note:
        - The function uses the 'WebMapService' class from the 'owslib' library, so make sure to install
          the 'owslib' library before using this function.
        - The 'bounding_box' parameter should be a Shapely 'box' object that defines the extent of the area
          you want to download tiles for.
        - Commented part of the code is for loop which was replaced by multiprocessing part below it.
    """
    # Create the output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)

    # Connect to the WMS service
    wms = WebMapService(wms_url, version="1.3.0")
    print('available layers from url: ')
    print(list(wms.contents))

    # Get the bounding box (extent) of the WMS layer
    bbox = bounding_box.bounds
    min_lon, min_lat, max_lon, max_lat = bbox

    # Calculate the grid of 1x1 degree tiles within the bounding box
    lon_step = tile_size
    lat_step = tile_size
    lon_range = range(int(min_lon), int(max_lon), lon_step)
    lat_range = range(int(min_lat), int(max_lat), lat_step)

    total_tiles = len(lon_range) * len(lat_range)

    # Iterate over the grid and create/download each tile
    # tiles_downloaded = 0
    # for lon in lon_range:
    #     for lat in lat_range:
    #         # Define the bounding box for the tile
    #         tile_bbox = (lon, lat, lon + lon_step, lat + lat_step)

    #         # Create and download the tile
    #         create_dem_tile(wms, layer_name, tile_bbox, output_directory)

    #         tiles_downloaded += 1
    #         progress = (tiles_downloaded / total_tiles) * 100
    #         print(f"Downloaded ({tiles_downloaded}/{total_tiles}) Progress: {progress:.2f}%")

    tile_args = []
    for lon in lon_range:
        for lat in lat_range:
            # Define the bounding box for the tile
            tile_bbox = (lon, lat, lon + lon_step, lat + lat_step)
            tile_args.append((wms_url, layer_name, tile_bbox, output_directory, crs, total_tiles))

    manager = Manager()
    counter = manager.Value('i', 0)
    lock = manager.Lock()

    with Pool(processes=12) as pool: 
        pool.starmap(create_WMS_tile, [(arg[0], arg[1], arg[2], arg[3], arg[4], arg[5], counter, lock) for arg in tile_args])


# Downlaod WMS flow:
def download_canada_WMS_flow():
    """
    Utilizes `download_WMS_tiles()` function and downloads tiles iteratively by accessing WMS.
    Multiprocessing is used for this purpose and timeout is set to 2500 since some tiles need more time to response.
    Configuration variables `wms_url`, `layer_name`, `canada_territory`, `output_directory` are listed below and taken from config file.
    Crs is `EPSG:4269`.
    """
    wms_url = config["canada_dem_WMS_url"]
    layer_name = config["canada_dem_wms_layer_name"]
    output_directory = config["canada_dem_WMS_dir"]

    # output_directory_tiles = f'{output_directory}/canada_tiles_boundaries'
    # canada_grid_url = config["canada_grid_url"]

    crs = "EPSG:4269"
    canada_territory = fiona.open(config["canada_territory"])
    download_WMS_tiles(wms_url, layer_name, canada_territory, output_directory, crs)
    print("Canadian DEM downloaded from WMS!")


