import rasterio
import subprocess
import re
import numpy as np
import xml.etree.ElementTree as ET
import os
import sys
import pandas as pd
import time
import math
import glob
from statistics import mean
import rasterio
from rasterio.transform import from_origin
from rasterio.warp import transform
from rasterio.windows import Window


sys.path.append('/usr/lib/python3/dist-packages')
from osgeo import gdal
from qgis.core import *
from qgis.utils import *
from PyQt5.QtCore import QVariant

os.environ["QT_QPA_PLATFORM"] = "offscreen"
QgsApplication.setPrefixPath("/usr/lib/qgis", True) # Supply path to qgis install location
qgs = QgsApplication([], False) # Create a reference to the QgsApplication. Setting the second argument to False disables the GUI.
qgs.initQgis() # Load providers

sys.path.append('/usr/share/qgis/python/plugins')
# import processing *after* initializing the application
import processing
from processing.core.Processing import Processing
Processing.initialize()

from osgeo import gdal

PATH_TMP = '/home/nikola/4_north_america/data/tmp_files'

def print_progress(subject_function:callable, function_args:list, count:int, total_count:int, fun_times):
    start_time = time.time()
    result = subject_function(*function_args)
    fun_time = round(time.time() - start_time, 2)
    # results.append((lat, lon, city, country, result, fun_time))
    fun_times.append(fun_time)

    # pretty print out for task completeness
    if count % 10 == 0:
        mean_time = mean(fun_times)
        time_left = round(((total_count-count)*mean_time)/60,2)

    percentage = round(count / total_count * 100, 2)
    sys.stdout.write('\r')
    sys.stdout.write(f'task for {subject_function} completeness: {percentage:.2f}% time left: {time_left:.2f} mins')
    sys.stdout.flush()
    return result, fun_times


def clip_extent(lat, lon, dis):
    #Latitude:  1 deg = 110.54 km #Longitude: 1 deg = 111.320*cos(latitude) km
    #dis = 10 #u km, konstantan radijus okolo mjesta od interesa
    R = 6371 #Radius of the Earth in [km]
    brng = 1.57 #Bearing is 90 degrees converted to radians.
    #0.01oO = 1.11km -> * 10 -> 0.1oO = 11.1km -> *10 -> 1oO = 111km, prosjek radiusa max vjetra je 25-500 km-> 250 km = 2.25oO
    lat_delta = abs(dis/111) #in degrees oO

    lat1 = math.radians(lat) #Current lat point converted to radians
    lon1 = math.radians(lon) #Current long point converted to radians

    lat2 = math.asin(math.sin(lat1)*math.cos(dis/R) + math.cos(lat1)*math.sin(dis/R)*math.cos(brng)) #lat2 potrebna za racunanje lon2
    lon_W = lon1 - (math.atan2(math.sin(brng)*math.sin(dis/R)*math.cos(lat1),math.cos(dis/R)-math.sin(lat1)*math.sin(lat2))) # zapadna granicna longituda
    lon_E = lon1 + (math.atan2(math.sin(brng)*math.sin(dis/R)*math.cos(lat1),math.cos(dis/R)-math.sin(lat1)*math.sin(lat2))) # istocna granicna longituda
    lon_W = math.degrees(lon_W) #zapadna granicna longituda u stupnjevima oO
    lon_E = math.degrees(lon_E) #istocna granicna longituda u stupnjevima oO

    return lon_W,lat-lat_delta,lon_E,lat+lat_delta  #x_min,y_min,x_max,y_max


def get_pixel_value_from_geotiff(geotiff_path, lat, lon):
    try:
        with rasterio.open(geotiff_path) as src:
            x, y = src.index(lon, lat)
            
            # Read the pixel value at the specified coordinates
            value = src.read(1, window=((int(y), int(y) + 1), (int(x), int(x) + 1)))
            
            # Extract the actual pixel value (assuming it's a single band)
            pixel_value = value[0, 0]
            
            return pixel_value
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None


def get_altitude_from_dem_30m(lat, lon, path_tmp = PATH_TMP, infile = '', interpolate = True):

    """
    18.07.2022. Dohvacanje podatka o nadmorskoj visini iz TIF DEM30 m (cijela Europa)
    12.09.2022. Napravljena korekcija ukoliko je tocka izvan domene jer domena TIF-a ne odgovara domeni modela VRACA SE nan
    19.09.2022. Dodana opcija interpolacije da interpolira TIF ako je interpolacija na True (True je definirano unaprijed)
    """

    dem, err = subprocess.Popen(f'exec gdallocationinfo -valonly -xml -wgs84 {infile} {lon} {lat}', stdout=subprocess.PIPE, shell=True, universal_newlines=True).communicate()

    if dem == "</Report>":
        # izvan domene europe
        return np.nan
    else:

        dem = dem.splitlines()

        if interpolate == True:

            pixel = int(dem[0][dem[0].find('pixel="')+len('pixel="'):dem[0].rfind('" line')])
            line = int(dem[0][dem[0].find('line="')+len('line="'):dem[0].rfind('">')])
            # subset_size = 12
            subset_size = 2 # -> cini se da je svejedno jel 2 ili 12

            out, err = subprocess.Popen('gdal_translate -outsize {} {} -r bilinear -srcwin {} {} {} {} {} {}/subset_eudem_dem_3035_europe.tif'.format(subset_size*4, subset_size*4, pixel-1, line-1, subset_size, subset_size, infile, path_tmp), stdout=subprocess.PIPE, shell=True, universal_newlines=True).communicate()
            dem_interpolate, err = subprocess.Popen('gdallocationinfo -xml -wgs84 {}/subset_eudem_dem_3035_europe.tif {} {}'.format(path_tmp, lon, lat), stdout=subprocess.PIPE, shell=True, universal_newlines=True).communicate()
            dem_interpolate = dem_interpolate.splitlines()

            return round(float(dem_interpolate[2][dem_interpolate[2].find('<Value>')+len('<Value>'):dem_interpolate[2].rfind('</Value>')]), 2)

        else:
            # vracaj neinterpoliranu vrijednost
            return round(float(dem[2][dem[2].find('<Value>')+len('<Value>'):dem[2].rfind('</Value>')]), 2)

def get_tif_pixel_value(lat, lon, input_raster, in_domain_nodata_code=None, out_of_domain_code=None):
    '''
    23.08.2023. (Lana) Funkcija koja prema lat, lon i kodova za in domain no data i out of domain
                        iz TIF rastera vadi vrijednost piksela i vraca njegovu vrijednost. 
                        Ako je piksel in domain, no data -> vraca None
                        Ako je piksel out of domain -> vraca np.nan
    '''
    
    args = ['gdallocationinfo', '-valonly', '-wgs84', '-b', '1', input_raster, str(lon), str(lat)]
    result = subprocess.run(args, capture_output=True, text=True).stdout
    result_clean = result.rstrip("\n")
    # result1 = result1.rstrip("\n255")
    if result_clean: # ako string nije prazan
        pixel_value=float(result_clean)
        
    else: # ako je (None je), izvan domene je -> vrati np.nan
        return np.nan

    # print('get_tif_pixel_value=', pixel_value, infile_tif)
    if pixel_value == in_domain_nodata_code:
        return None
    
    elif pixel_value == out_of_domain_code: # ako je izvan domene, nekad vrati sifru (npr -23767)
        return np.nan
    
    else:
        return np.round(pixel_value, 2)

# TODO: testing the landuse percentage??

def land_cover_percentage(lat, lon, input_path, dis=0.5):
    """
    30.8.2023 u funkciji se ne poziva više QGIS, već se većina temelji na gdalu. 
    """

    names_list = ['grassland','leaf_type','SWF'] 
    x_min, y_min, x_max, y_max = clip_extent(lat,lon,dis)
    count = 0 
    dictionary = {}

    xml_string, err = subprocess.Popen([f'gdallocationinfo -xml -wgs84 {input_path} {lon} {lat}'], 
                                  stdout=subprocess.PIPE, shell=True).communicate()

    root = ET.fromstring(xml_string)
    value_element = root.find('.//Value')
    try:
        result = int(value_element.text)
    except AttributeError as e:
        print(f"Error occurred: {e} therefore, results for given loc are {value_element}")
        result = None

    if result == None:
        dictionary = {'grassland': np.nan, 'leaf_type': np.nan, 'SWF': np.nan}
    
    else:
        if result == 255: 
            dictionary = {'grassland': None, 'leaf_type': None, 'SWF': None}

        else: 
            input_file = os.path.basename(input_path)
            output_path = f'{PATH_TMP}/subset_{input_file}.tif'
            out, err = subprocess.Popen(f'gdal_translate -projwin {x_min} {y_max} {x_max} {y_min} -projwin_srs EPSG:4326 {input_path} {output_path}', 
                                        stdout=subprocess.PIPE, shell=True, universal_newlines=True).communicate() 
            
            raster = gdal.Open(output_path)
            values, counts  = np.unique(raster.GetRasterBand(1).ReadAsArray(), return_counts=True)

            total_pixels = np.sum(counts)
            percentages = (counts / total_pixels) * 100

            pixel_percentage_dict = dict(zip(values, percentages))

            # check non_zero_sum values
            non_zero_sum = round(sum(percentages[1:]),2)

            # esa world cover classes: 
            # "Tree cover", "Shrubland", "Grassland", "Cropland", 
            # "Built-up", "Bare / sparse vegetation”, “Snow and Ice”, 
            # “Permanent water bodies”, “Herbaceous Wetland”, “Mangrove” and “Moss and lichen"

            key_mapping = {
                30: "grassland",
                10: "leaf_type",
                20: "SWF"
            }
            for old_key, new_key in key_mapping.items():
                dictionary[new_key] = round(pixel_percentage_dict.get(old_key, 0),2)   

    return dictionary

def clip_by_extent_gdal(lat, lon, output_raster, input_raster, dis = 0.5):
    
    # options = gdal.TranslateOptions(xRes=1, yRes=1, projWin = [lon+dis, lat-dis, lon-dis, lat+dis], projWinSRS = 4326)
    # options = (f"gdal_translate -projwin {lon+dis} {lat-dis} {lon-dis} {lat+dis} -projwin_srs EPSG:4326 {input_raster} {output_raster}")
    # options = dict(format="GTiff", xRes=1, yRes=1, projWin = [lon-dis, lat+dis, lon+dis, lat-dis], projWinSRS = "EPSG:4326")
    # output = gdal.Translate(output_raster, input_raster, **options)

    out, err = subprocess.Popen(f'gdal_translate -projwin {lon-dis} {lat+dis} {lon+dis} {lat-dis} -projwin_srs EPSG:4326 {input_raster} {output_raster}', 
                            stdout=subprocess.PIPE, shell=True, universal_newlines=True).communicate() 
    
    print(f"Gdal translate for {output_raster} done! coordinates are {lat}, {lon}")
    return out


def query_pixel_value_rasterio(geotiff_path, lat, lon):
    """
    Query pixel value for given latitude and longitude coordinates from a GeoTIFF file using rasterio.

    Parameters:
        geotiff_path (str): Path to the GeoTIFF file.
        lat (float): Latitude coordinate.
        lon (float): Longitude coordinate.

    Returns:
        pixel_value: Pixel value at the specified coordinates.
    """
    # Open the GeoTIFF file
    with rasterio.open(geotiff_path) as src:
        # Transform latitude and longitude to pixel coordinates
        col, row = transform(src.crs, {'init': 'EPSG:4326'}, [lon], [lat])
        print(col, row)
        row, col = int(row[0]), int(col[0])
        print(col, row)

        # Read the pixel value at the specified coordinates
        # pixel_value = src.read(1,window=(row, col))
        # print(pixel_value)
        # window = Window.from_slices((row, row+1), (col, col+1))
        pixel_value = src.read(1, window=((row, row + 0.5), (col, col + 0.5)))
        print(pixel_value)

    return pixel_value[0][0]  # Extract the pixel value from the array