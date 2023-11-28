import sys
import os
import time
import datetime
import yaml

import settings
config = settings.get_config()

sys.path.append('/usr/lib/python3/dist-packages')
from qgis.core import * 
from qgis.gui import *
from qgis.PyQt.QtGui import QImage, QPainter, QColor
from qgis.PyQt.QtCore import QSize, QFileInfo
from qgis.utils import *
from PyQt5.QtCore import QVariant


# # Initialize the QGIS application
# os.environ["QT_QPA_PLATFORM"] = "offscreen"
# QgsApplication.setPrefixPath("/usr/lib/qgis", True) # Supply path to qgis install location

# # initialize the app
# qgs = QgsApplication([], False) # Create a reference to the QgsApplication. Setting the second argument to False disables the GUI.
# qgs.initQgis() # Load providers

# sys.path.append('/usr/share/qgis/python/plugins')
# # import processing *after* initializing the application
# import processing
# from processing.core.Processing import Processing
# Processing.initialize()

from osgeo import gdal 

from .pipeline_transform_vrt_gdal import absolute_file_paths

import settings
config = settings.get_config()


def gdal_build_vrt(raster_files_dir:list, output_vrt:str, resample_method = "average"):
    """
    Create a Virtual Raster (VRT) dataset from a list of raster files.

    Args:
        raster_files (list): The path to locate the input raster files
        output_vrt (str)
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

    input_files = absolute_file_paths(raster_files_dir, ".tif") # returns generator object
    files_list = sorted(list(input_files))

    start_time = time.time()
    vrt_options = gdal.BuildVRTOptions(resampleAlg=resample_method, addAlpha=True, hideNodata=False)
    my_vrt = gdal.BuildVRT(output_vrt, files_list, options=vrt_options)
    my_vrt = None
    fun_time = round(time.time()-start_time, 2)
    return print(f"Gdal building vrt done for file: {output_vrt} in time: {datetime.timedelta(seconds=fun_time)}")