import sys
import os
import time
import yaml
import datetime
import subprocess
import re
# from .pipeline_transform_vrt_gdal import transform_raster
from statistics import mean
# from . import pipeline_transform_vrt_gdal

from .pipeline_transform_vrt_gdal import make_file_path
import settings
config = settings.get_config()

sys.path.append('/usr/lib/python3/dist-packages')
from qgis.core import * 
from qgis.gui import *
from qgis.PyQt.QtGui import QImage, QPainter, QColor
from qgis.PyQt.QtCore import QSize, QFileInfo
from qgis.utils import *
from PyQt5.QtCore import QVariant


# Initialize the QGIS application
os.environ["QT_QPA_PLATFORM"] = "offscreen"
QgsApplication.setPrefixPath("/usr/lib/qgis", True) # Supply path to qgis install location

# initialize the app
qgs = QgsApplication([], False) # Create a reference to the QgsApplication. Setting the second argument to False disables the GUI.
qgs.initQgis() # Load providers

sys.path.append('/usr/share/qgis/python/plugins')
# import processing *after* initializing the application
import processing
from processing.core.Processing import Processing
from qgis.analysis import QgsNativeAlgorithms
Processing.initialize()


def transform_geomorphon_qgis(dem_file, output_file, search = 90, skip = 3, flat = 5, dist = 6, excplicit_save = False):
    """
    Function uses GRASS GIS geomorphon algorithm for transforming raster tiles or virtual raster.
    Output is raster file that categorises relief into 8 geomorphology types.
    """
    project = QgsProject.instance()
    layer = QgsRasterLayer(dem_file, "Layer")

    if not layer.isValid():
        print("Failed to load the layer!")
        # qgs.exitQgis()
        return
    else:
        project.addMapLayer(layer)
        geomorphon_layer = processing.run("grass7:r.geomorphon", {'elevation':layer,
                                            'search':search,
                                            'skip':skip,
                                            'flat':flat,
                                            'dist':dist,
                                            'forms':output_file,
                                            '-m':False,
                                            '-e':False,
                                            'GRASS_REGION_PARAMETER':None,
                                            'GRASS_REGION_CELLSIZE_PARAMETER':0,
                                            'GRASS_RASTER_FORMAT_OPT':'',
                                            'GRASS_RASTER_FORMAT_META':''})
        
        if excplicit_save == True:
            # add layer for faster processing
            geomorphon_layer = QgsRasterLayer(geomorphon_layer['forms'], "geomorphon")
            project.addMapLayer(geomorphon_layer)
            # create exporter and export raster file using pipe
            file_writer = QgsRasterFileWriter(output_file)
            pipe = QgsRasterPipe()
            provider = geomorphon_layer.dataProvider()

            if not pipe.set(provider.clone()):
                print("Cannot set pipe provider")

            file_writer.writeRaster(
                pipe,
                provider.xSize(),
                provider.ySize(),
                provider.extent(),
                provider.crs()
                )
            
    return
    

def raster_rescale_percent_to_degree(input_raster, output_raster):
    project = QgsProject.instance()
    layer = QgsRasterLayer(input_raster, "Layer_to_be_rescaled")

    if not layer.isValid():
        print("Failed to load the layer!")
        # qgs.exitQgis()
        return
    else:
        project.addMapLayer(layer)
        start_time = time.time()
        rescaled_layer = processing.run("native:rescaleraster", {'INPUT':input_raster,
                                                'BAND':1,
                                                'MINIMUM':0,
                                                'MAXIMUM':90,
                                                'NODATA':None,
                                                'OUTPUT':config["NA_slope_tmp"]})
        rescaled_layer = QgsRasterLayer(rescaled_layer['OUTPUT'], "rescaled_layer")
        project.addMapLayer(rescaled_layer)
        # create exporter and export raster file using pipe
        file_writer = QgsRasterFileWriter(output_raster)
        pipe = QgsRasterPipe()
        provider = rescaled_layer.dataProvider()

        if not pipe.set(provider.clone()):
            print("Cannot set pipe provider")

        file_writer.writeRaster(
            pipe,
            provider.xSize(),
            provider.ySize(),
            provider.extent(),
            provider.crs()
            )

        end_time = time.time()
        total_time = round((end_time-start_time)/60, 2)
        print(f"Layer processing to geomorphon done! time (dd:hh:mm.ss): {datetime.timedelta(seconds=total_time)} - raw time:({total_time})")
        print(f"results saved to {output_raster}")
        qgs.exitQgis()


def transform_slope_flow():
    """    
    Compute slope using the gdaldem command from shell. 
    Rescaling algorithm is run later to adjust values from percentages into degrees.
    """
    start_time = time.time()
    dem_file = config["esa_na_dem_90_vrt"]
    slope_file = config["NA_slope"]
    rescaled_slope_file = config["NA_slope_rescaled"]

    # transform_raster(dem_file, slope_file, transform_type="slope")
    raster_rescale_percent_to_degree(slope_file, rescaled_slope_file)

    end_time = time.time()
    total_time = round((end_time-start_time)/60, 2)

    print(f"Slope computing and rescaling done in time(dd:hh:mm.ss): {datetime.timedelta(seconds=total_time)} - raw time:({total_time})")

def geomorphon_chunky(central_file, neighbour_vrt, output_dir, search = 90, skip = 3, flat = 5, dist = 6):
    """
    Function calls geomorphon algorithm in qgis (from grass gis).
    This is variation of algorithm in qgis - input layers are neighbors DEMs also.
    
    Parameters:
        central_file(str): path to DEM tile for which geomorphon is calculated
        neighbour_vrt(str): path to created vrt consisting of `central_file` and it's 8 neighbouring tiles
    
    Standard parameters are:
        Search radius = 90
        Skip radius = 3
        flatness threshold = 5
        flatness distance = 6
    """
    project = QgsProject.instance()
    layer = QgsRasterLayer(neighbour_vrt, "Neighbour vrt layer")
    central_layer = QgsRasterLayer(central_file, "Tile dem layer")

    # get pixels sizes of raster so that later extent can be redefined
    pixel_size_x = central_layer.rasterUnitsPerPixelX()
    pixel_size_y = central_layer.rasterUnitsPerPixelY()

    if not layer.isValid():
        print(f"Failed to load the layer!")
        return
    else:
        project.addMapLayer(layer)
        algorithm_params = {'elevation':layer,
                            'search':search,
                            'skip':skip,
                            'flat':flat,
                            'dist':dist,
                            'forms':output_dir,
                            '-m':False,
                            '-e':False,
                            'GRASS_REGION_PARAMETER':None,
                            'GRASS_REGION_CELLSIZE_PARAMETER':0,
                            'GRASS_RASTER_FORMAT_OPT':'',
                            'GRASS_RASTER_FORMAT_META':''}

        # Get the extent of the input raster
        chunk_extent = central_layer.extent()

        # Set the chunk output path
        chunk_output_path = make_file_path(central_file, output_dir, sufix="_geomorphon.tif")
        algorithm_params['forms'] = chunk_output_path

        # Set the extent for the r.geomorphon command
        algorithm_params['GRASS_REGION_PARAMETER'] = f"{chunk_extent.xMinimum()-pixel_size_x*100},{chunk_extent.xMaximum()+pixel_size_x*100}," \
                                                        f"{chunk_extent.yMinimum()-pixel_size_y*100},{chunk_extent.yMaximum()+pixel_size_y*100}"

        # Run the r.geomorphon algorithm for the current dem chunk/tile
        context = QgsProcessingContext() 
        feedback = QgsProcessingFeedback()
        geomorphon_result = processing.run('grass7:r.geomorphon', algorithm_params, context=context, feedback=feedback)
        return chunk_output_path


if __name__ == "__main__":
    qgs.exitQgis()
