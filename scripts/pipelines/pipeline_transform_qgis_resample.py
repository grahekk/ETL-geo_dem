import sys
import os
import time
import yaml
import datetime
import subprocess
# from .pipeline_transform_vrt_gdal import transform_raster

sys.path.append("/home/nikola/4_north_america/scripts/")
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
Processing.initialize()


def transform_geomorphon_qgis(dem_file, output_file):
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
        start_time = time.time()
        geomorphon_layer = processing.run("grass7:r.geomorphon", {'elevation':layer,
                                            'search':90,
                                            'skip':3,
                                            'flat':5,
                                            'dist':6,
                                            'forms':config["NA_geomorphon_tmp"],
                                            '-m':False,
                                            '-e':False,
                                            'GRASS_REGION_PARAMETER':None,
                                            'GRASS_REGION_CELLSIZE_PARAMETER':0,
                                            'GRASS_RASTER_FORMAT_OPT':'',
                                            'GRASS_RASTER_FORMAT_META':''})
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

        end_time = time.time()
        total_time = round((end_time-start_time)/60, 2)
        print(f"Layer processing to geomorphon done! time (dd:hh:mm.ss): {datetime.timedelta(seconds=total_time)} - raw time:({total_time})")
        print(f"results saved to {output_file}")
    

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


def transform_slope_flow():
    """    
    Compute slope using the gdaldem command from shell. 
    Rescaling algorithm is run later to adjust values from percentages into degrees.
    """
    start_time = time.time()
    dem_file = config["esa_na_dem_90_vrt"]
    slope_file = config["NA_slope"]
    rescaled_slope_file = config["NA_slope_rescaled"]

    # dem_file = "/mnt/volume-nbg1-1/shared/nikola/ESA_global_dem_90m/Copernicus_DSM_30_N45_00_E015_00_DEM.tif"
    # slope_file = "/mnt/volume-nbg1-1/shared/nikola/NA_slope_N45_E015.tif"

    # transform_raster(dem_file, slope_file, transform_type="slope")
    raster_rescale_percent_to_degree(slope_file, rescaled_slope_file)

    end_time = time.time()
    total_time = round((end_time-start_time)/60, 2)

    print(f"Slope computing and rescaling done in time(dd:hh:mm.ss): {datetime.timedelta(seconds=total_time)} - raw time:({total_time})")

if __name__ == "__main__":
    dem_file = config["esa_na_dem_90_vrt"]
    output_file = config["NA_geomorphon"]
    # dem_file = "/mnt/volume-nbg1-1/shared/nikola/ESA_global_dem_90m/Copernicus_DSM_30_N45_00_E015_00_DEM.tif"
    # output_file = "/mnt/volume-nbg1-1/shared/nikola/NA_geomorphon_N45_E015.tif"
    transform_geomorphon_qgis(dem_file, output_file)
    # transform_slope_flow()

    qgs.exitQgis()
