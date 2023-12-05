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

sys.path.append("/home/nikola/4_north_america/GeoDataPump/scripts/")
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
        
        input(f"please press any key to continue, time is {round((time.time()-start_time)/60, 2)}")

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

    # dem_file = "/mnt/volume-nbg1-1/shared/nikola/ESA_global_dem_90m/Copernicus_DSM_30_N45_00_E015_00_DEM.tif"
    # slope_file = "/mnt/volume-nbg1-1/shared/nikola/NA_slope_N45_E015.tif"

    # transform_raster(dem_file, slope_file, transform_type="slope")
    raster_rescale_percent_to_degree(slope_file, rescaled_slope_file)

    end_time = time.time()
    total_time = round((end_time-start_time)/60, 2)

    print(f"Slope computing and rescaling done in time(dd:hh:mm.ss): {datetime.timedelta(seconds=total_time)} - raw time:({total_time})")


# Define a custom processing algorithm to run r.geomorphon
class RGeomorphonAlgorithm(QgsProcessingAlgorithm):
    INPUT_RASTER = 'INPUT_RASTER'
    OUTPUT_RASTER = 'OUTPUT_RASTER'

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return RGeomorphonAlgorithm()

    def name(self):
        return 'r_geomorphon_algorithm'

    def displayName(self):
        return self.tr('r.geomorphon Algorithm')

    def group(self):
        return self.tr('Example Group')

    def groupId(self):
        return 'example_group'

    def shortHelpString(self):
        return self.tr('Example algorithm to run r.geomorphon on a raster layer.')

    def initAlgorithm(self, config=None):
        # Input raster layer
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.INPUT_RASTER,
                self.tr('Input Raster Layer')
            )
        )

        # Output raster layer
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.OUTPUT_RASTER,
                self.tr('Output Raster Layer')
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        # Get input and output raster layers
        input_raster = self.parameterAsRasterLayer(parameters, self.INPUT_RASTER, context)
        output_raster = self.parameterAsRasterLayer(parameters, self.OUTPUT_RASTER, context)


        algorithm_params = {'elevation':input_raster.source(),
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
                                            'GRASS_RASTER_FORMAT_META':''}


       # Get the extent of the input raster
        extent = input_raster.extent()

        # Specify the chunk size (adjust as needed)
        chunk_size = QSize(1201, 1201)

        # Loop through the extent in chunks
        for x in range(int(extent.xMinimum()), int(extent.xMaximum()), chunk_size.width()):
            for y in range(int(extent.yMinimum()), int(extent.yMaximum()), chunk_size.height()):
                # Create a new extent for the chunk
                chunk_extent = QgsRectangle(x, y, x + chunk_size.width(), y + chunk_size.height())

                # Set the chunk output path
                chunk_output_path = os.path.join(algorithm_params['output'], f'chunk_{x}_{y}.tif')

                # Set the extent for the r.geomorphon command
                algorithm_params['GRASS_REGION_PARAMETER'] = f"{chunk_extent.xMinimum()},{chunk_extent.xMaximum()}," \
                                                             f"{chunk_extent.yMinimum()},{chunk_extent.yMaximum()}"

                # Run the r.geomorphon algorithm for the current chunk
                processing.run('grass7:r.geomorphon', algorithm_params, context=context, feedback=feedback)

                # Load the processed chunk as a QgsRasterLayer
                chunk_raster = QgsRasterLayer(chunk_output_path, 'chunk_raster')

                # Save the processed chunk to the output file
                writer = QgsRasterFileWriter(chunk_output_path)
                writer.writeRaster(chunk_raster.renderer(), chunk_extent, chunk_size)

        return {self.OUTPUT_RASTER: output_raster}

def create_virtual_raster(self, extent, output_folder):
        # Get the bounds of the current extent
        xmin, ymin, xmax, ymax = extent.xMinimum(), extent.yMinimum(), extent.xMaximum(), extent.yMaximum()

        # Get the tile coordinates for the current extent
        tile_x = int(xmin)
        tile_y = int(ymin)

        # Create a list to store the paths of the input tiles
        input_tiles = []

        # Iterate over the surrounding 8 tiles
        for i in range(tile_x - 1, tile_x + 1):
            for j in range(tile_y - 1, tile_y + 1):
                # Construct the file name for each surrounding tile
                tile_name = f'Copernicus_DSM_30_S{j:02d}_W{i:03d}_DEM.tif'
                tile_path = os.path.join(output_folder, tile_name)

                # Add the tile path to the input_tiles list
                input_tiles.append(tile_path)
        
        # Create a virtual raster for the current tile and its surrounding 8 tiles
        vrt_path = os.path.join(output_folder, f'virtual_raster_{tile_x}_{tile_y}.vrt')
        QgsProcessingUtils.generateVRT(vrt_path, input_tiles, False)


def geomorphon_chunky(central_file, neighbour_vrt, context = QgsProcessingContext(), feedback = QgsProcessingFeedback()):
    """
    Variation of using geomorphon algorithm in qgis.
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
                            'GRASS_RASTER_FORMAT_META':''}

        # Get the extent of the input raster
        chunk_extent = central_layer.extent()

        # Set the chunk output path
        output_file = f'{os.path.splitext(os.path.basename(central_file))[0]}_geomorphon.tif'
        chunk_output_path = os.path.join(config['tmp_dir'], output_file)
        algorithm_params['forms'] = chunk_output_path

        # Set the extent for the r.geomorphon command
        algorithm_params['GRASS_REGION_PARAMETER'] = f"{chunk_extent.xMinimum()-pixel_size_x*4},{chunk_extent.xMaximum()+pixel_size_x*4}," \
                                                        f"{chunk_extent.yMinimum()-pixel_size_y*4},{chunk_extent.yMaximum()+pixel_size_y*4}"

        # Run the r.geomorphon algorithm for the current chunk
        geomorphon_result = processing.run('grass7:r.geomorphon', algorithm_params, context=context, feedback=feedback)
        # return print(f"processing done for {output_file}")


if __name__ == "__main__":
    dem_file = config["esa_na_dem_90_vrt"]
    output_file = config["NA_geomorphon"]
    # dem_file = "/mnt/volume-nbg1-1/shared/nikola/ESA_global_dem_90m/Copernicus_DSM_30_N45_00_E015_00_DEM.tif"
    # output_file = "/mnt/volume-nbg1-1/shared/nikola/NA_geomorphon_N45_E015.tif"
    # transform_geomorphon_qgis(dem_file, output_file)
    # transform_slope_flow()

    file = "Copernicus_DSM_30_S20_00_W060_00_DEM.tif"

    qgs.exitQgis()
