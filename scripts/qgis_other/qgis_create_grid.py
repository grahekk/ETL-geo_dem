import sys
import os

sys.path.append('/usr/lib/python3/dist-packages')
from qgis.core import (QgsApplication, 
QgsVectorLayer, 
QgsVectorFileWriter, 
QgsCoordinateReferenceSystem, 
QgsProject, 
QgsCoordinateTransform, 
QgsVectorLayerExporter, 
QgsGeometry, 
QgsProcessingFeatureSourceDefinition, 
QgsProcessingContext, 
QgsProcessingFeedback,
QgsExpression, 
QgsFeatureRequest,
QgsExpressionContext,
QgsExpressionContextUtils,
edit,
QgsField)

from qgis.utils import *
from PyQt5.QtCore import QVariant

# Initialize the QGIS application
os.environ["QT_QPA_PLATFORM"] = "offscreen"
QgsApplication.setPrefixPath("/usr/lib/qgis", True) # Supply path to qgis install location
qgs = QgsApplication([], False) # Create a reference to the QgsApplication. Setting the second argument to False disables the GUI.
qgs.initQgis() # Load providers

sys.path.append('/usr/share/qgis/python/plugins')
# import processing *after* initializing the application
import processing
from processing.core.Processing import Processing
Processing.initialize()


def create_polygon_grid(output:str, bounding_box_layer:str, epsg='4326', grid_type = 0, spacing = 5):
    """
    Create a polygon grid layer within a specified bounding box and save it to file.

    This function generates a polygon grid layer based on the provided bounding box and
    specified parameters. The grid can be of different types (point, line, or polygon) and
    is saved to a file with the specified filename.

    Parameters:
    - output (str): The filename of the vector file where the grid will be saved.
    - bounding_box_layer (str): The path to the bounding box layer used to define the extent.
    - epsg (str, optional): The EPSG code of the coordinate reference system (CRS). Default is '4269'.
    - grid_type (int, optional): The type of grid to create (0 for points, 1 for lines, 2 for polygons). Default is 0.

    Returns:
    - str: The filename of the saved grid vector layer.
    """

    bounding_box_layer = QgsVectorLayer(bounding_box_layer, "bounding_box_layer", "ogr")
    QgsProject.instance().addMapLayer(bounding_box_layer)

    target_crs = QgsCoordinateReferenceSystem(f'EPSG:{epsg}')
    
    # Get the extent and CRS of the bounding box layer
    bounding_box = bounding_box_layer.extent()

    # Create the grid using processing algorithms, also in EPSG 102008
    grid = processing.run("native:creategrid", {
        'TYPE': grid_type,
        'EXTENT': f'{bounding_box.xMinimum()},{bounding_box.xMaximum()},{bounding_box.yMinimum()},{bounding_box.yMaximum()} [EPSG:{epsg}]',
        'HSPACING': spacing, 
        'VSPACING': spacing,  
        'HOVERLAY': 0,
        'VOVERLAY': 0,
        'CRS': target_crs,
        'OUTPUT': output
    })['OUTPUT']

    print(f"Grid completed. Saved to the file {output}")
    return output


def clip_by_layer(layer_to_be_clipped, overlay_layer, output:str):
    """
    Clip a vector layer using another vector layer and save the result to a file.

    This function clips a vector layer (`layer_to_be_clipped`) using another vector layer
    (`overlay_layer`) as the clipping boundary. The clipped result is saved to a file
    with the specified filename.

    Parameters:
    - layer_to_be_clipped (QgsVectorLayer): The vector layer to be clipped.
    - overlay_layer (str): The path to the overlay vector layer used as the clipping boundary.
    - output (str): The filename for the clipped vector layer output.

    Returns:
    - str: The filename of the saved clipped vector layer.
    """
    layer_clip = processing.run("native:clip", {'INPUT':layer_to_be_clipped,
                                   'OVERLAY': overlay_layer,
                                   'OUTPUT':'TEMPORARY_OUTPUT'})['OUTPUT']
    writer = QgsVectorFileWriter.writeAsVectorFormat(layer_clip, output, fileEncoding = "utf-8", driverName = "ESRI shapefile")
    return output


if __name__ == "__main__":
    create_polygon_grid("output.gpkg", "4_north_america/data/USA_dataset_extent_OT.012021.4269.1.kml")

    # qgis.exit()