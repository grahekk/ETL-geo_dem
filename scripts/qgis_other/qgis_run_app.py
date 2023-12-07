import sys
import os
import time
import yaml
import psycopg2

sys.path.append("/home/nikola/4_north_america/scripts/")
import settings

# from .qgis_create_graphics import create_and_export_map
# from .qgis_create_grid import create_polygon_grid, clip_by_layer

config = settings.get_config()
conn_parameters = settings.get_conn_parameters()
schema = settings.get_schema()

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


# Use the app:
def run_qgis_app(fun:callable, args:list):
    """
        Run a specified function within a QGIS application instance.

        This function initializes a QGIS application without the GUI, calls the specified
        function with the provided arguments, and exits the QGIS application when the
        function has completed its tasks.

        Parameters:
            fun (callable): The function to be executed within the QGIS application instance.
            args (list): A list of arguments to be unpacked and passed to the function.
    """
    # initialize the app
    # qgs = QgsApplication([], False) # Create a reference to the QgsApplication. Setting the second argument to False disables the GUI.
    # qgs.initQgis() # Load providers
    
    # call the function and unpack the arguments provided in the list
    fun(*args)
    # Exit the QGIS application when done
    qgs.exitQgis()
    print(f"Task '{fun.__name__}' done and qgis app closed!")


def create_and_export_map(conn_parameters, schema:str, layer_path:str, output_image_path:str, layout_path:str, image_size_mm=(300,300), dpi=300, original_color=True, zoom = "query"):
    """
    Load a raster layer from a VRT (Virtual Raster) or other supported format,
    create a map layout, and export it as an image.

    Args:
        layer_path (str): The file path to the raster layer or VRT file to be loaded.
        output_image_path (str): The file path where the exported image will be saved.
        layout_path (str): The file path where the layout template will be saved.
        image_size_mm (tuple): A tuple (width_mm, height_mm) specifying the size of the layout in millimeters.
        dpi (int): The resolution in dots per inch (DPI) for the exported image.
        original_color (bool): to keep or to try to repaint the layer in new color scale
        zoom (str): choose type of zoom, "query" zoom takes extent from layer in local pgdb

    Returns:
        None: This function does not return a value but exports the map layout and layout as an image.
    """
    project = QgsProject.instance()

    # load the file as layer into project
    fileInfo = QFileInfo(layer_path)
    path = fileInfo.filePath()

    # Check if the file is a raster or a vector
    if path.lower().endswith(('.tif', '.tiff', '.asc', '.vrt')):
        # It's a raster
        layer = QgsRasterLayer(path, "Layer")
    else:
        # It's a vector
        layer = QgsVectorLayer(path, "Layer", "ogr")

    if not layer.isValid():
        print("Failed to load the layer!")
        # qgs.exitQgis()
        return
    else:
        project.addMapLayer(layer)

    if original_color == False:
        # change_layer_color(layer) # if raster layer
        # layer = change_vector_layer_color(layer, "values")

        style_file_path = "/home/nikola/4_north_america/scripts/qgis_other/values_test_dem_colored.qml"

        with open(style_file_path) as style_file:
            style_content = style_file.read()

        document = QDomDocument()
        document.setContent(style_content)
        # layer.loadNamedStyle(style_file)
        # layer.readXml(document, QgsReadWriteContext())
        # layer.readStyle(document, context = QgsReadWriteContext())
        # layer.readSld(document)
        layer.readSymbology(document, errorMessage = "error me", context=QgsReadWriteContext())
        layer.triggerRepaint()
        layer.reload()

    # sizes setting
    image_width_mm=image_size_mm[0]
    image_height_mm=image_size_mm[1]
    map_size = QgsLayoutSize(image_width_mm, image_height_mm)

    # Create a layout
    layout = QgsLayout(project)
    layout.initializeDefaults()

    # Set the layout size
    layout.pageCollection().pages()[0].setPageSize(map_size) #layout size in milimeters

    # create map layout item
    map_frame = QgsLayoutItemMap(layout)
    map_frame.attemptResize(map_size)
    # add layer to the map and then fix the zoom
    map_frame.setLayers([layer])

    if zoom == "query":
        with psycopg2.connect(**conn_parameters) as conn:
            with conn.cursor() as cursor:
                cursor.execute("SELECT ST_asText(ST_Envelope(geometry)) FROM osm.world_continents_boundaries as cna WHERE cna.continent = 'North America'")
                wkt_string = cursor.fetchone()[0]
        
        # zooming_layer = QgsVectorLayer(f"?query=SELECT ST_GeomFromText('{wkt_string}')", "extent", "virtual")
            zoom = QgsRectangle.fromWkt(wkt_string)

    elif zoom == "layer":
        zooming_layer = layer
        zoom = zooming_layer.extent()
    
    map_frame.zoomToExtent(zoom)
    map_frame.create(layout)
    # add map to the print layout and set size
    layout.addLayoutItem(map_frame)
    map_frame.setFixedSize(map_size)

    # Set the settings to export the layout as an image
    exporter = QgsLayoutExporter(layout)
    context = QgsReadWriteContext()
    image_settings = QgsLayoutExporter.ImageExportSettings()
    image_settings.dpi = dpi
    image_settings.imageDpi = dpi
    image_settings.imageQuality = 100

    # exporting the products
    # layout.saveAsTemplate(layout_path, context)
    exporter.exportToImage(output_image_path, image_settings)
    return print(f"Exporting done for image: {output_image_path} and print layout: {layout_path}!")


# flows here
def grapical_map_flow(file_path:str, zoom = "query"):
    """
    Generate a graphical map using QGIS with specified parameters.

    This function generates a graphical map by utilizing QGIS. It selects a file
    from a list of available files and retrieves the file path from a configuration.
    The map is created based on a predefined layout, and the resulting image is
    saved to a specified location.
    """
    files = ["usa_dem_vrt", "cec_lc_north_america_tif", "canada_dem_vrt", "canada_WMS_dem_vrt"]
    # choose a file from the list of files, and the rest is done by pulling path from config
    # file = files[0]
    # file_path = config[file]

    # or reference it by hand
    file = os.path.basename(file_path)

    output_image_path = f"./{file}_graphics.jpeg"
    layout_path = "./layout.qpt"
    image_size_mm=(300,300)
    dpi=300 
    original_color=True
    args_graphics = [conn_parameters, schema, file_path, output_image_path, layout_path, image_size_mm, dpi, original_color, zoom]
    return run_qgis_app(create_and_export_map, args_graphics)


def grid_creation_flow(output:str, boundary:str, extension = "shp"):
    """
    Create a polygon grid within the specified extent and save it to a file.

    This function generates a polygon grid within the specified extent and saves
    it to a file with the specified format (default is GeoPackage, "gpkg"). The grid
    is created using the provided extent from a configuration file.

    Parameters:
    - output (str): The base name for the output file, excluding the extension.
    - extent (str): The key for the extent in the configuration file.
    - extension (str, optional): The file format extension (e.g., "gpkg", "shp").

    Returns:
    - str: The filename of the saved grid layer.
    """
    grid_path = f"{output}.{extension}"
    layer_extent = config[boundary]
    grid = create_polygon_grid(grid_path, layer_extent)
    clipped_grid = clip_by_layer(grid, layer_extent, f"{output}.{extension}")
    print(f"grid layer clipped and saved into {clipped_grid} !")
    return f"{clipped_grid}"


if __name__=="__main__":
    # grid_file = grid_creation_flow("esa_usa_grid", "NA_aspect_categorized")
    # grid_file = grid_creation_flow("continents_grid_clipped", "world_continents_boundaries_file")
    # testing_file = "/home/nikola/4_north_america/scripts/testing/test_results.shp"
    # graphical_file = "/mnt/volume-nbg1-1/shared/nikola/ESA_global_dem_90m_vrt/ESA_global_dem_90m.vrt_0.vrt"
    # graphical_file = "/home/nikola/continents_grid_clipped_enrich.shp"
    graphical_file = config["NA_slope"]
    graphical_file = "/mnt/volume-nbg1-1/shared/nikola/ESA_global_dem_90m_vrt/ESA_usa_dem_90m.vrt"
    grapical_map_flow(graphical_file)

