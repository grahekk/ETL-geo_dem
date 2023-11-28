import sys
import os
import time
import yaml
import psycopg2
import settings

# load local config settings
config = settings.get_config()
conn_parameters = settings.get_conn_parameters()

sys.path.append('/usr/lib/python3/dist-packages')
from qgis.core import * 
from qgis.gui import *
from qgis.PyQt.QtGui import QImage, QPainter, QColor
from qgis.PyQt.QtCore import QSize, QFileInfo
from qgis.utils import *
from PyQt5.QtCore import QVariant
from PyQt5.QtXml import QDomDocument, QDomElement


# Initialize the QGIS application
os.environ["QT_QPA_PLATFORM"] = "offscreen"
QgsApplication.setPrefixPath("/usr/lib/qgis", True) # Supply path to qgis install location
# initialize the app
qgs = QgsApplication([], False) # Create a reference to the QgsApplication. Setting the second argument to False disables the GUI.
qgs.initQgis() # Load providers

sys.path.append('/usr/share/qgis/python/plugins')
# # import processing *after* initializing the application
# import processing
# from processing.core.Processing import Processing
# Processing.initialize()

#TODO: add zooming feature to lat lot coords for specified zoom
#TODO: fix/finish feature for adding colors (graduated colors) to raster and vector

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

    if zoom == "table": # query the layer to zoom to
        tablename = "world_continents_boundaries"
        continent = 'North America'
        uri = QgsDataSourceUri()
        uri.setConnection(*conn_parameters.values())
        uri.setDataSource (schema, tablename, "geometry", "continent = 'North America", '')
        
        zooming_layer=QgsVectorLayer(uri.uri(False), tablename, "postgres")
        if not zooming_layer.isValid():
            print(f"Failed to load the layer {zooming_layer}!")
            return
        QgsProject.instance().addMapLayer(zooming_layer)

    elif zoom == "layer":
        zooming_layer = config[""]
        QgsProject.instance().addMapLayer(zooming_layer)

    elif zoom == "query":
        with psycopg2.connect(**conn_parameters) as conn:
            with conn.cursor() as cursor:
                cursor.execute("SELECT ST_asText(ST_Envelope(geometry)) FROM osm.world_continents_boundaries as cna WHERE cna.continent = 'North America'")
                wkt_string = cursor.fetchone()[0]
        
        # zooming_layer = QgsVectorLayer(f"?query=SELECT ST_GeomFromText('{wkt_string}')", "extent", "virtual")
            zoom = QgsRectangle.fromWkt(wkt_string)

    elif zoom == "":
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

def change_raster_layer_color(layer):
    provider = layer.dataProvider()
    extent = layer.extent()
    stats = provider.bandStatistics(1, QgsRasterBandStats.All, extent, 0)

    min_value = stats.minimumValue
    max_value = stats.maximumValue

    # Define color ramp
    # color_ramp = [QgsColorRampShader.ColorRampItem(min_value, QColor(0, 0, 0, 0), f"Min ({min_value})"),QgsColorRampShader.ColorRampItem(max_value, QColor(255, 255, 255, 255), f"Max ({max_value})")]
    color_ramp = [QgsColorRampShader.ColorRampItem(min_value, QColor(0,255,0)),QgsColorRampShader.ColorRampItem(max_value, QColor(255,255,0))]
    # Define the color ramp
    color_ramp_shader = QgsColorRampShader()
    color_ramp_shader.setColorRampType(QgsColorRampShader.Interpolated) # QgsColorRampShader.Interpolated means 0
    color_ramp_shader.setColorRampItemList(color_ramp)

    raster_shader = QgsRasterShader()
    raster_shader.setRasterShaderFunction(QgsRasterShaderFunction(color_ramp_shader))

    # Apply the shader to the raster layer
    renderer = QgsSingleBandPseudoColorRenderer(layer.dataProvider(), 1, raster_shader)
    layer.setRenderer(renderer)
    
    # Refresh the layer to apply the changes
    # layer.triggerRepaint()
    return layer

def change_vector_layer_color(layer, field_name="values"):
    renderer = layer.renderer()
    renderer = QgsGraduatedSymbolRenderer(field_name).convertFromRenderer(renderer)

    # Define symbol categories (e.g., color and label for each category)
    categories = []
    # Replace these values and colors with your own
    categories.append(QgsRendererCategory([-99999, 0], QgsSymbol.defaultSymbol(layer.geometryType()).setColor(QColor("blue")), 'Category 1'))
    categories.append(QgsRendererCategory([1, 1000], QgsSymbol.defaultSymbol(layer.geometryType()).setColor(QColor("blue")), 'Category 2'))
    categories.append(QgsRendererCategory([1001, 4000], QgsSymbol.defaultSymbol(layer.geometryType()).setColor(QColor("blue")), 'Category 3'))

    # renderer.setSource(field_name)
    # renderer.setMode(QgsGraduatedSymbolRenderer.Custom)
    renderer.setClassAttribute(field_name)
    # renderer.updateClasses(categories)
    renderer.updateClasses(layer, nclasses=15)
    # color_ramp = QgsColorRamp().clone()
    # print(color_ramp.properties())
    renderer.updateColorRamp()
    # Apply the renderer to the layer
    layer.setRenderer(renderer)
    # Refresh the layer to see the changes
    layer.triggerRepaint()
    return layer
