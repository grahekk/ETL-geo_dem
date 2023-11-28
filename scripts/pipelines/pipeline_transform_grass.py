import sys
import os
from grass_session import Session
import grass.script as gscript
import grass.script.setup as gsetup
from grass.script import core as gcore
from grass.pygrass.modules.shortcuts import general as g
from grass.pygrass.modules.shortcuts import raster as r


sys.path.append("/home/nikola/4_north_america/GeoDataPump/scripts/")
import settings
config = settings.get_config()

def run_geomorphon(input_map, output_map, search_radius=90, skip_flat=0):
    """
    Run the r.geomorphon algorithm in GRASS GIS.

    Parameters:
    - input_map (str): Name of the input elevation map.
    - output_map (str): Name of the output geomorphon map.
    - search_radius (int): Search radius for the algorithm.
    - skip_flat (bool): Skip flat areas when computing landforms.

    Returns:
    - None
    """
    gcore.run_command('r.geomorphon', elevation=input_map, forms=output_map,
                      search=search_radius)
    

if __name__ == '__main__':
    # Set your GRASS GIS database and mapset
    grassdata = config["esa_na_dem_90_vrt"]
    gisdb = config["general_dir"]
    tmp_location = config["tmp_dir"]
    # mapset = grassdata

    input_vrt = grassdata
    output_prefix = "Geomorphon_output_grass"
    output_format = "GTiff"

    # # Initialize GRASS GIS session
    # session = gsetup.init(grassdata, location, mapset)
    # print(gs.read_command("g.gisenv"))
    # grass.init(grass_db, mapset)

    with Session(gisdb=gisdb, location=gisdb, create_opts="EPSG:4326", mapset="PERMANENT") as PERMANENT:
        # run something in PERMANENT mapset:
        print(gcore.parse_command("g.gisenv", flags="s"))

        # Specify input and output map names
        input_map = config["esa_na_dem_90_vrt"]
        input_map = os.path.basename(input_map)
        output_map = config["NA_geomorphon_tmp"]
        output_map = "esa_na_dem_90_vrt_geomorphon.tif"

        # Specify algorithm parameters (adjust as needed)
        search_radius = 90
        skip_flat = True

        # Run the geomorphon algorithm
        # run_geomorphon(input_map, output_map, search_radius, skip_flat)

         # Import the virtual raster into GRASS GIS
        # gscript.run_command('r.in.gdal', input=input_vrt, output='input_raster', overwrite=True)
        r.in_gdal(input=input_map, output = output_map)

        # Set the region to match the extent of the input raster
        gscript.run_command('g.region', raster='input_raster')

        # Run the r.geomorphon algorithm
        gscript.run_command('r.geomorphon', elevation='input_raster', forms=output_prefix, format=output_format)

        # Export the result to a GeoTIFF file (optional, adjust output path as needed)
        gscript.run_command('r.out.gdal', input=output_prefix + '_landforms', output='/path/to/your/output_landforms.tif', format='GTiff', overwrite=True)

        # Print a message indicating successful completion
        print("GRASS GIS script completed successfully.")


    # Exit GRASS GIS session
    gscript.run_command('g.remove', flags='f', type='raster', name=output_map)
    gscript.run_command('g.remove', flags='f', type='raster', name=f'{output_map}_forms')
    gscript.run_command('g.remove', flags='f', type='raster', name=f'{output_map}_tcf')

    PERMANENT.close()
    gscript.exit()
