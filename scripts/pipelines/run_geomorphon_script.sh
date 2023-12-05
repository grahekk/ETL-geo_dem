#!/bin/bash

# Define input parameters
input_vrt="/mnt/volume-nbg1-1/shared/nikola/ESA_NA_dem_90m.vrt"
output_prefix="output_prefix"
output_format="GTiff"  # Choose your desired output format

# Set the GRASS GIS environment variables
export GISBASE=/usr/bin/grass
export PATH=$PATH:$GISBASE/bin:$GISBASE/scripts
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GISBASE/lib

# Set the location and mapset
location_name="/mnt/volume-nbg1-1/shared/nikola"
mapset_name="PERMANENT"

# Create a new GRASS GIS location (skip if the location already exists)
# grass -c $input_vrt $location_name

# Set the region to match the extent of the input virtual raster
g.region raster=input_vrt

# Import the virtual raster into GRASS GIS
r.in.gdal input=$input_vrt output=input_raster --overwrite

# Run the r.geomorphon algorithm
r.geomorphon elevation=input_raster forms=$output_prefix format=$output_format

# Export the result to a GeoTIFF file (optional, adjust output path as needed)
r.out.gdal input=$output_prefix\_landforms output=/mnt/volume-nbg1-1/shared/nikola/geomorphon_output_landforms.tif format=GTiff --overwrite

# Exit GRASS GIS
grass -q -e

# Print a message indicating successful completion
echo "GRASS GIS script completed successfully."
