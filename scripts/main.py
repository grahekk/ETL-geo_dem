import sys
import os
from pipelines import (
    pipeline_download_Canada_DEM,
    pipeline_download_USGS_DEM,
    pipeline_download_utils_soils,
    pipeline_transform_vrt_gdal,
    # pipeline_download_WMS_canada,
    pipeline_download_s3_global,
    pipeline_flows,
    gdal_build_vrt,
    # pipeline_transform_qgis_resample
)
import settings
from qgis_other import qgis_run_app, qgis_create_graphics, qgis_create_grid

conn_parameters = settings.get_conn_parameters()
schema = settings.get_schema()
config = settings.get_config()

def main():
    # # Execute download functions
    # pipeline_download_Canada_DEM.download_canadian_dem_flow()
    # pipeline_download_WMS_canada.download_canada_WMS_flow()
    # pipeline_download_USGS_DEM.download_us_dem_flow()
    # pipeline_download_utils_soils.download_us_soil_flow()
    # pipeline_download_utils_soils.download_canada_soil_flow()
    # pipeline_download_utils_soils.download_cec_lc_flow()
    # pipeline_download_utils_soils.download_world_countries_boundaries_flow()
    # pipeline_download_utils_soils.download_world_continent_boundaries_flow()
    # pipeline_download_utils_soils.download_ESA_global_dem_grid_flow()
    # pipeline_download_s3_global.download_ESA_worldcover_flow()
    # pipeline_download_s3_global.download_ESA_global_dem_90_flow()
    # pipeline_download_s3_global.download_ESA_global_dem_30_flow()
    
    # Add code for transformation and loading steps when ready.
    # pipeline_transform_vrt.usa_dem_vrt_ovr_flow()
    # pipeline_transform_vrt.canada_dem_vrt_ovr_flow()
    # pipeline_transform_vrt.canada_WMS_dem_vrt_ovr_flow()
    # pipeline_transform_vrt.NASA_global_dem_vrt_ovr_flow()
    # pipeline_flows.transform_ESA_global_dem_90_vrt_flow()
    # pipeline_flows.transform_ESA_NA_dem_30_vrt_flow()
    # pipeline_flows.transform_ESA_NA_dem_90_vrt_flow()
    # pipeline_flows.transform_ESA_world_cover_vrt_flow()
    # pipeline_flows.transform_aspect_flow()
    pipeline_flows.transform_geomorphon_flow()

    # pipeline_flows.load_countries_boundaries_flow()
    # pipeline_flows.load_continents_boundaries_flow()
    # pipeline_flows.load_ESA_dem_grid_flow()

    # print(pipeline_download_s3_global.filter_by_continents())

    # creating global dem90 only for usa
    # pipeline_transform_vrt.ESA_usa_vrt_flow()

    # NA geomorphon
    # dem_file = config["esa_usa_dem_90_vrt"]
    # output_file = config["NA_geomorphon"]
    # dem_file = "/mnt/volume-nbg1-1/shared/nikola/ESA_global_dem_90m/Copernicus_DSM_30_N45_00_E015_00_DEM.tif"
    # output_file = "/mnt/volume-nbg1-1/shared/nikola/NA_geomorphon_N45_E015"
    # pipeline_transform_resample.transform_geomorphon_qgis(dem_file, output_file)

    # pipeline_transform_vrt.transform_aspect_flow()

    # eu_slope = "/mnt/volume-nbg1-1/satellite/eu_slope/eu_slope_100m_3035.tif"
    # na_slope = config["NA_slope"]
    # na_aspect_categorized = config["NA_aspect_categorized"]
    dem_90 = config["esa_global_dem_90_vrt"]
    # esa_world_cover_vrt = config["esa_world_cover_vrt"]
    # qgis_run_app.grapical_map_flow(esa_world_cover_vrt)

    # pipeline_flows.load_aspect_flow()
    # pipeline_flows.load_slope_flow()

if __name__ == "__main__": 
    print("Main ETL script running")
    main()
