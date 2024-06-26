import sys
import os
from pipelines import (
    # pipeline_download_Canada_DEM,
    # pipeline_download_USGS_DEM,
    # pipeline_download_utils_soils,
    pipeline_transform_vrt_gdal,
    # pipeline_download_WMS_canada,
    pipeline_download_s3_global,
    pipeline_flows,
    # pipeline_transform_qgis_resample,
    pipeline_transform_sea_level,
    model_data,
    pixel_conectedness
)
import settings
from qgis_other import qgis_run_app, qgis_create_graphics, qgis_create_grid
from docs import compile_json_metadata
# from line_profiler import LineProfiler
from tqdm import tqdm
from multiprocessing import Pool

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
    # pipeline_flows.download_tree_cover_density_flow()

    # Transform functions
    # pipeline_transform_vrt.usa_dem_vrt_ovr_flow()
    # pipeline_transform_vrt.canada_dem_vrt_ovr_flow()
    # pipeline_transform_vrt.canada_WMS_dem_vrt_ovr_flow()
    # pipeline_transform_vrt.NASA_global_dem_vrt_ovr_flow()
    # pipeline_flows.transform_ESA_global_dem_90_vrt_flow()
    # pipeline_flows.transform_ESA_NA_dem_30_vrt_flow()
    # pipeline_flows.transform_ESA_NA_dem_90_vrt_flow()
    # pipeline_flows.transform_esa_usa_dem_vrt_flow()
    # pipeline_flows.transform_ESA_world_cover_vrt_flow()
    # pipeline_flows.transform_aspect_flow()
    # pipeline_flows.transform_geomorphon_flow()
    # pipeline_flows.transform_geomorphon_chunky()
    # pipeline_flows.transform_geomorphon_multiprocess_flow()
    # pipeline_flows.transform_geomorphon_vrt_flow()
    # pipeline_flows.transform_tree_cover_density_vrt_flow()
    # pipeline_flows.eu_resample_1km_flow()
    # pipeline_flows.coastal_flooding_30m_NA_flow()
    # pipeline_flows.coastal_flooding_30m_NA_multiprocessing_flow()
    # pipeline_flows.coastal_flooding_10m_USA_flow()
    # pipeline_flows.coastal_flooding_10m_USA_multiprocessing_flow()
    # pipeline_flows.load_coastal_flooding_30m_to_aws_flow()
    # pipeline_flows.coastal_flooding_10m_Europe_flow()
    # pipeline_flows.load_coastal_flooding_10m_EU_to_mnt_flow()
    # pipeline_flows.load_coastal_flooding_10m_EU_to_aws_flow()
    pipeline_flows.load_coastal_flooding_10m_USA_to_aws_flow()

    # pipeline_flows.transform_coastal_flooding_vrt()
    # pipeline_flows.coastal_flooding_flow()

    # pipeline_flows.tree_cover_filter_flow()
    # pipeline_flows.shrubland_filter_flow()
    # pipeline_flows.grassland_filter_flow()

    # pipeline_flows.hawaii_repair_flow()
    # pipeline_flows.load_countries_boundaries_flow()
    # pipeline_flows.load_continents_boundaries_flow()
    # pipeline_flows.load_ESA_dem_grid_flow()
    # pipeline_flows.load_to_mnt_flow()
    # pipeline_flows.load_to_aws_flow()


    # pipeline_flows.coastal_flooding_90m_Europe_flow()
    # pipeline_flows.coastal_flooding_90m_pixels_flow()
    # pipeline_flows.coastal_flooding_90m_pixels_multiprocessing_flow()
    
    # pipeline_flows.transform_ESA_world_cover_vrt_flow()
    # compile_json_metadata.main()

if __name__ == "__main__": 
    print("Main ETL script running")
    main()
