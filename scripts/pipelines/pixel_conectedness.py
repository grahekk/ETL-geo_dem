import geopandas as gpd
import scipy
import numpy as np
import rasterio
from rasterio.features import geometry_mask
from scipy import ndimage
from skimage import measure
from shapely.geometry import shape, LineString
from shapely.geometry import shape

raster_path = "/tmp/tmpwhycykz2/Copernicus_DSM_30_N45_00_E014_00_DEM_sea_level_0.tif"
output_tif_path = "/tmp/tmpwhycykz2/Intersecting_pixels_for_level_0.tif"
vector_path = "/tmp/tmpwhycykz2/Copernicus_DSM_30_N45_00_E014_00_DEM_polygonized_0_queried_coastlines.shp"

# extract connected pixels

def extract_intersecting_groups(main_raster, mask_raster):
    """
    Extract pixel groups from the main raster based on their intersection with a given mask raster.

    Parameters:
    - main_raster (numpy.ndarray): The main raster containing labeled connected components.
    - mask_raster (numpy.ndarray): The rasterized line mask used for intersection.

    Returns:
    - numpy.ndarray: A new raster where pixel values are retained only for groups that intersect with the mask.
    
    Example:
    >>> main_raster = np.array([[1, 1, 0, 0],
                                [1, 0, 2, 0],
                                [0, 0, 2, 2],
                                [3, 3, 3, 0]])
    >>> mask_raster = np.array([[0, 1, 0, 0],
                                [0, 1, 1, 0],
                                [0, 0, 0, 0],
                                [0, 0, 1, 0]])
    >>> result = extract_intersecting_groups(main_raster, mask_raster)
    >>> print(result)
    array([[0, 0, 0, 0],
           [1, 0, 2, 0],
           [0, 0, 2, 2],
           [0, 0, 0, 0]])

    Note:
    - The function uses skimage.measure.label to identify connected components in the main raster.
    - It iterates through labeled regions and checks for intersection with the provided mask raster.
    - The resulting raster retains only pixel values for groups that intersect with the mask.
    """

    # label connected components
    labeled_image = measure.label(main_raster)

    # Initialize a list to store intersecting groups
    intersecting_groups = []

    # Iterate through labeled regions and check for intersection
    for region in measure.regionprops(labeled_image):
        label_value = region.label
        region_mask = (labeled_image == label_value).astype(np.uint8)

        # Check if the region intersects with the rasterized line
        if np.any(region_mask * mask_raster):
            intersecting_groups.append(label_value)

    result_raster = np.where(np.isin(labeled_image, intersecting_groups), main_raster, 0)

    return result_raster

def get_intersection_between_raster_and_geometry(raster_path, output_tif_path, geometry_path):
    """
    Extract pixel groups from a raster based on their intersection with a vector geometry.

    Parameters:
    - raster_path (str): The file path to the input raster TIFF file.
    - output_tif_path (str): The file path to save the output GeoTIFF file.
    - geometry_path (str): The file path to the vector geometry (e.g., Shapefile, GeoJSON).

    Returns:
    None

    Example:
    >>> raster_path = "path/to/your/raster.tif"
    >>> output_tif_path = "path/to/your/output.tif"
    >>> geometry_path = "path/to/your/geometry.shp"
    >>> get_intersection_between_raster_and_geometry(raster_path, output_tif_path, geometry_path)

    Note:
    - The function reads the raster data and vector geometry using rasterio and geopandas, respectively.
    - It creates a mask based on the vector geometry using rasterio.features.geometry_mask.
    - The extract_intersecting_groups function is then used to obtain pixel groups that intersect with the mask.
    - The result is saved to a new GeoTIFF file with the specified output path.

    """
    with rasterio.open(raster_path) as src:
        # Read geometry from the vector file
        gdf = gpd.read_file(geometry_path)
        # Create a mask based on the geometry
        mask = geometry_mask(gdf.geometry, out_shape=src.shape, transform=src.transform, invert=True)
        raster_data = src.read(1)

        result = extract_intersecting_groups(raster_data, mask)

        # Create the output GeoTIFF file
        with rasterio.open(output_tif_path, 'w', driver='GTiff', 
                           width=src.width, height=src.height,
                           count=1, dtype=np.int8,
                           crs=src.crs, transform=src.transform) as dst:
            # Write the labeled data to the output raster
            dst.write(result, 1)

# Example usage:
get_intersection_between_raster_and_geometry(raster_path, output_tif_path, vector_path)
