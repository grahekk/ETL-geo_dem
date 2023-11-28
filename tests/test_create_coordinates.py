import requests
import pandas as pd
import random
import geopandas as gpd
from shapely.geometry import Point
import os
import sys
import time
import statistics

column_names = ["lat", "lon", "city", "country"]

# Define the bounding box for testing points to be randomized in it (latitude and longitude range)
europe_bounds = {
    "min_lat": 35.0,   
    "max_lat": 71.0,   
    "min_lon": -24.0,  
    "max_lon": 45.0    
}

north_american_bounds = {
    "min_lat": 14.559032,  # Southernmost point (Santo Tomas, Chiapas, Mexico)
    "max_lat": 84.057411,  # Northernmost point (Kaktovik, Alaska)
    "min_lon": -168.101265,  # Westernmost point (Adak, Alaska)
    "max_lon": -34.633207  # Easternmost point (St. John's, Newfoundland, Canada)
}

def generate_random_points(num_points, bounds):
    """
    Generate random latitude and longitude points within specified bounds.

        This function generates a specified number of random latitude and longitude points
        within the provided bounds. The generated points are stored in a list and returned.

        Parameters:
        - num_points (int): The number of random points to generate.
        - bounds (dict): A dictionary containing the minimum and maximum latitude and longitude
        values as "min_lat," "max_lat," "min_lon," and "max_lon."
    """
    random_points = []
    for _ in range(num_points):
        latitude = random.uniform(bounds["min_lat"], bounds["max_lat"])
        longitude = random.uniform(bounds["min_lon"], bounds["max_lon"])
        random_points.append((latitude, longitude))
    return random_points


def get_city_country(point):
    """
    Retrieve city and country information for a latitude and longitude point.

    This function sends a request to the Nominatim API to retrieve city and country
    information for a given latitude and longitude point. The response data is parsed,
    and the city and country information is extracted and returned.

    Parameters:
    - point (tuple): A tuple containing latitude and longitude values.
    """
    url = f"https://nominatim.openstreetmap.org/reverse?lat={point[0]}&lon={point[1]}&format=json"
    response = requests.get(url)
    data = response.json()
    city = data.get("address", {}).get("city", "")
    country = data.get("address", {}).get("country", "")
    return city, country


def add_city_country(points):
    """
    Add city and country information to a list of latitude and longitude points.

    This function takes a list of latitude and longitude points and enriches it with
    corresponding city and country information obtained using the Nominatim API.
    The resulting list contains tuples with the original point followed by city and country.

    Parameters:
    - points (list of tuple): A list of latitude and longitude points
    """
    total_count = len(points)
    points_with_city_country = []
    fun_times = []
    for count, point in enumerate(points):

        start_time = time.time()
        city, country = get_city_country(point)
        points_with_city_country.append((*point, city, country))
        fun_time = time.time()-start_time
        fun_times.append(fun_time)

        # print out the task completeness
        if count % 10 == 0:
            mean_time = statistics.mean(fun_times)
            time_left = round(((total_count-count)*mean_time)/60,2)
        percentage = round(count / total_count * 100, 2)
        sys.stdout.write('\r')
        sys.stdout.write(f'add_city_country task completeness: {percentage:.2f}% time left: {time_left:.2f} min')
        
    print("\nTask add_city_country complete!")
    return points_with_city_country


def save_testing_points(points, output="test_coordinates.csv"):
    """
    Save latitude, longitude, city, and country data to csv or to a shapefile.

    This function takes a list of latitude and longitude points enriched with city and country
    information and saves it to a file with specified column names.

    Parameters:
    - points (list of tuple): A list of tuples containing latitude, longitude, city, and country data.
    - output (str): The filename for the output file.
    """
    df = pd.DataFrame(points)
    df.columns = column_names
    if output.lower().endswith("csv"):
        df.to_csv(output, index=False)
    elif output.lower().endswith("shp") or output.lower().endswith("gpkg"):
        gdf = gpd.GeoDataFrame(points, columns=["latitude", "longitude", "city", "country"],
                            geometry=[Point(lon, lat) for lat, lon, _, _ in points])
        gdf.to_file(output)


def read_coordinates(file_name:str):
    """
    Read coordiantes from existing files. 
    Function recognizes whether coords are saved into csv, shp or gpkg format.
    """
    if file_name.lower().endswith("csv"):
        df_test_coordinates = pd.read_csv(file_name)
    elif file_name.lower().endswith("shp") or file_name.lower().endswith("gpkg"):
        gdf = gpd.read_file(file_name)
        # try to use the explode method to separate MultiPoints into individual points 
        try:
            gdf_exploded = gdf.explode(index_parts = True)
            # gdf explode also creates duplicates so they need to be removed
            gdf_exploded.drop_duplicates()
        except:
            gdf_exploded = gdf
        latitudes = gdf_exploded['geometry'].y
        longitudes = gdf_exploded['geometry'].x
        try: 
            city = gdf_exploded['city']
            country = gdf_exploded['country']
            df_test_coordinates = pd.DataFrame({'lat': latitudes, 'lon': longitudes, "city": city, "country": country})
        except:
            df_test_coordinates = pd.DataFrame({'lat': latitudes, 'lon': longitudes})

    test_coordinates = [tuple(row) for row in df_test_coordinates.to_numpy()]
    return test_coordinates


# Utilize all the functions and create points from zero 
def create_testing_points():
    """
    Utilize all the functions from this module and create testing points with lat, lon coordinates from the scratch
    """
    random_points = generate_random_points(1000, north_american_bounds)
    points = add_city_country(random_points)
    save_testing_points(points)
    return print("testing points created and saved into file")

# or use existing coordinates
def enrich_existing_shp(file_path:str):
    """
    Add city and country attributes to existing file with test coordinates.

    Args:
        file_path (str): path to file that needs city and country appending
    Returns:
        filename (str): a name of new file (enriched with cty, country)
    """
    test_coordinates = read_coordinates(file_path)
    points = add_city_country(test_coordinates)
    filename = os.path.splitext(file_path)[0]
    new_file = f"{filename}_enrich.shp"
    save_testing_points(points, f"{filename}_enrich.shp")
    # save_testing_points(points)
    return new_file

if __name__ == "__main__":
    # create_testing_points()
    enrich_existing_shp("/home/nikola/continents_grid_clipped.shp")

