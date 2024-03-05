import pytest
import sys
import time
import pandas as pd
import csv

from .test_pixel_utils import get_tif_pixel_value

sys.path.append("/home/nikola/4_north_america/GeoDataPump/scripts/")
from scripts.pipelines import model_data
import settings
config = settings.get_config()

GDAL_DATASET_PATH = config["NA_coastal_flooding_90m_tif"]

@pytest.fixture
def lat():
    lat = 30.52807498421037
    return lat

@pytest.fixture
def lon():
    lon = -89.30395396106162
    return lon


class TestCoastalFlood:
    def test_get_na_coastal_flooding_90(self, lat, lon):
        # Add more test cases as needed
        result = model_data.get_na_coastal_flooding_90(lat, lon)
        assert isinstance(result, (int, float)) or result is None

    def test_gdal_location_info(self, lat, lon):
        # Add more test cases as needed
        result = get_tif_pixel_value(lat, lon, GDAL_DATASET_PATH)
        assert isinstance(result, (int, float)) or result is None

    def test_integration(database_session):
        # Test the integration of both functions and CSV output
        from test_create_coordinates import generate_random_points, north_american_bounds

        fields = ['Latitude','Longitude','Database_Result', 'db_time', 'GDAL_Result', 'tif_time']
        csv_file_path = 'test_output.csv'
        points = generate_random_points(50, north_american_bounds)
        results = []

        for p in points:
            lat, lon = p

            start_time = time.time()
            result_tif = get_tif_pixel_value(lat, lon, GDAL_DATASET_PATH)
            fun_time_tif = round(time.time() - start_time, 2)

            start_time = time.time()
            result_db = model_data.get_na_coastal_flooding_90(lat, lon)
            fun_time_db = round(time.time() - start_time, 2)

            results.append((lat, lon, result_db, fun_time_db, result_tif, fun_time_tif))

        # export
        df = pd.DataFrame(results, columns=fields)
        df.to_csv(csv_file_path, index=False)

        # Check if the CSV file is created and has the expected structure
        with open(csv_file_path, 'r') as csvfile:
            lines = csvfile.readlines()
            assert len(lines) > 0
            assert lines[0].strip().split(",") == fields

            reader = csv.DictReader(csvfile)
            for row in reader:
                # Extract values from the row
                database_result = row['Database_Result']
                gdal_result = row['GDAL_Result']

                # Convert values to integers (assuming they are numeric)
                database_result = int(database_result) if database_result.isdigit() else None
                gdal_result = int(gdal_result) if gdal_result.isdigit() else None

                # Check the conditions
                if gdal_result == 0:
                    assert database_result is None
                else:
                    assert database_result == gdal_result - 1

        # Clean up: Remove the generated CSV file
        # import os
        # os.remove(csv_file_path)

class TestCogSpeed:


