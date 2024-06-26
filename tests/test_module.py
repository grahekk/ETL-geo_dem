import pytest
import os
import sys
import time
import pandas as pd
import csv

from .test_pixel_utils import get_tif_pixel_value
# from .test import measure_execution_time
from .create_testing_coordinates import generate_random_points, north_american_bounds, europe_bounds


from scripts.pipelines import model_data
from scripts import settings
config = settings.get_config()

test_output = config['test_output']

na_coastal_flooding_90m_tif = config["NA_coastal_flooding_90m_tif"]
na_coastal_flooding_30m_tif = config["NA_coastal_flooding_30m_tif"]
EU_coastal_flooding_10m_vrt = config["EU_coastal_flooding_10m_vrt"]

esa_world_cover_vrt = config["esa_world_cover_vrt"]
esa_world_cover_aws = config["esa_world_cover_aws_url"]

NA_aspect = config["NA_aspect"]
NA_aspect_table = config["NA_aspect_table"]
NA_geomorphon = config["NA_geomorphon"]

NA_slope = config["NA_slope"]
NA_slope_table = config["NA_slope_table"]
NA_slope_aws = config["NA_slope_aws"]

NA_tree_cover_density_tif = config["NA_tree_cover_density_tif"]
NA_tree_cover_density_cog = config["NA_tree_cover_density_cog"]
NA_tree_cover_density_aws = config["NA_tree_cover_density_aws"]

esa_global_dem_90_vrt = config["esa_global_dem_90_vrt"]
esa_na_dem_90_cog = config["esa_global_dem_90_cog"]
esa_global_dem_90_aws = config["esa_global_dem_90_aws_url"]

esa_na_dem_30_vrt = config["esa_na_dem_30_vrt"]
esa_na_dem_30_cog = config["esa_na_dem_30_cog"]
esa_na_dem_30_aws_url = config["esa_na_dem_30_aws_url"]

usgs_dem_vrt = config["usgs_dem_vrt"]
usgs_dem_cog = config["usgs_dem_cog"]
usgs_dem_10_aws_url = config["usgs_dem_10_aws_url"]

flood_jrc = config["flood_jrc_dir"]
rps = 10,20,50,100,200,500
flood_jrc_paths = []
flood_jrc_cog_paths = []

for r in rps:
    flood_jrc_rp = os.path.join(flood_jrc, f"floodMapGL_rp{r}y", f"floodMapGL_rp{r}y.tif")
    flood_jrc_paths.append(flood_jrc_rp)
    flood_jrc_rp = os.path.join(flood_jrc, f"floodMapGL_rp{r}y", f"floodMapGL_rp{r}y_COG.tif")
    flood_jrc_cog_paths.append(flood_jrc_rp)


def measure_execution_time(func):
    """
    Decorator function to measure the execution time of another function.

    Parameters:
    - func (callable): The function to be decorated.
    """
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        execution_time = round(time.time() - start_time, 2)
        return result, execution_time
    return wrapper

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
        result = get_tif_pixel_value(lat, lon, na_coastal_flooding_90m_tif)
        assert isinstance(result, (int, float)) or result is None

    def test_integration(self):
        # Test the integration of both functions and CSV output
        # from test_create_coordinates import generate_random_points, north_american_bounds

        fields = ['Latitude','Longitude','Database_Result', 'db_time', 'GDAL_Result', 'tif_time']
        points = generate_random_points(10, north_american_bounds)
        results = []

        # the test
        for lat, lon in points:

            result_tif, fun_time_tif = measure_execution_time(get_tif_pixel_value)(lat, lon, na_coastal_flooding_90m_tif)
            result_db, fun_time_db = measure_execution_time(model_data.get_na_coastal_flooding_90)(lat, lon)
            results.append((lat, lon, result_db, fun_time_db, result_tif, fun_time_tif))

        # export
        df = pd.DataFrame(results, columns=fields)
        df.to_csv(test_output, index=False)

        # Check if the CSV file is created and has the expected structure
        with open(test_output, 'r') as csvfile:
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
        # os.remove(test_output)
    

    def test_loaded_eu_10m(self):
        """
        Test eu coastal flooding on 10m dem after it has been loaded into database. 
        Compare it with tif files that were predecestors to 
        """

        fields = ['Latitude','Longitude','Database_Result', 'db_time', 'GDAL_Result', 'tif_time']
        points = generate_random_points(10, europe_bounds)
        results = []

        # the test
        for lat, lon in points:

            result_tif, fun_time_tif = measure_execution_time(get_tif_pixel_value)(lat, lon, EU_coastal_flooding_10m_vrt)
            result_db, fun_time_db = measure_execution_time(model_data.get_EU_coastal_flooding_10)(lat, lon)
            results.append((lat, lon, result_db, fun_time_db, result_tif, fun_time_tif))

        # export
        df = pd.DataFrame(results, columns=fields)
        df.to_csv(test_output, index=False)

        # Check if the CSV file is created and has the expected structure
        with open(test_output, 'r') as csvfile:
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
        # os.remove(test_output)
    

class TestSpeed:
    def first_test():
        pass

    def test_integration(self):
        # Test the integration of both functions and CSV output

        fields = ['lat', 'lon', 
          'aspect_result', 
          'fun_time_aspect', 
          'aspect_table_result', 'fun_time_aspect_table', 
          'aws_aspect_table_result', 'fun_time_aws_aspect_table',
          'slope_result', 
          'fun_time_slope', 
          'slope_table_result', 'fun_time_slope_table', 
        #   'aws_slope_table_result', 'fun_time_aws_slope_table',
          'geomorphon_result', 'fun_time_geomorphon',
          'geomorphon_table_result', 'fun_time_geomorphon_table',
          'tcd_result', 'fun_time_tcd',
          'na_dem_90_vrt_result', 'fun_time_na_dem_90_vrt',
          'na_dem_90_local_90_dem_result', 'fun_time_na_dem_90_cog',
          'na_dem_30_vrt_result', 'fun_time_na_dem_30_vrt',
          'na_dem_30_cog_result', 'fun_time_na_dem_30_cog',
          'usgs_dem_vrt_result', 'fun_time_usgs_dem_vrt',
          'usgs_dem_cog_result', 'fun_time_usgs_dem_cog',
        #   'cog_result', 'fun_time_cog', 
          'vrt_result', 'fun_time_vrt']
        points = generate_random_points(60, north_american_bounds)
        results = []

        # the test
        for lat, lon in points:

            # TODO: 
            # JRC flood
            # convert USGS to COG
            aspect_result, fun_time_aspect = measure_execution_time(get_tif_pixel_value)(lat, lon, NA_aspect)
            aspect_table_result, fun_time_aspect_table = measure_execution_time(model_data.get_na_aspect_proc)(lat, lon)
            aws_aspect_table_result, fun_time_aws_aspect_table = measure_execution_time(model_data.get_na_aspect_proc_psycopg)(lat, lon, "AWS")

            slope_result, fun_time_slope = measure_execution_time(get_tif_pixel_value)(lat, lon, NA_slope)
            slope_table_result, fun_time_slope_table = measure_execution_time(model_data.get_na_slope_proc)(lat, lon)
            # aws_slope_table_result, fun_time_aws_slope_table = measure_execution_time(model_data.get_na_slope_proc)(lat, lon, "AWS")

            geomorphon_result, fun_time_geomorphon = measure_execution_time(get_tif_pixel_value)(lat, lon, NA_geomorphon)
            geomorphon_table_result, fun_time_geomorphon_table = measure_execution_time(model_data.get_geomorphon_proc_psycopg)(lat, lon)

            tcd_result, fun_time_tcd = measure_execution_time(get_tif_pixel_value)(lat, lon, NA_tree_cover_density_tif)
            na_dem_90_vrt_result, fun_time_na_dem_90_vrt = measure_execution_time(get_tif_pixel_value)(lat, lon, esa_global_dem_90_vrt)
            na_dem_90_cog_result, fun_time_na_dem_90_cog = measure_execution_time(get_tif_pixel_value)(lat, lon, esa_na_dem_90_cog)
            na_dem_30_vrt_result, fun_time_na_dem_30_vrt = measure_execution_time(get_tif_pixel_value)(lat, lon, esa_na_dem_30_vrt)
            na_dem_30_cog_result, fun_time_na_dem_30_cog = measure_execution_time(get_tif_pixel_value)(lat, lon, esa_na_dem_30_cog)

            # add usgs COG
            usgs_dem_vrt_result, fun_time_usgs_dem_vrt = measure_execution_time(get_tif_pixel_value)(lat, lon, usgs_dem_vrt)
            usgs_dem_cog_result, fun_time_usgs_dem_cog = measure_execution_time(get_tif_pixel_value)(lat, lon, usgs_dem_cog)

            # cog_result, fun_time_cog = measure_execution_time(get_tif_pixel_value)(lat, lon, esa_world_cover_cog)
            vrt_result, fun_time_vrt = measure_execution_time(get_tif_pixel_value)(lat, lon, esa_world_cover_vrt)

            results.append((lat, lon, 
                            aspect_result, 
                            fun_time_aspect, 
                            aspect_table_result, fun_time_aspect_table,
                            aws_aspect_table_result, fun_time_aws_aspect_table,
                            slope_result, 
                            fun_time_slope, 
                            slope_table_result, fun_time_slope_table,
                            # aws_slope_table_result, fun_time_aws_slope_table,
                            geomorphon_result, 
                            fun_time_geomorphon,
                            geomorphon_table_result, fun_time_geomorphon_table,
                            tcd_result, fun_time_tcd,
                            na_dem_90_vrt_result, fun_time_na_dem_90_vrt,
                            na_dem_90_cog_result, fun_time_na_dem_90_cog,
                            na_dem_30_vrt_result, fun_time_na_dem_30_vrt,
                            na_dem_30_cog_result, fun_time_na_dem_30_cog,
                            usgs_dem_vrt_result, fun_time_usgs_dem_vrt,
                            usgs_dem_cog_result, fun_time_usgs_dem_cog,
                            # cog_result, 
                            # fun_time_cog, 
                            vrt_result, 
                            fun_time_vrt))

        # export
        df = pd.DataFrame(results, columns=fields)
        df.to_csv(test_output, index=False)

        # Check if the CSV file is created and has the expected structure
        # with open(test_output, 'r') as csvfile:
        #     lines = csvfile.readlines()
        #     assert len(lines) > 0

        #     reader = csv.DictReader(csvfile)
        #     for row in reader:
        #         # Extract values from the row
        #         cog_result = row['cog_result']
        #         vrt_result = row['vrt_result']

                # Convert values to integers (assuming they are numeric)
                # cog_result = int(cog_result) if cog_result.isdigit() else None
                # vrt_result = int(vrt_result) if vrt_result.isdigit() else None

                # Check the conditions
                # assert cog_result == vrt_result


class TestWorldCover:
    def first_test():
        pass

    def test_integration(self):
        # Test the integration of both functions and CSV output

        fields = ['Latitude','Longitude','vrt_result', 'vrt_time', 'cog_result', 'cog_time']
        points = generate_random_points(100, north_american_bounds)
        results = []

        # the test
        for lat, lon in points:

            # cog_result, fun_time_cog = measure_execution_time(get_tif_pixel_value)(lat, lon, esa_world_cover_cog)
            vrt_result, fun_time_vrt = measure_execution_time(get_tif_pixel_value)(lat, lon, esa_world_cover_vrt)
            results.append((lat, lon, 
                            vrt_result, 
                            fun_time_vrt 
                            # cog_result, 
                            # fun_time_cog,
                            ))

        # export
        df = pd.DataFrame(results, columns=fields)
        df.to_csv(test_output, index=False)

        # Check if the CSV file is created and has the expected structure
        with open(test_output, 'r') as csvfile:
            lines = csvfile.readlines()
            assert len(lines) > 0
            assert lines[0].strip().split(",") == fields

            reader = csv.DictReader(csvfile)
            for row in reader:
                # Extract values from the row
                cog_result = row['cog_result']
                vrt_result = row['vrt_result']

                # Convert values to integers (assuming they are numeric)
                cog_result = int(cog_result) if cog_result.isdigit() else None
                vrt_result = int(vrt_result) if vrt_result.isdigit() else None

                # Check the conditions
                assert cog_result == vrt_result


class TestJrc:
    def first_test():
        pass

    def test_jrc_integration(self):
        # Test the integration of both functions and CSV output

        fields = ['Latitude','Longitude','tif_result', 'tif_time', 'cog_result', 'cog_time']
        points = generate_random_points(100, europe_bounds)
        results = []

        # the test
        for lat, lon in points:

            cog_result, fun_time_cog = measure_execution_time(get_tif_pixel_value)(lat, lon, eu_flood_jrc_path)
            tif_result, fun_time_tif = measure_execution_time(get_tif_pixel_value)(lat, lon, eu_flood_jrc_path)
            results.append((lat, lon, tif_result, fun_time_tif, cog_result, fun_time_cog))

        # export
        df = pd.DataFrame(results, columns=fields)
        df.to_csv(test_output, index=False)

        # Check if the CSV file is created and has the expected structure
        with open(test_output, 'r') as csvfile:
            lines = csvfile.readlines()
            assert len(lines) > 0
            assert lines[0].strip().split(",") == fields

            reader = csv.DictReader(csvfile)
            for row in reader:
                # Extract values from the row
                cog_result = row['cog_result']
                tif_result = row['tif_result']

                # Check the conditions
                assert cog_result == tif_result


class TestTreeCover:
    def first_test():
        pass

    def test_tcd_integration(self):
        # Test the integration of both functions and CSV output

        fields = ['Latitude','Longitude','tif_result', 'tif_time', 'cog_result', 'cog_time', 'aws_result', "aws_time"]
        points = generate_random_points(100, north_american_bounds)
        results = []

        # the test
        for lat, lon in points:

            cog_result, fun_time_cog = measure_execution_time(get_tif_pixel_value)(lat, lon, NA_tree_cover_density_cog)
            tif_result, fun_time_tif = measure_execution_time(get_tif_pixel_value)(lat, lon, NA_tree_cover_density_tif)
            aws_result, fun_time_aws = measure_execution_time(get_tif_pixel_value)(lat, lon, NA_tree_cover_density_aws)
            results.append((lat, lon, tif_result, fun_time_tif, cog_result, fun_time_cog, aws_result, fun_time_aws))

        # export
        df = pd.DataFrame(results, columns=fields)
        df.to_csv(test_output, index=False)

        # Check if the CSV file is created and has the expected structure
        with open(test_output, 'r') as csvfile:
            lines = csvfile.readlines()
            assert len(lines) > 0
            assert lines[0].strip().split(",") == fields

            reader = csv.DictReader(csvfile)
            for row in reader:
                # Extract values from the row
                cog_result = row['cog_result']
                tif_result = row['tif_result']
                aws_result = row['aws_result']

                # Check the conditions
                assert cog_result == tif_result == aws_result


class TestAws:
    """
    Test files on aws s3 storage
    """

    def test_90m_dem(self):
         # Test the integration of both functions and CSV output

        fields = ['Latitude','Longitude','dem_90_result', 'dem_90_time', 'local_90_dem_result', 'local_90_dem_time']
        points = generate_random_points(20, north_american_bounds)
        results = []

        # the test
        for lat, lon in points:

            dem_90_result, fun_time_dem_90 = measure_execution_time(get_tif_pixel_value)(lat, lon, esa_global_dem_90_aws)
            local_90_dem_result, fun_local_90_dem_time = measure_execution_time(get_tif_pixel_value)(lat, lon, esa_na_dem_90_cog)
            results.append((lat, lon, dem_90_result, fun_time_dem_90, local_90_dem_result, fun_local_90_dem_time))

        # export
        df = pd.DataFrame(results, columns=fields)
        df.to_csv(test_output, index=False)

        # Check if the CSV file is created and has the expected structure
        with open(test_output, 'r') as csvfile:
            lines = csvfile.readlines()
            assert len(lines) > 0
            assert lines[0].strip().split(",") == fields

            reader = csv.DictReader(csvfile)
            for row in reader:
                # Extract values from the row
                dem_90_result = row['dem_90_result']
                local_90_dem_result = row['local_90_dem_result']

                # Check the conditions
                assert local_90_dem_result == dem_90_result

    def test_world_cover_aws(self):
         # Test the integration of both functions and CSV output

        fields = ['Latitude','Longitude','esa_world_cover_aws_result', 'world_cover_aws_time', 'local_world_cover_result', 'local_world_cover_time']
        points = generate_random_points(20, north_american_bounds)
        results = []

        # the test
        for lat, lon in points:

            esa_world_cover_aws_result, fun_time_world_cover_aws = measure_execution_time(get_tif_pixel_value)(lat, lon, esa_world_cover_aws)
            local_world_cover_result, fun_local_world_cover_time = measure_execution_time(get_tif_pixel_value)(lat, lon, esa_world_cover_vrt)
            results.append((lat, lon, esa_world_cover_aws_result, fun_time_world_cover_aws, local_world_cover_result, fun_local_world_cover_time))

        # export
        df = pd.DataFrame(results, columns=fields)
        df.to_csv(test_output, index=False)

        # Check if the CSV file is created and has the expected structure
        with open(test_output, 'r') as csvfile:
            lines = csvfile.readlines()
            assert len(lines) > 0
            assert lines[0].strip().split(",") == fields

            reader = csv.DictReader(csvfile)
            for row in reader:
                # Extract values from the row
                esa_world_cover_aws_result = row['esa_world_cover_aws_result']
                local_world_cover_result = row['local_world_cover_result']

                # Check the conditions
                assert local_world_cover_result == esa_world_cover_aws
    
    def test_usgs_na_dem(self):
         # Test the integration of both functions and CSV output

        fields = ['Latitude','Longitude','usgs_dem_10_aws_url_result', 'time_usgs_10', 'esa_na_dem_30_aws_result', 'time_esa_30_dem_aws']
        points = generate_random_points(20, north_american_bounds)
        results = []

        # the test
        for lat, lon in points:

            usgs_dem_10_aws_url_result, time_usgs_10 = measure_execution_time(get_tif_pixel_value)(lat, lon, usgs_dem_10_aws_url)
            esa_na_dem_30_aws_result, time_esa_30_dem_aws = measure_execution_time(get_tif_pixel_value)(lat, lon, esa_na_dem_30_aws_url)
            results.append((lat, lon, usgs_dem_10_aws_url_result, time_usgs_10, esa_na_dem_30_aws_result, time_esa_30_dem_aws))

        # export
        df = pd.DataFrame(results, columns=fields)
        df.to_csv(test_output, index=False)

        # Check if the CSV file is created and has the expected structure
        with open(test_output, 'r') as csvfile:
            lines = csvfile.readlines()
            assert len(lines) > 0
            assert lines[0].strip().split(",") == fields

            reader = csv.DictReader(csvfile)
            for row in reader:
                # Extract values from the row
                usgs_dem_10_aws_url_result = row['usgs_dem_10_aws_url_result']
                esa_na_dem_30_aws_result = row['esa_na_dem_30_aws_result']

                # Check the conditions
                # assert usgs_dem_10_aws_url_result == esa_na_dem_30_aws_result
    
    def test_slope_aws(self):
         # Test the integration of both functions and CSV output

        fields = ['Latitude','Longitude','NA_slope_aws_result', 'time_NA_slope_aws', 'NA_slope', 'time_NA_slope']
        points = generate_random_points(20, north_american_bounds)
        results = []

        # the test
        for lat, lon in points:

            NA_slope_aws_result, time_NA_slope_aws = measure_execution_time(get_tif_pixel_value)(lat, lon, NA_slope_aws)
            NA_slope_result, time_NA_slope = measure_execution_time(get_tif_pixel_value)(lat, lon, NA_slope)
            results.append((lat, lon, NA_slope_aws_result, time_NA_slope_aws, NA_slope_result, time_NA_slope))

        # export
        df = pd.DataFrame(results, columns=fields)
        df.to_csv(test_output, index=False)

        # Check if the CSV file is created and has the expected structure
        with open(test_output, 'r') as csvfile:
            lines = csvfile.readlines()
            assert len(lines) > 0
            assert lines[0].strip().split(",") == fields

            reader = csv.DictReader(csvfile)
            for row in reader:
                # Extract values from the row
                NA_slope_aws_result = row['NA_slope_aws_result']
                NA_slope_result = row['NA_slope']

                # Check the conditions
                assert NA_slope_aws_result == NA_slope_result


class TestEUCog:
    def test_speed(self):
        fields = ['Latitude','Longitude','tif_result', 'tif_time', 'cog_result', 'cog_time']
        points = generate_random_points(100, europe_bounds)
        results = []

        # the test
        for lat, lon in points:

            cog_result, fun_time_cog = measure_execution_time(get_tif_pixel_value)(lat, lon, eu_vrt_dem)
            tif_result, fun_time_tif = measure_execution_time(get_tif_pixel_value)(lat, lon, eu_cog_dem)
            results.append((lat, lon, tif_result, fun_time_tif, cog_result, fun_time_cog))

        # export
        df = pd.DataFrame(results, columns=fields)
        df.to_csv(test_output, index=False)

        # Check if the CSV file is created and has the expected structure
        with open(test_output, 'r') as csvfile:
            lines = csvfile.readlines()
            assert len(lines) > 0
            assert lines[0].strip().split(",") == fields

            reader = csv.DictReader(csvfile)
            for row in reader:
                # Extract values from the row
                cog_result = row['cog_result']
                tif_result = row['tif_result']

                # Check the conditions
                assert cog_result == tif_result