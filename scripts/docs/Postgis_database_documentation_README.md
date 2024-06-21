# Postgres / PostGIS climatig Database - tables and functions in use

## General purpose

General purpose of this document is to show which tables and functions are used in the PG database. They are shown below, at the bottom.
To update this readme file or to get the tables in csv you can execute the code.

You can execute python code which executes SQL queries in our PostgreSQL database and display the results in tables mentioned down below in the middle section of this doc.

Part of the code below is in the separate `update_README_tables.py` script that updates the tables at the end of the document.

## Prerequisites for executing python code for fetching the tables and functions

Before running the Python code, make sure you have the following prerequisites installed:

- PostgreSQL Python library (`psycopg2`)
- Pandas, os
- Library for reading env files for db connection parameters (`dotenv`)

You can also show results interactively in the jupiter notebook, without having to run the `update_README_tables.py` script to update the data. Just run the code below.

## Python Code for SQL queries

```python
import psycopg2
import pandas as pd
import os
from dotenv import load_dotenv


# Define the PostgreSQL connection parameters
load_dotenv()

host_name = os.getenv('host_name')
port = os.getenv('port')
database_name = os.getenv('database_name')
username = os.getenv('username')
password = os.getenv('password')
schema = 'osm'

db_params = {
    "host": host_name,
    "database": database_name,
    "user": username,
    "password": password,
}

# Function to execute a SQL query and return the results as a Pandas DataFrame
def execute_sql_query(query):
    connection = psycopg2.connect(**db_params)
    cursor = connection.cursor()
    cursor.execute(query)
    data = cursor.fetchall()
    columns = [desc[0] for desc in cursor.description]
    cursor.close()
    connection.close()
    return pd.DataFrame(data, columns=columns)

# Execute the first query to get matching functions
query_functions = """
SELECT 
    routine_name AS function_name,
    routines.data_type AS return_type,
    array_agg(parameters.data_type) AS argument_types
FROM 
    information_schema.routines
LEFT JOIN
    information_schema.parameters
ON 
    routines.specific_name = parameters.specific_name
WHERE 
	routine_schema LIKE 'osm'
GROUP BY
    routine_name, routines.data_type
ORDER BY
    routine_name;
"""

# Execute the query and display the result as a table
result_functions = execute_sql_query(query_functions)


query_tables = """
SELECT 
   table_name, 
   ARRAY_AGG(ROW(column_name, data_type)::text) AS columns_and_data_types
FROM 
   information_schema.columns
WHERE 
   table_schema = 'osm'
GROUP BY table_name;
"""

result_tables = execute_sql_query(query_tables)

# these results can be then saved to csv or/and directly populate the tables in this markdown document

```

Tables that were produced in the queries, mentioned above in the python code, are then insertend into this markdown document using `update_README_tables.py`.

# Postgres/PostGIS sctructures used in the app
## Tables


<!-- START TABLES TABLE -->
| table_name                              | columns_and_data_types                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
|:----------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| world_lakes                             | ['(geom,USER-DEFINED)', '(fid,integer)', '(osm_id,bigint)', '(waterway,"character varying")', '(landuse,"character varying")', '(amenity,"character varying")', '(leisure,"character varying")', '(layer,"character varying")', '(path,"character varying")', '(name,"character varying")', '(natural,"character varying")']                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| europe_copernicus_clc2018_legend        | ['(grid_code,integer)', '(clc_code,"character varying")', '(label_corine,"character varying")', '(rgb,"character varying")', '(climatig_code,"character varying")', '(label_climatig,"character varying")']                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| europe_copernicus_forest_types_10m_3035 | ['(rid,integer)', '(rast,USER-DEFINED)']                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| world_rivers_lines                      | ['(fid,integer)', '(osm_id,bigint)', '(geom,USER-DEFINED)', '(barrier,"character varying")', '(name,"character varying")', '(waterway,"character varying")', '(layer,"character varying")', '(path,"character varying")']                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| world_rivers_polygons                   | ['(geom,USER-DEFINED)', '(fid,integer)', '(osm_id,bigint)', '(name,"character varying")', '(waterway,"character varying")', '(landuse,"character varying")', '(amenity,"character varying")', '(leisure,"character varying")', '(natural,"character varying")', '(layer,"character varying")', '(path,"character varying")']                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| world_rivers_difference                 | ['(fid,integer)', '(osm_id,"double precision")', '(geom,USER-DEFINED)', '(barrier,"character varying")', '(name,"character varying")', '(waterway,"character varying")', '(layer,"character varying")', '(path,"character varying")']                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| europe_copernicus_clc2018_100m          | ['(objectid,integer)', '(area_ha,"double precision")', '(geom,USER-DEFINED)', '(CLC_CODE,integer)', '(code_18,"character varying")', '(remark,"character varying")', '(id,"character varying")']                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| na_aspect_90m_4326                      | ['(rid,integer)', '(rast,USER-DEFINED)']                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| europe_coastal_flood_100m               | ['(fid,integer)', '(id,integer)', '(max,"double precision")', '(sea_level,"double precision")', '(geom,USER-DEFINED)', '(layer,"character varying")', '(path,"character varying")']                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| europe_coastal_flooding_90m_4326        | ['(gid,integer)', '(level,"double precision")', '(geom,USER-DEFINED)']                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| world_coastlines                        | ['(ogc_fid,integer)', '(fid,numeric)', '(geom,USER-DEFINED)']                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| world_landuse                           | ['(fid,integer)', '(geom,USER-DEFINED)', '(osm_id,"character varying")', '(man_made,"character varying")', '(barrier,"character varying")', '(landuse,"character varying")', '(tourism,"character varying")', '(amenity,"character varying")', '(natural,"character varying")', '(leisure,"character varying")', '(place,"character varying")', '(military,"character varying")', '(aeroway,"character varying")', '(name,"character varying")', '(highway,"character varying")', '(osm_way_id,"character varying")', '(type,"character varying")', '(admin_level,"character varying")', '(boundary,"character varying")', '(building,"character varying")', '(craft,"character varying")', '(geological,"character varying")', '(historic,"character varying")', '(land_area,"character varying")', '(office,"character varying")', '(shop,"character varying")', '(sport,"character varying")', '(other_tags,"character varying")', '(layer,"character varying")', '(path,"character varying")'] |
| world_seas                              | ['(fid,integer)', '(x,integer)', '(y,integer)', '(geom,USER-DEFINED)']                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| na_coastal_flooding_90m_4326            | ['(gid,integer)', '(level,"double precision")', '(geom,USER-DEFINED)']                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| europe_coastal_flood_10m                | ['(fid,integer)', '(sea_level,"double precision")', '(geom,USER-DEFINED)']                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| na_geomorphon_90m_4326                  | ['(rid,integer)', '(rast,USER-DEFINED)']                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
<!-- END TABLES TABLE -->

## Postgis functions used on mentioned tables

<!-- START FUNCTIONS TABLE -->
| function_name                      | return_type   | argument_types                                                                                                              |
|:-----------------------------------|:--------------|:----------------------------------------------------------------------------------------------------------------------------|
| calculate_coastal_flood            | record        | {numeric,"double precision",numeric,text,text,integer,integer,numeric}                                                      |
| calculate_distance_func            | numeric       | {numeric,text,text,integer,numeric,numeric,numeric,numeric,numeric}                                                         |
| calculate_distance_func_3035       | numeric       | {numeric,numeric,numeric,integer,text,text,numeric,numeric,numeric}                                                         |
| calculate_forest_distance_raster   | numeric       | {integer,numeric,numeric,numeric,numeric,numeric,numeric,text,text,integer}                                                 |
| calculate_get_pixel_value          | numeric       | {text,numeric,numeric,integer,integer,text}                                                                                 |
| calculate_is_within                | integer       | {numeric,integer,text,text,numeric}                                                                                         |
| calculate_landuse_area_percentage  | record        | {"double precision","character varying",numeric,text,text,integer,"double precision","double precision","double precision"} |
| get_closest_point_coordinates      | record        | {numeric,numeric,numeric,integer,text,text,numeric,numeric,numeric,"double precision","double precision"}                   |
| get_closest_point_coordinates_3035 | record        | {numeric,"double precision","double precision",numeric,numeric,numeric,numeric,text,text,integer,numeric}                   |
<!-- END FUNCTIONS TABLE -->