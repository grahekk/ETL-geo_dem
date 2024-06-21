from sqlalchemy import create_engine
from sqlalchemy import Column, Integer, String, MetaData, Table
from geoalchemy2 import Geometry, Raster
from geoalchemy2.functions import ST_Intersects, ST_SetSRID, ST_MakePoint
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import declarative_base

import sys
from sqlalchemy import func
import psycopg2

sys.path.append("/home/nikola/4_north_america/GeoDataPump/scripts/")
import settings
database_url = settings.get_database_url()
aws_database_url = settings.get_aws_database_url()
schema = settings.get_schema()
config = settings.get_config()
conn_parameters = settings.get_conn_parameters()
aws_conn_parameters = settings.get_aws_conn_parameters()

aws_engine = create_engine(aws_database_url)
engine = create_engine(database_url)
metadata_obj = MetaData(schema=schema)

Base = declarative_base()

world_countries_boundaries = Table("world_countries_boundaries", 
                        metadata_obj,
                        Column("geom", Geometry('POLYGON')),
                        autoload_with = engine)

world_continents_boundaries = Table("world_continents_boundaries", 
                        metadata_obj,
                        Column("geometry", Geometry('POLYGON')),
                        autoload_with = engine)

esa_grid = Table("esa_global_dem_grid", 
                        metadata_obj,
                        Column("geometry", Geometry('POLYGON')),
                        autoload_with = engine)

na_coastal_flooding_90 = Table(config["NA_coastal_flooding_90m_table"], 
                        metadata_obj,
                        Column("geom", Geometry('POLYGON')),
                        autoload_with = engine)

na_coastal_flooding_30 = Table(config["NA_coastal_flooding_30m_table"], 
                        metadata_obj,
                        Column("geom", Geometry('POLYGON')),
                        autoload_with = engine)

eu_coastal_flooding_10m = Table(config["EU_coastal_flooding_10m_table"], 
                        metadata_obj,
                        Column("geom", Geometry('MULTIPOLYGON')),
                        autoload_with = engine)


Base.metadata.create_all(bind=engine)
Session = sessionmaker(bind=engine)
AwsSession = sessionmaker(bind=aws_engine)

def get_countries_boundaries(country = "United States"):
    """
    function to return boundary/geometry of a specified country

    Parameters:
        country (str): country name to be queried from database
    """
    # session is for interacting with db
    session = Session()
    country_boundaries = session.query(world_countries_boundaries).filter_by(country=country).first()

    for boundary in country_boundaries:
        print(f"ID: {boundary.ogc_fid}, Shapename: {boundary.country}, Geometry: {boundary.geom}")
        geometry = boundary.geom
        return_data = geometry
    # Close the session
    session.close()
    return return_data


def get_geocellid(c_name = "United States", by = "country", extent = (-25, 35, 35, 72)):
    """
    Retrieve a list of geocellids based on specified criteria (extent or c_name) from the ESA grid in a PostgreSQL database.

    Parameters:
    - c_name (str): Name of the country or continent to be queried from the database, depending on the 'by' argument.
    - by (str): Specify the category by which geocellids are filtered ('country', 'continent', or 'extent').
    - extent (tuple): Tuple specifying the bounding box extent (lon_min, lon_max, lat_min, lat_max) when 'by' is set to 'extent'.

    Returns:
    list: A list of geocellids based on the specified criteria.

    Example:
    ```
    # Retrieve geocellids for the United States by country name
    geocellids_by_country = get_geocellid(c_name="United States", by="country")

    # Retrieve geocellids for a specified extent, ex. Europe
    geocellids_by_extent = get_geocellid(by="extent", extent=(-25, 35, 35, 72))
    ```

    Notes:
    - The function uses GeoAlchemy for spatial queries on the ESA grid.
    - The 'by' parameter can be set to 'country', 'continent', or 'extent'.
    - When 'by' is set to 'country' or 'continent', the 'c_name' parameter specifies the country or continent name.
    - When 'by' is set to 'extent', the 'extent' parameter specifies the bounding box extent.
    - The returned list contains geocellids based on the specified criteria.
    - The session is automatically closed after retrieving the data.
    """
    assert by in ("country", "continent", "extent")
    assert len(extent) == 4

    # session is for interacting with db
    session = Session()

    if by == "country":
        boundary = world_countries_boundaries
        geocells = session.query(esa_grid).join(boundary, esa_grid.c.geometry.intersects(boundary.c.geom)).filter(boundary.c.country == c_name).all()
        return_data = [i.geocellid for i in geocells]

    elif by == "extent":
        lon_min, lon_max, lat_min, lat_max = extent
        bounding_box = func.ST_MakeEnvelope(lon_min, lat_min, lon_max, lat_max, 4326)
        geocells = session.query(esa_grid).filter(func.ST_Intersects(esa_grid.c.geometry, bounding_box)).all()
        return_data = [i.geocellid for i in geocells]

    elif by == "continent":
        boundary = world_continents_boundaries
        geocells = session.query(esa_grid).join(boundary, esa_grid.c.geometry.intersects(boundary.c.geometry)).filter(boundary.c.continent == c_name).all()
        return_data = [i.geocellid for i in geocells]
        
    # Close the session
    session.close()
    return return_data


def get_product_name(geocell_ids_list:list):
    """
    Retrieve product names associated with given Geocell IDs from the ESA grid.

    Parameters:
    - geocell_ids_list (list): List of Geocell IDs for which product names are to be retrieved.

    Returns:
    list: A list of product names corresponding to the provided Geocell IDs.

    Example:
    ```
    geocell_ids_list = ['N00E006', 'N00E007', 'N00E008']
    product_names = get_product_name(geocell_ids_list)
    print(product_names)
    ```

    Notes:
    - The function utilizes SQLAlchemy for database interaction.
    - Ensure that the `esa_grid` table is properly defined in your SQLAlchemy setup.
    - The returned list contains product names associated with the provided Geocell IDs.
    - The session is automatically closed after retrieving the data.
    """
    assert len(geocell_ids_list)>0 and len(geocell_ids_list[0])>5

    session = Session()
    product_names = session.query(esa_grid).filter(esa_grid.c.geocellid in geocell_ids_list).all()
    return_data = [i.product30 for i in product_names]

    session.close()
    return return_data

def get_na_coastal_flooding_90(lat:float, lon:float):
    """
    Query the information from `na_coastal_flooding_90` table by given point

    Arguments:
        lat(float): latitude
        lon(float): longitude 
    """
    session = Session()
    point = func.ST_SetSrid(func.ST_MakePoint(lon, lat), 4326)
    sea_level = session.query(na_coastal_flooding_90).filter(func.ST_Intersects(na_coastal_flooding_90.c.geom, point)).first()
    session.close()
    try:
        result = sea_level.sea_level
    except:
        result = None
    return result

def get_EU_coastal_flooding_10(lat:float, lon:float):
    """
    Query the information from `eu_coastal_flooding_10m_4326` table by given point

    Arguments:
        lat(float): latitude
        lon(float): longitude 
    """
    session = Session()

    # Create the point geometry

    # Build the query
    # query = select([data_table]).where(ST_Intersects(data_table.c.geom, point))

    # Execute the query and fetch results
    # results = session.execute(query).fetchall()

    point = ST_SetSRID(ST_MakePoint(lon, lat), 4326) # either this way, without use of func.
    # point = func.ST_SetSrid(func.ST_MakePoint(lon, lat), 4326)
    sea_level = session.query(eu_coastal_flooding_10m).where(func.ST_Intersects(eu_coastal_flooding_10m.c.geom, point)).first()
    session.close()
    try:
        result = sea_level.sea_level
    except:
        result = None
    return result


def get_na_slope_proc(lat:float, lon:float, base="local"):
    """
    Query info about slope by using stored procedure inside pgdb.
    """
    if base == "AWS":
        session = AwsSession()
    elif base == "local":
        session = Session()

    session = Session()
    slope_value = session.query(func.data.calculate_get_pixel_value(lon, lat, 4326, 4326, schema, config["NA_slope_table"])).first()
    result = slope_value[0]
    session.close()
    return result

def get_geomorphon(lat:float, lon:float, base="AWS"):
    """
    Query info about slope by using stored procedure inside pgdb.
    """
    if base == "AWS":
        session = AwsSession()
    elif base == "local":
        session = Session()

    session = Session()
    slope_value = session.query(func.data.calculate_get_pixel_value(lon, lat, 4326, 4326, schema, config["NA_geomorphon_table"])).first()
    result = slope_value[0]
    session.close()
    return result

def get_na_aspect_proc(lat:float, lon:float, base="local"):
    """
    Query info about aspect by using stored procedure inside pgdb.
    """
    if base == "AWS":
        session = AwsSession()
    elif base == "local":
        session = Session()

    aspect_value = session.query(func.data.calculate_get_pixel_value(lon, lat, 4326, 4326, schema, config["NA_aspect_table"])).first()
    result = aspect_value[0]
    session.close()
    return result

def get_na_aspect_proc_psycopg(lat:float, lon:float, base="local"):
    """
    Query info about aspect by using stored procedure inside pgdb.
    """
    with psycopg2.connect(**aws_conn_parameters) as conn:
        with conn.cursor() as curs:
            # sql = f'SELECT * FROM data.calculate_get_pixel_value({lon}, {lat}, {4326}, {4326}, {schema}, {config["NA_aspect_table"]})'
            curs.callproc(f'{schema}.calculate_get_pixel_value', [lon, lat, 4326, 4326, schema, config["NA_aspect_table"]])
            res = curs.fetchone()[0]
            # aspect_value = session.query(func.data.calculate_get_pixel_value(lon, lat, 4326, 4326, schema, config["NA_aspect_table"])).first()

    # leaving contexts doesn't close the connection
    conn.close()

    result = res
    return result

def get_geomorphon_proc_psycopg(lat:float, lon:float):
    """
    Query info about aspect by using stored procedure inside pgdb.
    """
    with psycopg2.connect(**aws_conn_parameters) as conn:
        with conn.cursor() as curs:
            # sql = f'SELECT * FROM data.calculate_get_pixel_value({lon}, {lat}, {4326}, {4326}, {schema}, {config["NA_aspect_table"]})'
            curs.callproc(f'{schema}.calculate_get_pixel_value', [lon, lat, 4326, 4326, schema, config["NA_aspect_table"]])
            res = curs.fetchone()[0]
            # aspect_value = session.query(func.data.calculate_get_pixel_value(lon, lat, 4326, 4326, schema, config["NA_aspect_table"])).first()

    # leaving contexts doesn't close the connection
    conn.close()

    result = res
    return result


def min_required_sea_level_for_submersion(lat, lon, table=config['NA_coastal_flooding_90m_table'], conn_params=conn_parameters):
    """
    Calculate the minimum required sea level for submersion of a specified location.

    This function connects to a PostgreSQL database using the provided connection parameters
    and calls a stored procedure 'data.calculate_coastal_flood' to determine the minimum sea level
    required for a specified location to be submerged.

    Args:
        lat (float): The latitude of the location for which the submersion is calculated.
        lon (float): The longitude of the location for which the submersion is calculated.
        conn_params (dict): A dictionary containing database connection arguments (values for host, port, dbname, user, and password).

    Returns:
        dict: Contains 'distance' (distance to the 10 m layer (in meters) if the location is out of geographical bounds of layer, up to 100 m)
        and 'level' (the minimum required sea level (in meters) for the specified location to be submerged).
    Note:
        - The function queries a PostGIS database using the 'data.calculate_coastal_flood' stored procedure.
    """
    EPSG_in = 4326
    EPSG_out = 3035 # must be metric CRS
    schema=table['schema']
    table_name=table['table_name']

    with psycopg2.connect(**conn_params) as conn:
        with conn.cursor() as cursor: 
            cursor.callproc(f'{schema}.calculate_coastal_flood', (lon, lat, EPSG_in, EPSG_out, schema, table_name))
            sea_level_distance = cursor.fetchone() # tuple - first is distance, second is sea_level

    # context manager doesn't close connection, unfortunatelly
    conn.close()

    if isinstance(sea_level_distance,tuple):
        distance, level = float(sea_level_distance[0]), sea_level_distance[1]
    else:
        distance, level = None, None
    out = {'distance_10m': distance,
        'level': level}
    return out