from sqlalchemy import create_engine
from sqlalchemy import Column, Integer, String, MetaData, Table
from geoalchemy2 import Geometry
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base
import sys
from sqlalchemy import func

sys.path.append("/home/nikola/4_north_america/GeoDataPump/scripts/")
import settings
database_url = settings.get_database_url()
schema = settings.get_schema()

engine = create_engine(database_url)
metadata_obj = MetaData(schema=schema)

Base = declarative_base()

# class CountryBoundary(Base):
#     __tablename__ = f'{schema}.world_countries_boundaries'
#     id = Column(Integer, primary_key=True)
#     shapename = Column(String)
#     geometry = Column(Geometry(geometry_type='GEOMETRY', srid=4326))    


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


engine = create_engine(database_url)
Base.metadata.create_all(bind=engine)

Session = sessionmaker(bind=engine)

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