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


def get_geocellid(country = "United States", coc = "country"):
    """
    function to return geocellid of a specified country

    Parameters:
        country (str): country name to be queried from database
        coc (str): (coc - Continent or country) either "continent" or "country" - a switch whether to filter by continent or a country
    """
    # session is for interacting with db
    session = Session()

    if coc == "country":
        boundary = world_countries_boundaries

    else:
        boundary = world_continents_boundaries
        
    geocells = session.query(esa_grid).join(boundary, esa_grid.c.geometry.intersects(boundary.c.geom)).filter(boundary.c.country == country).all()
    return_data = [i.geocellid for i in geocells]

    # Close the session
    session.close()
    return return_data
