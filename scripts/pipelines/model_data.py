from sqlalchemy import create_engine
from sqlalchemy import Column, Integer, String, MetaData, Table
from geoalchemy2 import Geometry
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base
import sys

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

world_continents_boundaries = Table("world_continents_boundaries", 
                        metadata_obj,
                        Column("geometry", Geometry('POLYGON')),
                        autoload_with = engine)


engine = create_engine(database_url)
Base.metadata.create_all(bind=engine)

Session = sessionmaker(bind=engine)

def get_countries_boundaries(country = "United States"):
    """
    function to return boundary of a country
    """
    # session is for interacting with db
    session = Session()
    united_states_boundaries = session.query(world_countries_boundaries).filter_by(country=country).first()
    for boundary in united_states_boundaries:
        print(f"ID: {boundary.ogc_fid}, Shapename: {boundary.country}, Geometry: {boundary.geom}")
        geometry = boundary.geom
        return_data = geometry

    # Close the session
    session.close()
    return return_data