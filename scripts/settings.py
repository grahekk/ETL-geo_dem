import os
from dotenv import load_dotenv
import yaml
import sys
from sqlalchemy import create_engine

scripts_directory = os.path.dirname(os.path.abspath(__file__))
parent_directory = os.path.dirname(scripts_directory)
sys.path.append(parent_directory)
config_path = os.path.join(scripts_directory, 'config.yaml')
env_path = os.path.join(scripts_directory, '.env')

with open(config_path, 'r') as file:
    config = yaml.safe_load(file)

load_dotenv(env_path)

class ConnectionParameters():
    """
    A class to store connection parameters for a database connection.

    Attributes:
        host_name (str): The host name or IP address for the database server.
        port (str): The port number for the database server.
        database_name (str): The name of the database to connect to.
        username (str): The username for authentication.
        password (str): The password for authentication.
        schema (str): The default schema to use in the database.
        conn_parameters (dict): A dictionary containing connection parameters with keys:
            - 'host': The host name or IP address.
            - 'port': The port number.
            - 'database': The database name.
            - 'user': The username.
            - 'password': The password.

    Note:
    - The values of the attributes are typically retrieved from environment variables using 'os.getenv'.
    - You can access the connection parameters as class attributes, e.g., 'ConnectionParameters.host_name'.
    """
    host_name = os.getenv('host_name')
    port = os.getenv('port')
    database_name = os.getenv('database_name')
    username = os.getenv('username')
    password = os.getenv('password')
    schema = 'data'

    aws_host_name = os.getenv('aws_host_name')
    aws_port = os.getenv('aws_port')
    aws_database_name = os.getenv('aws_database_name')
    aws_username = os.getenv('aws_username')
    aws_password = os.getenv('aws_password')

    conn_parameters = {'host':host_name,
                        'port':port,
                        'database':database_name,
                        'user':username,
                        'password':password}
    
    aws_conn_parameters = {'host':aws_host_name,
                        'port':aws_port,
                        'database':aws_database_name,
                        'user':aws_username,
                        'password':aws_password}
        
    database_url: str = f"postgresql://{username}:{password}@{host_name}:{port}/{database_name}"
    aws_database_url: str = f"postgresql://{aws_username}:{aws_password}@{aws_host_name}:{aws_port}/{aws_database_name}"
    engine = create_engine(database_url)

class ConfigurationParameters():
    """
    A class for storing configuration parameters.

    Attributes:
        config (config object): The configuration object that holds various settings and options, mostly paths to files.
    """
    config = config
    data_directory_path = os.path.join(parent_directory, 'data/')
    parent_directory = parent_directory
    log_file_path = os.path.join(parent_directory, 'pipeline.log')
    ESA_dem_files_regex_match = r'_(S|N)(\d+)_00_(W|E)(\d+)'
    USGS_dem_files_regex_match = r'_(s|n)(\d+)(w|e)(\d+)_'
    
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

    global_bounds = {
    "min_lat": -90.0,   
    "max_lat": 90.0,   
    "min_lon": -180.0,  
    "max_lon": 180.0    
    }


def get_conn_parameters():
    """
    Retrieve the database connection parameters.

    Returns:
        dict: A dictionary containing connection parameters with keys:
            - 'host': The host name or IP address for the database server.
            - 'port': The port number for the database server.
            - 'database': The name of the database to connect to.
            - 'user': The username for authentication.
            - 'password': The password for authentication.

    This function retrieves the connection parameters for the database from the 'ConnectionParameters' class
    and returns them as a dictionary. The connection parameters are typically obtained from environment variables.
    """
    conn_parameters = ConnectionParameters.conn_parameters
    return conn_parameters

def get_aws_conn_parameters():
    """
    Similar to `get_conn_parameters`, only with conn parameters for connecting to aws db.
    """
    return ConnectionParameters.aws_conn_parameters

def get_database_url():
    """
    get database url made up of connection parameters
    """
    database_url = ConnectionParameters.database_url
    return database_url

def get_aws_database_url():
    """
    get database url made up of connection parameters
    """
    aws_database_url = ConnectionParameters.aws_database_url
    return aws_database_url


def get_schema():
    """
    Retrieve the default schema for the database.

    Returns:
        str: The name of the default schema to use in the database.

    This function retrieves and returns the name of the default schema to be used in database operations.
    The schema name is obtained from the 'ConnectionParameters' class.
    """
    return ConnectionParameters.schema


def get_config():
    """
    Retrieve the configuration object that holds yaml attributes.

    Returns:
        config object: The configuration object that holds various settings and options.

    This function retrieves and returns the configuration object, which contains various settings and options
    for application configuration. The configuration object is obtained from the 'ConfigurationParameters' class.
    """
    return ConfigurationParameters().config

def get_config_params():
    """
    initialize ConfigurationParameters object
    """
    return ConfigurationParameters()



def get_data_path():
    """
    Retrieve path to the data folder where temporary files are stashed

    Returns:
        data_directory_path: variable that holds string to data folder
    """
    return ConfigurationParameters.data_directory_path


def get_parent_dir_path():
    """
    Retrieve path to the parent folder where all files are stashed

    Returns:
        parent_directory(str): variable that holds string to parent folder
    """
    return ConfigurationParameters.parent_directory

def get_engine():
    """
    Get the engine initialised in ConnectionParameters class
    """
    return ConnectionParameters.engine