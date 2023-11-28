import subprocess
import sys
import psycopg2


# load gpkg to local db with ogr2ogr

def import_to_local_db(schema:str, table_name:str, source_path:str, dbname:str, format_type:str):
    """
    Import data from a source file (raster, geopackage, or shapefile) into a local PostgreSQL database.

    Args:
        schema (str): The name of the schema in the local database.
        table_name (str): The name of the table to store the imported data.
        source_path (str): The file path to the source data file (raster, geopackage, or shapefile).
        dbname (str): The name of the local PostgreSQL database.
        format_type (str): The format type ('raster', 'geopackage', or 'shapefile').

    This function uses the appropriate utility (e.g., 'raster2pgsql', 'ogr2ogr', or 'shp2pgsql') to import data
    from a source file into a local PostgreSQL database. The imported data is stored in the specified schema and table.

    Note:
    - Ensure that the necessary utility for the specified format is installed and available in the system.
    - For 'format_type', use 'raster' for raster data, 'geopackage' for geopackage data, and 'shapefile' for shapefile data.
    """
    if format_type == 'raster':
        # import_command = f'raster2pgsql -f "PostgreSQL" PG:"dbname={dbname}" -nln {schema}.{table_name} -nlt PROMOTE_TO_MULTI {source_path}'
        tile_size = "80x80"
        import_command = f'raster2pgsql -s {4326} -C -M -I -Y {source_path} -t {tile_size} | psql -d {dbname}'
    elif format_type == 'geopackage':
        import_command = f'ogr2ogr -f "PostgreSQL" PG:"dbname={dbname}" -nln {schema}.{table_name} -nlt PROMOTE_TO_MULTI {source_path}'
    elif format_type == 'shapefile':
        import_command = f'shp2pgsql -I -s {4326} -g geom {source_path} {schema}.{table_name} | psql -d {dbname}'
    else:
        raise ValueError("Invalid format_type. Use 'raster', 'geopackage', or 'shapefile'.")

    p = subprocess.Popen(import_command, stdout=subprocess.PIPE, shell=True).communicate()
    print(p)


# dump locally
def dump_local_db(schema:str, table_name:str, dump_database_name:str, dump_file:str):
    """
    Create a database dump file for a specific table in a local PostgreSQL database.

    Args:
        schema (str): The name of the schema in the local database.
        table_name (str): The name of the table to be dumped.
        dump_database_name (str): The name of the local PostgreSQL database to be dumped.
        dump_file (str): The file path to save the database dump.
    
    """
    pg_dump_command = f"pg_dump -d {dump_database_name} -t {schema}.{table_name} -Fc -f {dump_file}"
    p2 = subprocess.Popen(pg_dump_command, stdout=subprocess.PIPE, shell=True).communicate()
    print(f"pg dump went: {p2}")


# load to aws, it doesn't work with -W (password needs to be written like input in shell)
def load_to_aws(dump_file:str, host_name:str, restore_database_name:str, username:str):
    """
    Restore a database dump file to a remote PostgreSQL database on AWS.

    Args:
        dump_file (str): The file path to the database dump file.
        host_name (str): The host name or IP address of the remote database on AWS.
        restore_database_name (str): The name of the remote PostgreSQL database to be restored.
        username (str): The username for authentication.

    Note:
        pg_restore doesn't work with -W option from python, password needs to be like input in shell
    """
    pg_restore_command = f'pg_restore -j 4 -h {host_name} -d {restore_database_name} -U {username} {dump_file} -I --no-acl --no-owner'
    p3 = subprocess.Popen(pg_restore_command, stdout=subprocess.PIPE, shell=True).communicate()
    print(f"pg restore code: {p3}")


# if __name__== "__main__":
#     gpkg_to_local_db(schema, table_name)
#     dump_local_db(schema, table_name, dump_database_name, dump_file)
#     load_to_aws(dump_file, host_name, restore_database_name, username)
