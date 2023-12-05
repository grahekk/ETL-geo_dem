import logging
import time
import datetime
from tqdm import tqdm


from .pipeline_download_s3_global import download_tree_cover_density
from .pipeline_load_localPG import import_to_local_db
from .pipeline_transform_vrt_gdal import gdal_build_vrt, absolute_file_paths, geofilter_paths_list, create_neighbour_vrt
from .pipeline_transform_qgis_resample import transform_geomorphon_qgis, geomorphon_chunky
import settings


config = settings.get_config()
schema = settings.get_schema()

logging.basicConfig(filename=config["log_file_dir"], level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Decorator function
def log_execution_time_and_args(func):
    """
    A decorator that logs the execution time, arguments, and result of a function. 
    Logs are saved into the file provided in the `logging.basicConfig`.

    Args:
        func (callable): The function to be decorated.

    Returns:
        callable: The decorated function.
    """
    def wrapper(*args, **kwargs):
        """
        Wrapper function for the decorator.

        Args:
            *args: Positional arguments passed to the decorated function.
            **kwargs: Keyword arguments passed to the decorated function.

        Returns:
            Any: The result of the decorated function.
        """
        start_time = time.time()
        logging.info(f"Started execution of {func.__name__}")
        try:
            result = func(*args, **kwargs)
            total_time = time.time() - start_time
            logging.info(
                f"{func.__name__} executed with arguments {args}, {kwargs} "
                f"in {total_time:.2f} seconds or (dd:hh:mm.ss): {datetime.timedelta(seconds=total_time)}. Result: {result}"
            )
            return result

        except Exception as e:
            logging.info(f"{func.__name__} failed with arguments {args}, {kwargs}: {str(e)}")
            raise
    return wrapper

def loop_progress(iterable, desc="Processing", unit="files"):
    total_iterations = len(iterable)
    progress_bar = tqdm(iterable, total=total_iterations, desc=desc, unit=unit)

    def decorator(func):
        def wrapper(*args, **kwargs):
            start_time = time.time()
            for item in progress_bar:
                func(item, *args, **kwargs)
                current_iteration = progress_bar.n
                elapsed_time = time.time() - start_time
                time_per_iteration = elapsed_time / current_iteration if current_iteration > 0 else 0
                time_left = (total_iterations - current_iteration) * time_per_iteration

                progress_bar.set_postfix(
                    TimeLeft=f"{time_left:.2f}s",
                    Elapsed=f"{elapsed_time:.2f}s",
                    Iteration=f"{current_iteration}/{total_iterations}"
                )
        return wrapper
    return decorator


class GisDataDownloader:
    """
    Class for downloading GIS data.

    Attributes:
        download_url (str): The URL from which to download GIS data.

    Methods:
        download(): Downloads GIS data from the specified URL. Decorator `@log_execution_time_and_args` is used for logging
    """
    def __init__(self, download_urls):
        self.download_urls = download_urls

    @log_execution_time_and_args
    def download(self):
        """
        Downloads GIS data from the specified URL.

        Returns:
            Any: The downloaded GIS data.
        """
        logging.info(f"Downloading data from {self.download_urls}")
        # Logic for downloading GIS data
        # ...

    @log_execution_time_and_args
    def download_tree_cover_density(self):
        """
        Download raster tiles for tree cover density data.
        """
        download_tree_cover_density(extent="North America")


class DataTransformer:
    """
    Class for transforming downloaded gis data.

    Attributes:
        data (Any): File path pointing to the data to be transformed.

    Methods:
        transform(): Transforms the data (generic).
        build_vrt(): Builds virtual raster
    """
    def __init__(self, data:str):
        self.data = data

    @log_execution_time_and_args
    def transform(self):
        """
        Transforms the data using some of methods.

        Returns:
            Any: The transformed data.
        """
        # Logic for transforming data
        logging.info("Transforming data")
        # ...

    @log_execution_time_and_args
    def build_vrt(self, output_vrt:str, extent = "standard"):
        """
        builds virtual raster from set of files

        Args:
            output_vrt (str): path to the output virtual raster (in config file)
            extent (str): choose whether geographical extent should or shouldn't be reshaped
            
        Returns:
            Any: The transformed data, vrt built and saved into `output_vrt` path.
        """
        input_dir = self.data
        if extent == "standard":
            input_files = absolute_file_paths(input_dir, ".tif") # returns generator object
            files_list = sorted(list(input_files))
        else:
            files_list = geofilter_paths_list(input_dir, coc = "country")
        gdal_build_vrt(files_list, output_vrt)
        return print("Task done!")
    
    @log_execution_time_and_args
    def geomorphon(self, output_file:str):
        """
        Transform given raster dataset using geomorphon algorithm in `transform_geomorphon_qgis`.

        Args:
            output_file(str): file path to write output raster to
        """
        dem_file = self.data
        transform_geomorphon_qgis(dem_file, output_file)

    @log_execution_time_and_args
    def geomorphon_class(self):
        """
        Transform given raster dataset using geomorphon algorithm in `transform_geomorphon_qgis`.
        Chunky style.

        Args:
            output_file(str): file path to write output raster to
        """
        dem_files = self.data
        for file in tqdm(dem_files):
            vrt_path = create_neighbour_vrt(file, config["NA_geomorphon_dir"])
            geomorphon_chunky(file, vrt_path)
            logging.info(f"vrt file {vrt_path} created as well as geomorphon for it")


class DataLoader:
    """
    Class for loading data into a database.

    Attributes:
        table_name (str): 
        data (str): File path pointing to the data to be transformed.
        database_name (str): The name of the database.
        data_type (str): The type of the data to be loaded (raster(tif) or vector (gpkg, shp)).

    Methods:
        load(): Loads data into the specified database.
    """
    def __init__(self, table_name:str, data:str, database_name:str, data_type:str):
        self.schema = schema
        self.table_name = table_name
        self.data = data
        self.database_name = database_name
        self.data_type = data_type

    @log_execution_time_and_args
    def load_local(self):
        """
        Loads data into the specified local database.

        Returns:
            Any: The result of the loading operation.
        """
        return import_to_local_db(self.schema, self.table_name, self.data, self.database_name, self.data_type)


class MainPipeline:
    """
    Main pipeline class that orchestrates the download, transform, and load processes.

    Attributes:
        download_pipeline (GisDataDownloader): The download pipeline instance.
        transform_pipeline (DataTransformer): The transform pipeline instance.
        load_pipeline (DataLoader): The load pipeline instance.

    Methods:
        execute(): Executes the entire pipeline, including download, transform, and load.
    """
    def __init__(self, download_url, database_url):
        self.download_pipeline = GisDataDownloader(download_url)
        self.transform_pipeline = DataTransformer(None)  # Pass actual data during execution
        self.load_pipeline = DataLoader(None, database_url)  # Pass actual data during execution

    @log_execution_time_and_args
    def execute(self):
        """
        Executes the entire pipeline, including download, transform, and load.

        Returns:
            Any: The result of the entire pipeline execution.
        """
        downloaded_data = self.download_pipeline.download()
        transformed_data = self.transform_pipeline.transform(downloaded_data)
        self.load_pipeline.load(transformed_data)
