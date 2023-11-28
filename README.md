# North America - Data Pipeline and ETL Project

This project is designed to help you organize and automate the process of downloading, transforming, and loading spatial data for North America from various sources into a PostgreSQL database. It follows a modular structure for better maintainability and scalability. 
It contains scripts for downloading spatial files for North America


## Directory Structure

- **data**: This directory is intended for storing both raw and transformed data. Raw data should be placed in the `raw_data` directory, and transformed data should be placed in the `transformed_data` directory. Also, it contains helper files contaning lists of urls etc.

    - **url_txt_files**: `USGS_dem_urls.txt` is a file created from USGS web page after selecting all available products, as described on the page: https://apps.nationalmap.gov/uget-instructions/. There is also `NASA_global_dem_urls.txt` file is downloaded from nasa's page as described on : https://nsidc.org/data/user-resources/help-center/nasa-earthdata-cloud-data-access-guide
as well as downoad.sh file with urls pointing to dem tiles

    - **dem files**: sources for this pipelines are US, NASA and Canadian official pages. More is described in *Data Sources* headline.

- **scripts**: This directory contains Python scripts that make up your data pipeline and ETL process.

    - **pipeline_download**: This subdirectory contains modules related to data download.

        - `pipeline_download_Canada_DEM.py`: Functions for downloading Canadian DEM data.
        - `pipeline_download_USGS_DEM.py`: Functions for downloading USGS DEM data.
        - `pipeline_download_utils.py`: Utility functions shared among download modules.
        - `pipeline_download_WMS.py`: Functions for downloading data from a Web Map Service (WMS).

    - `main.py`: The main script that orchestrates the entire data pipeline. It imports functions from the pipeline modules and executes them in the desired order.

## How to Use

1. Set up a virtual environment and install the required dependencies.

   ```bash
   # Create a virtual environment
   python -m venv venv

   # Activate the virtual environment
   source venv/bin/activate  # On Windows, use `venv\Scripts\activate`

   # Install dependencies
   pip install -r requirements.txt
    ```

2. Customize the download modules (pipeline_download_*.py) to suit your specific data sources and requirements.

3. Define transformation and loading modules as needed and import them into main.py for orchestration. Check and adjust the `config.yaml` file.

4. Run the data pipeline using the following command:

    ```
    python scripts/main.py
    ```
5. Monitor and maintain the pipeline by adding logging and error handling as necessary.

## Data Sources
The pipelines in this project serve for collecting data generally for North America from the following sources:

- USGS: For USA DEM (Digital Elevation Model) and forest resources data.
- CEC (Commission for Environmental Cooperation): For land cover data of North America.
- Canadian Government's Website: For DEM (Digital Elevation Model) data of Canada, including HRDEM and for soil data of Canada
- Canadian Government's Website: For soil data of Canada
- US Web Soil Survey: For soil data of the USA.
- NASA's Earth Data: For global DEM (Digital Elevation Model) data.
- ESA world cover: land cover data for whole world, available on: https://registry.opendata.aws/esa-worldcover-vito/
- ESA global dem: global dem in 90m and 30m resolution, downloaded from https://sentinels.copernicus.eu/web/sentinel/-/copernicus-dem-new-direct-data-download-access
- Eurostat: countries boundaries from https://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/administrative-units-statistical-units/countries 
- Geoboundaries: countries boundaries from https://www.geoboundaries.org/index.html#features-wrapper 
- Geoboundaries: continent boundaries from https://hub.arcgis.com/datasets/esri::world-continents/about


## Target PostgreSQL Database
The collected data is intended to be transformed and loaded into a given PostgreSQL database. You should configure the database connection details and schema as needed in your transformation and loading scripts.

## Final words

This README file provides an overview of your project's structure and instructions on how to use it. You can further customize it to include specific details about your project's goals, data sources, and any additional setup or configuration steps.





