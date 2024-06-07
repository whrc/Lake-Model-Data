import datetime
from google.auth import compute_engine
import ee
import time

# Initialize the Earth Engine module
ee.Authenticate()
ee.Initialize()

# Filter the image collections
def filterCollections(ERA5, ERA5_land, buffer, date_start, date_end):

    ERA5_filtered = (ERA5
                    .filterBounds(buffer)
                    .filterDate(date_start, date_end)
                    .select(['mean_2m_air_temperature', 'dewpoint_2m_temperature', 
                            'total_precipitation', 'surface_pressure', 
                            'u_component_of_wind_10m', 'v_component_of_wind_10m']))

    ERA5_land_filtered = (ERA5_land
                        .filterBounds(buffer)
                        .filterDate(date_start, date_end)
                        .select(['temperature_2m', 'dewpoint_temperature_2m', 
                                'u_component_of_wind_10m', 'v_component_of_wind_10m', 
                                'surface_pressure', 'total_precipitation_sum', 
                                'surface_solar_radiation_downwards_sum', 
                                'surface_thermal_radiation_downwards_sum']))
    
    return ERA5_filtered, ERA5_land_filtered

# Function to reduce region and create feature for ERA5 data
def create_ERA5_feature(img, buffer):
    stats = img.reduceRegion(reducer=ee.Reducer.mean(), geometry=buffer, scale=27830, crs='EPSG:4326')
    return ee.Feature(None, {
        'date': ee.Date(img.get('system:time_start')).format('YYYY-MM-dd'),
        'mean_2m_air_temperature': stats.get('mean_2m_air_temperature', -9999),
        'dewpoint_2m_temperature': stats.get('dewpoint_2m_temperature', -9999),
        'total_precipitation': stats.get('total_precipitation', -9999),
        'surface_pressure': stats.get('surface_pressure', -9999),
        'u_component_of_wind_10m': stats.get('u_component_of_wind_10m', -9999),
        'v_component_of_wind_10m': stats.get('v_component_of_wind_10m', -9999)
    })

# Function to reduce region and create feature for ERA5_land data
def create_ERA5_land_feature(img, buffer):
    stats = img.reduceRegion(reducer=ee.Reducer.mean(), geometry=buffer, scale=11132, crs='EPSG:4326')
    return ee.Feature(None, 
        {
        'date': ee.Date(img.get('system:time_start')).format('YYYY-MM-dd'),
        'temperature_2m': stats.get('temperature_2m', -9999),
        'dewpoint_temperature_2m': stats.get('dewpoint_temperature_2m', -9999),
        'u_component_of_wind_10m': stats.get('u_component_of_wind_10m', -9999),
        'v_component_of_wind_10m': stats.get('v_component_of_wind_10m', -9999),
        'surface_pressure': stats.get('surface_pressure', -9999),
        'total_precipitation_sum': stats.get('total_precipitation_sum', -9999),
        'surface_solar_radiation_downwards_sum': stats.get('surface_solar_radiation_downwards_sum', -9999),
        'surface_thermal_radiation_downwards_sum': stats.get('surface_thermal_radiation_downwards_sum', -9999)
        })



# Define the image collections
ERA5 = ee.ImageCollection("ECMWF/ERA5/DAILY")
ERA5_land = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")

#coords in lon, lat
#indices of array = lat, lon, buffer

Brazil_sites = {'DosMeinacos':[-53.125, -12.555, 2000] #verified in ee
                }
	
NSF_sites = {'TKL873':[-162.5878, 67.6070, 1800],
             'TKL524':[-162.5756, 67.5003, 1800], #verified in ee
             'TKL917':[], 
             'TKL884':[], 
             'TKL807':[], 
             'SitukuyukBP':[], 
             'RabitCreekBP':[-163.6089, 67.4964, 1800], #verified in ee
             'EliRiverBP':[-161.8608, 67.8064, 1800], #verified in ee
             'DosMeinacos':[-53.125, -12.555, 2000] #verified in ee
             }

site = 'EliRiverBP'
# Define the date range
date_start = '2000-01-01'
date_end = '2024-01-01'

point = ee.Geometry.Point(NSF_sites[site][:2])
buffer = point.buffer(NSF_sites[site][2])

# Filter the image collections
ERA5_filtered, ERA5_land_filtered = filterCollections(ERA5, ERA5_land, buffer, date_start, date_end)

ERA5_data = ERA5_filtered.map(lambda img: create_ERA5_feature(img, buffer))
ERA5_land_data = ERA5_land_filtered.map(lambda img: create_ERA5_land_feature(img, buffer))

# Export data to Google Cloud Storage
task1 = ee.batch.Export.table.toCloudStorage(
    collection=ERA5_data,
    description=f'{site}_ERA5',
    bucket='lake-model-data',
    fileNamePrefix=f'{site}/GEE/{site}_ERA5',
    fileFormat='CSV'
)
task1.start()

task2 = ee.batch.Export.table.toCloudStorage(
    collection=ERA5_land_data,
    description=f'{site}_ERA5_Land',
    bucket='lake-model-data',
    fileNamePrefix=f'{site}/GEE/{site}_ERA5_Land',
    fileFormat='CSV'
)
task2.start()

def wait_for_tasks(task_list):
    """Wait for tasks to complete and print their status."""
    for task in task_list:
        while True:
            status = task.status()
            state = status['state']
            if state in ['COMPLETED', 'FAILED']:
                print(f"Task {task.id} finished with state: {state}")
                if state == 'FAILED':
                    print(f"Error message: {status['error_message']}")
                break
            else:
                print(f"Task {task.id} is in state: {state}")
                time.sleep(10)  # Wait for 10 seconds before checking again

# Wait for both tasks to complete
wait_for_tasks([task1, task2])
