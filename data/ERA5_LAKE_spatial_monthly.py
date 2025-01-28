import datetime
from google.auth import compute_engine
import ee
import time
import os
from google.cloud import storage
from google.auth.transport.requests import AuthorizedSession
import google.auth
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import xarray as xr
from datetime import datetime
import glob
from functools import reduce
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
import numpy as np
import json

# This script downloads ERA5 data needed to run the LAKE model from Google Earth Engine, saves to google storage bucket, 
# reads the Google Earth Engine exported data, and formats it for LAKE. This script has been created and tested only in a Google Cloud VM.
# Users may encounter Earth Engine and Cloud Storage authentication issues if running elsewhere.
# 
# Requires:
#   - Earth Engine API access for both VM and user/service account.
#   - Google cloud storage bucket and full editing access for user/service account.
#   - Lon/Lat (decimal degrees WGS84) of target point location
#   - Buffer size of area to average ERA5 data surrounding point location
#   - Start and end date


# Filter the image collections
def filterCollections(ERA5_land, geometry, date_start, date_end):
    ERA5_land_filtered = (ERA5_land
                        .filterBounds(geometry)
                        .filterDate(date_start, date_end)
                        .select(['mean_2m_air_temperature', 'dewpoint_2m_temperature', 
                                'u_component_of_wind_10m', 'v_component_of_wind_10m', 
                                'surface_pressure', 'total_precipitation']))
    
    return ERA5_land_filtered

def fetch_ERA5(geometry, date_start, date_end):

    # Define the image collections
    ERA5_land = ee.ImageCollection("ECMWF/ERA5/MONTHLY")

    # Filter the image collections
    ERA5_land_filtered = filterCollections(ERA5_land, geometry, date_start, date_end)

    return ERA5_land_filtered

def export_ERA5(ERA5_land_data, bucket_name, prefix, filename_noext, region_geometry):
    
    # Get list of dates from image collection
    dates = ERA5_land_data.aggregate_array('system:time_start').getInfo()
    
    tasks = []
    # Export each image in collection
    for i, date in enumerate(dates):
        # Convert timestamp to date string
        date_str = datetime.fromtimestamp(date/1000).strftime('%Y%m%d')
        
        # Get single image for this date
        image = ee.Image(ERA5_land_data.toList(ERA5_land_data.size()).get(i))
        
        # Export data to Google Cloud Storage as GeoTIFF
        task = ee.batch.Export.image.toCloudStorage(
            image=image.clip(region_geometry),
            description=f"{filename_noext}_{date_str}",
            bucket=bucket_name,
            fileNamePrefix=os.path.join(prefix, f"{filename_noext}_{date_str}"),
            scale=11132,  # Scale in meters
            crs='EPSG:4326',
            fileFormat='GeoTIFF',
            maxPixels=1e13
        )
        task.start()
        tasks.append(task)

    return tasks

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

def cancel_running_tasks():
  ee.Authenticate()
  ee.Initialize()
  tasks = ee.batch.Task.list()
  for task in tasks:
      task_id = task.status()['id']
      task_state = task.status()['state']
      if task_state == 'RUNNING' or task_state == 'READY':
          task.cancel()
          print('Task {} canceled'.format(task_id))
      else:
          print('Task {} state is {}'.format(task_id, task_state))

def ERA5_gee_pipeline(date_start, date_end, path_to_region_geojson, bucket_name):
    
    # Initialize the Earth Engine module
    ee.Authenticate()
    ee.Initialize()

    #coords in lon, lat
    #indices of array = lon, lat, buffer
    
    # Read in region geometry from geojson file
    with open(path_to_region_geojson) as f:
        region_geojson = json.load(f)
    
    # Convert geojson to Earth Engine geometry
    region_geometry = ee.Geometry(region_geojson['features'][0]['geometry'])
    
    print('fetching data')
    ERA5_land_data = fetch_ERA5(region_geometry, date_start, date_end)

    print('exporting to cloud storage')
    tasks = export_ERA5(ERA5_land_data, bucket_name, f'{site_name}/GEE/', f'{site_name}_ERA5_Land', region_geometry)

    # Wait for both tasks to complete
    wait_for_tasks(tasks)

    return os.path.join(f'{site_name}/GEE/', f'{site_name}_ERA5_Land')

def calc_vp(temp_K):
    
    a = 2.5008e6/461.2
    b = (1/273.16) - (1/temp_K)
    
    return 611.12*np.exp(a*b) #vapor pressure in Pa
    
def calc_vp_v2(temp_K):
    
    temp_C = temp_K - 273.15
    
    return 6.11*np.power(10,(7.5*temp_C)/(237.3+temp_C)) * 100 #vapor pressure in Pa
    
def calc_mr(pres, vapor_pres):
    #pressures in Pa
    return (vapor_pres * 0.622) / (pres - vapor_pres) #kg water / kg dry air

def calc_rh(temp, dewpoint_temp):
    
    return (calc_vp(dewpoint_temp) / calc_vp(temp)) * 100 #rh, in %

def format_df_for_LAKE(ERA5_land):
    
    ERA5_land['year'] = ERA5_land['date'].dt.year
    ERA5_land['month'] = ERA5_land['date'].dt.month
    ERA5_land['day'] = ERA5_land['date'].dt.day
    ERA5_land.columns

    data_df_LAKE= ERA5_land[['u_component_of_wind_10m', 'v_component_of_wind_10m',
                              'mean_2m_air_temperature', 'dewpoint_2m_temperature',  
                              'surface_pressure', 'total_precipitation']]
    data_df_LAKE.columns=['Year','Month','Day','Uspeed','Vspeed','Temp','Hum','Pres','Precip']

    return data_df_LAKE

def process_vars_for_LAKE(path_to_ERA5_land):

    ERA5_land = pd.read_csv(path_to_ERA5_land, parse_dates=['date'])
    ERA5_land['precip_rate'] = ERA5_land['total_precipitation_sum'] / 86400 #m to m s-1

    vp =  calc_vp(ERA5_land['dewpoint_temperature_2m'])
    ERA5_land['mixing_ratio'] = calc_mr(ERA5_land['surface_pressure'], vp)
    ERA5_land['relative_humidity'] = calc_rh(ERA5_land['temperature_2m'], ERA5_land['dewpoint_temperature_2m'])

    return ERA5_land

def plot_ERA5_timeseries(ERA5_land, save_fig=None):

    var_names = ['u_component_of_wind_10m', 'v_component_of_wind_10m',
                        'temperature_2m', 'mixing_ratio', 'surface_pressure', 'precip_rate']

    fig, axes = plt.subplots(len(var_names),1, figsize=(8, 2*len(var_names)), sharex=True)

    for i, var in enumerate(var_names):
        sns.lineplot(x=ERA5_land['date'], y=ERA5_land[var], ax=axes[i])
        
    fig.autofmt_xdate()
    fig.tight_layout()

    plt.show()

    if not save_fig:
        return
    else:
        plt.savefig(save_fig, dpi=300)

def gen_ERA5_spinup(ERA5_land_for_LAKE, save_fig):

    data_df_avg = ERA5_land_for_LAKE.loc[~((ERA5_land_for_LAKE['Month']==2) & (ERA5_land_for_LAKE['Day']==29))].groupby(by=['Month', 'Day']).mean().reset_index()
    data_df_avg['Precip'] = ERA5_land_for_LAKE['Precip'].mean()
    spin_years = 20

    data_df_avg = pd.concat([data_df_avg]*spin_years, ignore_index=True)
    data_df_avg['Year'] = [2023-spin_years+i for i in range(0,spin_years) for n in range(0,365)]

    var_names = ['Uspeed','Vspeed','Temp','Hum','Pres','Precip']
    fig, axes = plt.subplots(len(var_names),1, figsize=(8, 2*len(var_names)), sharex=True)

    for i, var in enumerate(var_names):
        sns.lineplot(x=data_df_avg.index, y=data_df_avg[var], ax=axes[i])
    fig.autofmt_xdate()
    fig.tight_layout()

    plt.show()

    if not save_fig:
        return
    else:
        plt.savefig(save_fig, dpi=300)

    data_df_avg = data_df_avg.reindex(columns=['Year','Month','Day','Uspeed','Vspeed','Temp','Hum','Pres','Precip'])

    return data_df_avg


def postprocess_ERA5(gs_full_path, out_transient_path, out_transient_figure, out_spinup_path, out_spinup_figure):

    ERA5_land = process_vars_for_LAKE(gs_full_path)
    plot_ERA5_timeseries(ERA5_land, save_fig=out_transient_figure)
    ERA5_land_for_LAKE = format_df_for_LAKE(ERA5_land)
    ERA5_land_spinup = gen_ERA5_spinup(ERA5_land_for_LAKE, save_fig=out_spinup_figure)

    ERA5_land_for_LAKE.to_csv(out_transient_path, index=False)
    ERA5_land_spinup.to_csv(out_spinup_path, index=False)

    
bucket_name='lake-model-data'
site_name='YKD_full_region_monthly'
date_start = '2000-01-01'
#date_end = '2000-01-05'
date_end = '2024-09-15'
region_geometry = '/home/amullen/Upscaling/region/ykd_region_clip.geojson'


#download ERA5 from Earth Engine
gs_path = ERA5_gee_pipeline(date_start = date_start, date_end = date_end, 
                            path_to_region_geojson=region_geometry, bucket_name=bucket_name)
#cancel_running_tasks()

#print('Processing locally for LAKE')
#gs_full_path = os.path.join('gs://', bucket_name, gs_path + '.csv')
#out_transient_path = os.path.join('gs://', bucket_name, f'{site_name}/LAKE/{site_name}_transient.dat')
#out_transient_figure = f'NSF_sites/{site_name}/{site_name}_transient.jpg'

#out_spinup_path = os.path.join('gs://', bucket_name, f'{site_name}/LAKE/{site_name}_spinup.dat')
#out_spinup_figure= f'NSF_sites/{site_name}/{site_name}_spinup.jpg'

#process ERA5 locally for LAKE
#postprocess_ERA5(gs_full_path, 
#                 out_transient_path,
#                 out_transient_figure, 
#                 out_spinup_path, 
#                 out_spinup_figure)