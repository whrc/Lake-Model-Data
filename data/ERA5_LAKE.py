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

SITES = {'test': [-113.9586, 46.9339, 1800],
             'TKL873':[-162.5878, 67.6070, 1800],
             'TKL524':[-162.5756, 67.5003, 1800], #verified in ee
             'TKL917':[-162.6062, 67.61894, 1800], 
             'TKL884':[-162.5845, 67.61177, 1800], 
             'TKL807':[], 
             'YKD_two_ponds': [-163.230, 61.260, 1800],
             'SitukuyukBP':[-163.17819, 67.154665, 1800], 
             'RabitCreekBP':[-163.6089, 67.4964, 1800], #verified in ee
             'EliRiverBP':[-161.8608, 67.8064, 1800], #verified in ee
             'DosMeinacos':[-53.125, -12.555, 2000] #verified in ee
            }


# Filter the image collections
def filterCollections(ERA5_land, buffer, date_start, date_end):
    ERA5_land_filtered = (ERA5_land
                        .filterBounds(buffer)
                        .filterDate(date_start, date_end)
                        .select(['temperature_2m', 'dewpoint_temperature_2m', 
                                'u_component_of_wind_10m', 'v_component_of_wind_10m', 
                                'surface_pressure', 'total_precipitation_sum', 
                                'surface_solar_radiation_downwards_sum', 
                                'surface_thermal_radiation_downwards_sum']))
    
    return ERA5_land_filtered

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

def fetch_ERA5_summary(geometry, date_start, date_end):

    # Define the image collections
    ERA5_land = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")

    # Filter the image collections
    ERA5_land_filtered = filterCollections(ERA5_land, geometry, date_start, date_end)

    # average over geometry
    ERA5_land_data = ERA5_land_filtered.map(lambda img: create_ERA5_land_feature(img, geometry))

    return ERA5_land_data

def export_ERA5(ERA5_land_data, bucket_name, prefix, filename_noext):
    
    # Export data to Google Cloud Storage
    task = ee.batch.Export.table.toCloudStorage(
        collection=ERA5_land_data,
        description=filename_noext,
        bucket=bucket_name,
        fileNamePrefix=os.path.join(prefix, filename_noext),
        fileFormat='CSV'
    )
    task.start()

    return task

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

def ERA5_gee_pipeline(date_start, date_end, site_name, bucket_name):
    
    # Initialize the Earth Engine module
    ee.Authenticate()
    ee.Initialize()

    #coords in lon, lat
    #indices of array = lon, lat, buffer
    
    if site_name not in SITES.keys():
        print('site not defined in SITES')
        return

    point = ee.Geometry.Point(SITES[site_name][:2])
    buffer = point.buffer(SITES[site_name][2])

    print('fetching data')
    ERA5_land_data = fetch_ERA5_summary(buffer, date_start, date_end)

    print('exporting to cloud storage')
    task = export_ERA5(ERA5_land_data, bucket_name, f'{site_name}/GEE/', f'{site_name}_ERA5_Land')

    # Wait for both tasks to complete
    wait_for_tasks([task])

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

    data_df_LAKE= ERA5_land[['year', 'month', 'day', 'u_component_of_wind_10m', 'v_component_of_wind_10m',
                        'temperature_2m', 'mixing_ratio', 'surface_pressure', 'shortwave', 'longwave', 'precip_rate']]
    data_df_LAKE.columns=['Year','Month','Day','Uspeed','Vspeed','Temp','Hum','Pres','SWdown','LWdown','Precip']

    return data_df_LAKE

def process_vars_for_LAKE(path_to_ERA5_land):

    ERA5_land = pd.read_csv(path_to_ERA5_land, parse_dates=['date'])

    ERA5_land['shortwave'] = ERA5_land['surface_solar_radiation_downwards_sum'] / 86400 #J m-2 (daily cumulative rad) -> daily average rad in w m-2
    ERA5_land['longwave'] = ERA5_land['surface_thermal_radiation_downwards_sum'] / 86400 #J m-2 (daily cumulative rad) -> daily average rad in w m-2
    ERA5_land['precip_rate'] = ERA5_land['total_precipitation_sum'] / 86400 #m to m s-1

    vp =  calc_vp(ERA5_land['dewpoint_temperature_2m'])
    ERA5_land['mixing_ratio'] = calc_mr(ERA5_land['surface_pressure'], vp)
    ERA5_land['relative_humidity'] = calc_rh(ERA5_land['temperature_2m'], ERA5_land['dewpoint_temperature_2m'])

    return ERA5_land

def plot_ERA5_timeseries(ERA5_land, save_fig=None):

    var_names = ['u_component_of_wind_10m', 'v_component_of_wind_10m',
                        'temperature_2m', 'mixing_ratio', 'surface_pressure', 'shortwave', 'longwave', 'precip_rate']

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

    var_names = ['Uspeed','Vspeed','Temp','Hum','Pres','SWdown','LWdown','Precip']
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

    data_df_avg = data_df_avg.reindex(columns=['Year','Month','Day','Uspeed','Vspeed','Temp','Hum','Pres','SWdown','LWdown','Precip'])

    return data_df_avg


def postprocess_ERA5(gs_full_path, out_transient_path, out_transient_figure, out_spinup_path, out_spinup_figure):

    ERA5_land = process_vars_for_LAKE(gs_full_path)
    plot_ERA5_timeseries(ERA5_land, save_fig=out_transient_figure)
    ERA5_land_for_LAKE = format_df_for_LAKE(ERA5_land)
    ERA5_land_spinup = gen_ERA5_spinup(ERA5_land_for_LAKE, save_fig=out_spinup_figure)

    ERA5_land_for_LAKE.to_csv(out_transient_path, index=False)
    ERA5_land_spinup.to_csv(out_spinup_path, index=False)

    
bucket_name='lake-model-data'
site_name='YKD_two_ponds'
#date_start = '2000-01-01'
#date_end = '2024-01-01'
date_start = '2019-01-01'
date_end = '2024-09-01'

# TKL873 x
# TKL524 x
# TKL917 x
# TKL884 x
# TKL807 x
# SitukuyukBP
# RabitCreekBP
# EliRiverBP
# DosMeinacos

if not os.path.exists(f'NSF_sites/{site_name}/'):
    os.mkdir(f'NSF_sites/{site_name}/')

#download ERA5 from Earth Engine
gs_path = ERA5_gee_pipeline(date_start = date_start, date_end = date_end, site_name=site_name, bucket_name=bucket_name)

print('Processing locally for LAKE')
gs_full_path = os.path.join('gs://', bucket_name, gs_path + '.csv')
out_transient_path = os.path.join('gs://', bucket_name, f'{site_name}/LAKE/{site_name}_transient.dat')
out_transient_figure = f'NSF_sites/{site_name}/{site_name}_transient.jpg'

out_spinup_path = os.path.join('gs://', bucket_name, f'{site_name}/LAKE/{site_name}_spinup.dat')
out_spinup_figure= f'NSF_sites/{site_name}/{site_name}_spinup.jpg'

#process ERA5 locally for LAKE
postprocess_ERA5(gs_full_path, 
                 out_transient_path,
                 out_transient_figure, 
                 out_spinup_path, 
                 out_spinup_figure)