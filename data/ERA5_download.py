import cdsapi
import requests
import os

# CDS API script to use CDS service to retrieve daily ERA5* variables and iterate over
# all months in the specified years.
 
# Requires:
# 1) the CDS API to be installed and working on your system
# 2) You have agreed to the ERA5 Licence (via the CDS web page)
# 3) Selection of required variable, daily statistic, etc
 
# Output:
# 1) separate netCDF file for chosen daily statistic/variable for each month
out_dir='/home/amullen/Lake-Model-Data/data/raw/TKL873'

c = cdsapi.Client(timeout=300, sleep_max = 20)
 
# Uncomment years as required
 
#years =  ['1982', '1983', '1984', '1985', '1986', '1987', '1988', '1989', '1990', '1991', '1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999', '2000']
years =  ['2021','2022', '2023']
months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
 
# For valid keywords, see Table 2 of:
# https://datastore.copernicus-climate.eu/documents/app-c3s-daily-era5-statistics/C3S_Application-Documentation_ERA5-daily-statistics-v2.pdf
 
# select your variable; name must be a valid ERA5 CDS API name.
varlist = ['2m_temperature', 
           '2m_dewpoint_temperature', 
           'surface_pressure',
           '10m_u_component_of_wind',
           '10m_v_component_of_wind',
           'surface_solar_radiation_downwards', 
           'surface_thermal_radiation_downwards',     
           'mean_total_precipitation_rate']

# Select the required statistic, valid names given in link above
stat = "daily_mean"
 
# Loop over years and months
for var in varlist:
    for yr in years:
        for mn in months:
            
            file_name = os.path.join(out_dir, "download_" + stat + "_" + var + "_" + yr + "_" + mn + ".nc")
            
            if os.path.isfile(file_name):
              print('file exists, skipping')
              continue
              
            result = c.service(
            "tool.toolbox.orchestrator.workflow",
            params={
                 "realm": "user-apps",
                 "project": "app-c3s-daily-era5-statistics",
                 "version": "master",
                 "kwargs": {
                     "dataset": "reanalysis-era5-single-levels",
                     "product_type": "reanalysis",
                     "variable": var,
                     "statistic": stat,
                     "year": yr,
                     "month": mn,
                     "time_zone": "UTC+00:0",
                     "frequency": "1-hourly",
    #
    # Users can change the output grid resolution and selected area
    #
                     "grid": "0.1/0.1",
                     "area":{"lat": [67.55, 67.65], "lon": [-162.63, -162.53]}
     
                     },
            "workflow_name": "application"
            })
            
            location=result[0]['location']
            res = requests.get(location, stream = True)
            print("Writing data to " + file_name)
            with open(file_name,'wb') as fh:
                for r in res.iter_content(chunk_size = 1024):
                    fh.write(r)
            fh.close()
