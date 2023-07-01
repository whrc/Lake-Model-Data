import pandas as pd
import datetime
import matplotlib.pyplot as plt
import numpy as np
import preproccess_USGS_data_utils as plu

def preprocess_Atqasuk(wprecip,p_reduce):
    # preprocess_Atgasuk(wprecip,p_reduce)
    # wprecip: depth of winter precipitation 
    # p_reduce (true, flase): reduce year1 winter precip
    
    # load South Meade USGS met data
    df_SM=plu.read_met_data('data/South Meade.txt')
    df_SM['snowD'].loc[df_SM['snowD']<0]=0
    #plt.plot(df_SM['Time'], df_SM['snowD'],alpha=0.5)

    # zero out suspecious snow depth
    early_snow = (df_SM['Time'] >= '2013-9-01 14:00:00') & \
            (df_SM['Time'] <= '2013-9-25 09:00:00')
    late_snow = (df_SM['Time'] >= '2014-5-11 14:00:00') & \
            (df_SM['Time'] <= '2014-8-10 09:00:00')
    yr1_snow = (df_SM['Time'] >= '2013-9-25 09:00:00') & \
            (df_SM['Time'] <= '2014-5-11 14:00:00')
    df_SM['snowD'].loc[late_snow]=0
    df_SM['snowD'].loc[early_snow]=0
    df_SM['snowD'].loc[yr1_snow] = df_SM['snowD'].loc[yr1_snow].interpolate()

    early_snow = (df_SM['Time'] >= '2014-9-01 14:00:00') & \
            (df_SM['Time'] <= '2014-10-3 09:00:00')
    late_snow = (df_SM['Time'] >= '2015-5-21 06:00:00') & \
            (df_SM['Time'] <= '2015-6-01 09:00:00')
    yr2_snow = (df_SM['Time'] >= '2014-10-3 09:00:00') & \
            (df_SM['Time'] <= '2015-5-21 06:00:00') 

    df_SM['snowD'].loc[late_snow]=0
    df_SM['snowD'].loc[early_snow]=0
    df_SM['snowD'].loc[yr2_snow] = df_SM['snowD'].loc[yr2_snow].interpolate()


    twoyears = (df_SM['Time'] >= '2013-8-12 14:00:00') & \
                (df_SM['Time'] <= '2015-8-10 09:00:00')
    dsnow_SM=df_SM['snowD'].loc[twoyears]
    dsnow_SM[dsnow_SM<0]=0
    print ('snow depth size:',dsnow_SM.size)
    sn=dsnow_SM.values-dsnow_SM.shift(periods=1).values

    #remove wind blown snow, assume min snow precip
    # wprecip: 0, 1, 2
    sn[sn<wprecip]=0 # assuming that it be could wind blown snow 
    swe=np.zeros(len(df_SM['rain'].loc[twoyears]))
    swe=np.abs(sn)*.3 # converting to SWE assuming density 0.3
    swe[np.isnan(swe)]=0
        
    # rain is in mm/h -> 1e-3/3600 m/s ~ 2.77e-07#
    #plt.plot(df_SM['Time'].loc[twoyears], swe*2.77e-07,alpha=1)
    #plt.plot(df_SM['Time'].loc[twoyears], df_SM['rain'].loc[twoyears]*2.77e-07,alpha=1)

    #plt.ylabel('precip. snow $[m SWE s^{-1}]$', fontsize=12)
    #plt.xlabel('Time [yr]', fontsize=12)
    #plt.xticks(rotation=45);

    # load Atgasuk met data and convert the date format for better match with SM USGS data
    df_met = pd.read_csv('data/atqasuk-atq-2014-2015-meteorology-timeseries-calon.csv', encoding = "ISO-8859-1",usecols=[i for i in range(0,10)])

    df1 = df_met.copy()
    df1['Year'] = pd.DatetimeIndex(df_met['Date']).year
    df1['Month'] = pd.DatetimeIndex(df_met['Date']).month
    df1['Day'] = pd.DatetimeIndex(df_met['Date']).day
    df1['Hour'] = pd.DatetimeIndex(df_met['Time, ADT']).hour
    [m,n]=df_met.shape
    print ('Atqasuk_lake_NS dimensions:',m,n)

    df_time = pd.DataFrame({'year': df1['Year'].values,
                        'month': df1['Month'].values,
                        'day': df1['Day'].values,
                        'hour': df1['Hour'].values})

    df_time=pd.to_datetime(df_time)
    frames = [df_time, df_met.iloc[:, 2:]]
    df_Atg = pd.concat(frames, axis=1)
    df_Atg.rename(columns = {list(df_Atg)[0]: 'Time'}, inplace = True)

    #combined precip data from two datasets
    twoyears = (df_SM['Time'] >= '2013-8-12 14:00:00') & \
                (df_SM['Time'] <= '2015-8-10 09:00:00')

    twoyears1 = (df_Atg['Time'] >= '2013-8-12 14:00:00') & \
                (df_Atg['Time'] <= '2015-8-10 09:00:00')
    
    swe=swe+df_SM['rain'].loc[twoyears].values
    swe_full=np.maximum(swe,df_Atg['Rain, mm (LBL: Prec)'].loc[twoyears1].values)
    swe_full[np.isnan(swe_full)]=0
    #plt.plot(df_SM['Time'].loc[twoyears], swe_full*2.77e-07,alpha=1)

    # reduce amount winter precip for year 1
    if p_reduce:
        for i in range(4000,7000):
            if swe_full[i]<1:
                swe_full[i]=0
        swe_full[800:1200]=0

    #combined wind speed data from two datasets by their max values
    df_SM['U']=df_SM['U'].loc[twoyears].fillna(0)
    df_Atg['Wind Speed, m/s (LBL: Wind)']=df_Atg['Wind Speed, m/s (LBL: Wind)'].loc[twoyears1].fillna(0)

    #plt.plot(df_SM['Time'].loc[twoyears], df_SM['U'].loc[twoyears],alpha=1)
    #plt.plot(df_Atg['Time'].loc[twoyears1], df_Atg['Wind Speed, m/s (LBL: Wind)'].loc[twoyears1],alpha=.6)
    #plt.xticks(rotation=45);
    wind_speed_full=np.maximum(df_SM['U'].loc[twoyears].values,df_Atg['Wind Speed, m/s (LBL: Wind)'].loc[twoyears1].values)
    #plt.plot(df_Atg['Time'].loc[twoyears1],wind_speed_full,alpha=0.5)

    #combined wind direction data from two datasets by their max values
    df_SM['Uang']=df_SM['Uang'].loc[twoyears].fillna(0)
    df_Atg['Wind Direction, Ã¸ (LBL: Dir)']=df_Atg['Wind Direction, Ã¸ (LBL: Dir)'].loc[twoyears1].fillna(0)
    #plt.plot(df_SM['Time'].loc[twoyears], df_SM['Uang'].loc[twoyears],alpha=1)
    #plt.plot(df_Atg['Time'].loc[twoyears1], df_Atg['Wind Direction, Ã¸ (LBL: Dir)'].loc[twoyears1],alpha=.6)
    #plt.xticks(rotation=45);
    wind_dir_full=np.maximum(df_SM['Uang'].loc[twoyears].values,df_Atg['Wind Direction, Ã¸ (LBL: Dir)'].loc[twoyears1].values)
    #plt.plot(df_Atg['Time'].loc[twoyears1],wind_dir_full, alpha=0.5)

    # update snow data 
    df_met['Wind Direction, Ã¸ (LBL: Dir)']=wind_dir_full
    df_met['Wind Speed, m/s (LBL: Wind)']=wind_speed_full
    swe_full=np.maximum(swe_full,df_Atg['Rain, mm (LBL: Prec)'].loc[twoyears1].values)


    df_Atgasuk_data=plu.USGS2LAKE_met_data(df_met,"data/Atqasuk.dat",swe_full)
