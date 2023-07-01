import h5py
import pandas as pd
import datetime
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal

def read_met_data(dfile):
    df_snow=pd.read_csv(dfile,sep='\s+')
    [m,n]=df_snow.shape
    print ('matrix dimensions:',m,n)

    df_time = pd.DataFrame({'year': df_snow.Year.values,
                            'month': df_snow.Mon.values,
                            'day': df_snow.Day.values,
                            'hour': df_snow.Hour.values})
    df_time=pd.to_datetime(df_time)
    frames = [df_time, df_snow.iloc[:, 4:]]
    df_data = pd.concat(frames, axis=1)
    df_data.columns = ["Time", "Tair", "U", "Uang", "So_d",   
                   "So_u", "rain", "snowD", "Tg5cm", "Tg10cm", "Tg15cm", "Tg20cm",
                   "Tg25cm", "Tg30cm", "Tg45cm", "Tg70cm", "Tg95cm", "Tg120cm", "Smoist"]
    return df_data

def plot_met_data(ddata):
    fig, axs = plt.subplots(4, 2, figsize=(10, 6))
    axs[0, 0].plot(ddata['Time'],ddata['Tair'])
    axs[0, 0].set_title('Tair')
    axs[0, 1].plot(ddata['Time'],ddata['So_d'], 'tab:orange')
    axs[0, 1].set_title('Downwelling Shortwave Flux')
    axs[1, 0].plot(ddata['Time'],ddata['So_u'], 'tab:green')
    axs[1, 0].set_title('Upwelling Shortwave Flux')
    axs[1, 1].plot(ddata['Time'],ddata['rain'], 'tab:red')
    axs[1, 1].set_title('rain [mm/h]')
    axs[2, 0].plot(ddata['Time'],ddata['snowD'], 'tab:red')
    axs[2, 0].set_title('snowD [cm]')
    axs[2, 1].plot(ddata['Time'],ddata['Smoist'], 'tab:red')
    axs[2, 1].set_title('Smoist [-]')
    axs[3, 0].plot(ddata['Time'],ddata['U'], 'tab:red')
    axs[3, 0].set_title('Wind speed [m/s]')
    axs[3, 1].plot(ddata['Time'],ddata['Uang'], 'tab:red')
    axs[3, 1].set_title('Wind Direction [degrees]')
    plt.tight_layout()
    

def longwave(Tskin):
    #NOTE: Tskin in Kelvin 
    eps=0.97 #bulk emissivity of land surface
    sigma=5.6686e-6 #Stefan‐Boltzmann constant 5.6686 × 10−6 (W m−2 K−4);
    #Tskin=df['Temp,°C (LBL: AirTemp)']+273.15
    Lw=eps*sigma*Tskin**4/100
    return Lw
    
def compare_data_column(data1,data2,Tcolumn,mask1,mask2):

    plt.plot(data1['Time'].loc[mask1], data1[Tcolumn].loc[mask1])
    plt.plot(data2['Time'].loc[mask2], data2[Tcolumn].loc[mask2],alpha=0.7)
    
    plt.title('Data series comparison')
    plt.legend(['DewPoint','Ikpikpuk']);
    plt.xlabel('Time [yr]')
    plt.ylabel(Tcolumn)

def plot_compare_data(data1,data2,mask1,mask2):
    
    #mask1 = (data2['Time'] >= '2006-9-1 14:00:00') & \
    #        (data2['Time'] <= '2017-9-1 09:00:00')
    #mask = (data1['Time'] >= '2006-9-1 14:00:00') & \
    #        (data1['Time'] <= '2017-9-1 09:00:00')
    fig, axs = plt.subplots(3, 2, figsize=(10, 6))
    axs[0, 0].plot(data1['Time'].loc[mask1], data1['Tair'].loc[mask1])
    axs[0, 0].plot(data2['Time'].loc[mask2], data2['Tair'].loc[mask2],alpha=0.7)
    axs[0, 0].set_title('Tair')
    axs[0, 1].plot(data1['Time'].loc[mask1], data1['So_d'].loc[mask1])
    axs[0, 1].plot(data2['Time'].loc[mask2], data2['So_d'].loc[mask2],alpha=0.7)   
    axs[0, 1].set_title('Downwelling Shortwave Flux')
    axs[1, 0].plot(data1['Time'].loc[mask1], data1['So_u'].loc[mask1])
    axs[1, 0].plot(data2['Time'].loc[mask2], data2['So_u'].loc[mask2],alpha=0.7)   
    axs[1, 0].set_title('Upwelling Shortwave Flux')
    axs[1, 1].plot(data1['Time'].loc[mask1], data1['rain'].loc[mask1])
    axs[1, 1].plot(data2['Time'].loc[mask2], data2['rain'].loc[mask2],alpha=0.7)  
    axs[1, 1].set_title('rain [mm/h]')
    #axs[2, 0].plot(ddata['Time'],ddata['snowD'], 'tab:red')
    axs[2, 0].plot(data1['Time'].loc[mask1], data1['snowD'].loc[mask1])
    axs[2, 0].plot(data2['Time'].loc[mask2], data2['snowD'].loc[mask2],alpha=0.7)  
    axs[2, 0].set_title('snowD [cm]')
    #axs[2, 1].plot(ddata['Time'],ddata['Smoist'], 'tab:red')
    axs[2, 1].plot(data1['Time'].loc[mask1], data1['Smoist'].loc[mask1])
    axs[2, 1].plot(data2['Time'].loc[mask2], data2['Smoist'].loc[mask2],alpha=0.7)  
    axs[2, 1].set_title('Smoist [-]')
    plt.legend(['DewPoint','Ikpikpuk']);
    plt.tight_layout()
    
def preprocess2ats(df,swe):
    # fill nans by interpolating and zeroing out
    [m,n]=df.shape
    data=np.zeros((m,8))
    df['Tair'] = df['Tair'].interpolate() # fill nans
    data[:,0]=df['Tair'].values+273.15
    data[:,1]=longwave(data[:,0])
    df['So_d'] = df['So_d'].interpolate() # fill nans
    data[:,2]=df['So_d'].values
    
    swe[np.isnan(swe)]=0
    snow = np.where(df["Tair"] < 0, swe, 0)
    rain = np.where(df["Tair"] >= 0, df['rain'].values, 0)
    rain[np.isnan(rain)]=0
    
    # rain is in mm/h -> 1e-3/3600 m/s ~ 2.77e-07
    data[:,3]=rain*2.77e-07 # precip rain
    data[:,4]=snow*2.77e-07 # precip snow
    data[:,5]=0.8*np.ones(m) # humidity
    data[:,6]=range(m) 
    df['U'] = df['U'].interpolate() # fill nans
    data[:,7]=df['U'].values
    
    df1=pd.DataFrame(data=data) 
    df1.columns = ["Temp [K]","LWdown [W/m2]","SWdown [W/m2]","rain [m/s]",
                        "snow [m/s]","Hum [kg/kg]","Time [s]","Wind speed [m/s]"]
    return df1

def hourly_to_daily(data):
    [m,n]=np.shape(data)
    h=list(data.columns.values)
    print ('hourly matrix:',m,n)
    ndays=m//24
    data_days=np.zeros((ndays,n))
    print ('days:',ndays)
    t_days=np.zeros(ndays)
    for j in range(n):
        for i in range(ndays):
            #data_days[i,j]=np.mean(data[i*24:24*(i+1),j])
            data_days[i,j]=np.mean(data[h[j]].loc[i*24:24*(i+1)])
    
    df1=pd.DataFrame(data=data_days) 
    df1.columns = ["Temp [K]","LWdown [W/m2]","SWdown [W/m2]","rain [m/s]",
                        "snow [m/s]","Hum [kg/kg]","Time [s]","Wind speed [m/s]"]
    return df1

def smooth_timeseries(df,w):
    df1=df.copy()
    df1["Temp [K]"] = scipy.signal.savgol_filter(df1["Temp [K]"], w, 2, mode='wrap')
    df1["LWdown [W/m2]"] = scipy.signal.savgol_filter(df1["LWdown [W/m2]"], w, 2, mode='wrap')
    df1["SWdown [W/m2]"] = scipy.signal.savgol_filter(df1["SWdown [W/m2]"], w, 2, mode='wrap')
    df1["Wind speed [m/s]"] = scipy.signal.savgol_filter(df1["Wind speed [m/s]"], w, 2, mode='wrap')
    m_precip_rain = df1["rain [m/s]"].mean()
    m_precip_snow = df1["snow [m/s]"].mean()
    snow = np.where(df1["Temp [K]"] < 273.15, m_precip_snow, 0)
    rain = np.where(df1["Temp [K]"] >= 273.15, m_precip_rain, 0)
    df1["rain [m/s]"] = rain
    df1["snow [m/s]"] = snow

    return df1
    
def save2ats_met_file(df,filename):
    [m,n]=df.shape
    print (m,n)
    time = np.arange(0, 86400.0*m, 86400.0)
    with h5py.File(filename,'w') as out:
        out.create_dataset("time [s]", data=time)
        out.create_dataset("air temperature [K]", data=df["Temp [K]"].values)
        out.create_dataset("relative humidity [-]", data=df["Hum [kg/kg]"].values)
        out.create_dataset("precipitation rain [m s^-1]", data=df["rain [m/s]"].values)
        out.create_dataset("precipitation snow [m SWE s^-1]", data=df["snow [m/s]"].values)
        out.create_dataset("incoming shortwave radiation [W m^-2]", data=df["SWdown [W/m2]"].values)
        out.create_dataset("incoming longwave radiation [W m^-2]", data=df["LWdown [W/m2]"].values)
        out.create_dataset("wind speed [m s^-1]", data=df["Wind speed [m/s]"].values)
        out.attrs.create("wind speed reference height [m]", data=3.0) # fix me

def plot_ats_h5_met_data(filename):
    d = h5py.File(filename,'r')
    print(d.keys())

    fax = plt.subplots(4,2,figsize=(10,20))
    ax = fax[1].ravel()

    time = d['time [s]'][:]
    print (d['time [s]'][:].shape)

    for i,k in enumerate([k for k in d.keys() if not "reference height" in k]):
        print (d[k])
        ax[i].plot(time/86400.0/365.0, d[k][:], 'b')
        ax[i].set_title(k);
        
def relative_humidity(Tair):
    # not working, need pressure data
    # relative humidity via Buck equation (eg.) https://en.wikipedia.org/wiki/Vapour_pressure_of_water
    Ta_c = Ta - 273.15
    vp_sat = 611.21 * np.exp((18.678 - Ta_c/234.5)*(Ta_c/(257.14-Ta_c)))
    rh = df['vp (Pa)']/vp_sat

    ka0_ = 16.635764
    ka_ = -6096.9385
    kb_ = -2.7111933e-2
    kc_ = 1.673952e-5
    kd_ = 2.433502
    vp_sat2 = 100.0*np.exp(ka0_ + ka_/Ta + (kb_ + kc_*Ta)*Ta + kd_*np.log(Ta))
    rh2 = df['vp (Pa)']/vp_sat2

    print (df['vp (Pa)'].min(), df['vp (Pa)'].max())
    print (vp_sat.min(), vp_sat.max())
    print (vp_sat2.min(), vp_sat2.max())

    print (rh.min(), rh.max())
    print (rh2.min(), rh2.max())

    #plt.plot(df['vp (Pa)'], 'k')
    #plt.plot(vp_sat,'b')
    #plt.plot(vp_sat2, 'r')
    #plt.show()

    # a few bad data points with 0 vapor pressure
    RH = np.where(rh2 > 1, 1, np.where(rh2 <= 0, np.nan, rh2))
    for i in range(len(RH)):
        if np.isnan(RH[i]):
            j = (j for j in range(i,len(RH)) if RH[j] > 0).next()
            RH[i] = RH[j]

    # spinup RH -- daily averages of the full dataset
    fax = plt.subplots(2,1)

    RH_raw = RH.reshape((-1,365)).transpose()
    for i in range(RH_raw.shape[1]):
        fax[1][0].plot(RH_raw[:,i])


    RH_m = RH_raw.mean(axis=1)
    fax[1][0].plot(RH_m,'k',linewidth=3)

    import scipy.signal
    RH_sm = scipy.signal.savgol_filter(RH_m, 61, 2, mode='wrap')
    fax[1][1].plot(RH_m, 'k')
    fax[1][1].plot(RH_sm, 'b')

def USGS2LAKE_met_data(df,fLAKEname,swe):
    swe[np.isnan(swe)]=0
    [m,n]=df.shape
    pi_180= np.pi/180
    df['Wind Direction, Ã¸ (LBL: Dir)']=df['Wind Direction, Ã¸ (LBL: Dir)'].interpolate()
    wind_dir=df['Wind Direction, Ã¸ (LBL: Dir)'].copy()
    wind_dir[wind_dir==360]=0
    df['Wind Speed, m/s (LBL: Wind)']=df['Wind Speed, m/s (LBL: Wind)'].interpolate()
    wind_speed=df['Wind Speed, m/s (LBL: Wind)'].copy()
    u=-wind_speed*np.cos(wind_dir*pi_180)
    v=-wind_speed*np.sin(wind_dir*pi_180)
    u[abs(u)<0.1]=0
    v[abs(v)<0.1]=0

    df['Temp,°C (LBL: AirTemp)']=df['Temp,°C (LBL: AirTemp)'].interpolate()
    TAVG=df['Temp,°C (LBL: AirTemp)'].copy() # this is valid for the Celsuis temperatures 
    RH=df['RH, % (LBL: RelHum)'].copy()
    ASP=df['Pressure, mbar (LBL: Press)'].copy()*100
    N=len(TAVG)
    hr=np.zeros(N)
    for i in range(N):
        if TAVG[i]>0:
            pws=610.94*np.exp(1)**(17.625*TAVG[i]/(TAVG[i]+243.04))
        else:
            pws=611.21*np.exp(1)**(22.587*TAVG[i]/(TAVG[i]+273.86))
        pw = (pws*RH[i])/100
        hr[i] = 0.62198*pw/(ASP[i] - pw)

    #Estimation of surface longwave radiation components from ground‐based historical net radiation and weather data
    eps=0.97 #bulk emissivity of land surface
    sigma=5.6686e-6 #Stefan‐Boltzmann constant 5.6686 × 10−6 (W m−2 K−4);
    Tskin=df['Temp,°C (LBL: AirTemp)']+273.15
    Lw=eps*sigma*Tskin**4/100
    # rain is in mm/h -> 1e-3/3600 m/s ~ 2.77e-07#
    precip=swe*2.77e-07

    pres=df['Pressure, mbar (LBL: Press)'].copy()*100
    pres[pres<97000]=101325 
    pres[pres>104500]=101325
    data=np.zeros((m,n-2))
    np.shape(data[:,0])
    data[:,0]=df['Temp,°C (LBL: AirTemp)']+273.15
    data[:,1]=pres
    data[:,2]=Lw
    data[:,3]=df['Solar Radiation, W/mÂ² (LBL: SolarUP)']
    data[:,4]=u
    data[:,5]=v
    data[:,6]=hr
    data[:,7]=precip

    df1=pd.DataFrame(data=data) 
    df1.columns = ["Temp [K]","Pres [Pa]","LWdown [W/m2]","SWdown [W/m2]",
                            "Uspeed [m/s]","Vspeed [m/s]","Hum [kg/kg]","Precip [m/s]"]
    df1.plot(subplots=True, layout=(4,2),figsize=(12, 10));
    plt.tight_layout()
    df1.to_csv(fLAKEname,index=False,header=False,float_format='%.3e')
    return df1