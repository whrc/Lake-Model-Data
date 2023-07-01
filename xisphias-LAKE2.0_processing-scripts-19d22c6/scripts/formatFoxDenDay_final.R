
#################################################################
##                 Pre-process Toolik met data                 ##
##                         Jason Clark                         ##
##                     jaclark2@alaska.edu                     ##
##                          2021-04-01                         ##
#################################################################

##---------------------------------------------------------------
##                import foxden met data                       --
##               fix data errors, missing values               --
##                    convert units for LAKE                   --
##                       format for LAKE                       --
##                      save LAKE met data                     --
##---------------------------------------------------------------

library(tidyverse)
library(RcppRoll)
library(ncdf4)

##................................................................
##                            Inputs                            ..
##................................................................

#met data folder
met_data<-'/Volumes/GoogleDrive/.shortcut-targets-by-id/176DvgNPhsKW7TBVpTCtZG1D-yOhdnWiv/FoxDen/met data/'
#plots are output to met data folder

#output dir for model data
outDir<-'/Volumes/GoogleDrive/.shortcut-targets-by-id/176DvgNPhsKW7TBVpTCtZG1D-yOhdnWiv/FoxDen/model infiles/'
dir.create(outDir)
#output met
outFile<-'FoxDenDay'

#lake temp data dates are: 
minDate<-as.POSIXct('2009-06-10',tz='UTC')
maxDate<-as.POSIXct('2013-06-11',tz='UTC')

##................................................................
##                         Begin script                         ..
##................................................................

## read noaa data
noaaWX<-read.csv(paste0(met_data,'2277028.csv'))
unique(noaaWX$STATION)
## select kotz station
##70133026616 Kotz
Station<-70133026616
kotz<-noaaWX %>% filter(STATION==Station)
kotz$DateTime<-as.POSIXct(strptime(kotz$DATE,format = '%Y-%m-%dT%H:%M:%S'),tz='UTC')
## subset to date range
kotz<-kotz %>% filter(DateTime>minDate,DateTime<maxDate)
## select relevant cols
dWXall<-kotz %>% select(DateTime,TAVG=DailyAverageDryBulbTemperature,
                        TMIN=DailyMinimumDryBulbTemperature,
                        TMAX=DailyMaximumDryBulbTemperature,
                        PRCP=DailyPrecipitation,
                        AWND=DailySustainedWindDirection,
                        AWNS=DailyAverageWindSpeed,
                        RH=DailyAverageRelativeHumidity,
                        AWBT=DailyAverageWetBulbTemperature,
                        ASP=DailyAverageStationPressure)
dWX<-dWXall %>% mutate(Date=(as.Date(DateTime,format='%Y-%m-%s')),Time=format(DateTime,format='%H:%M:%S')) %>% 
  filter(Time=='23:59:00')

#fix precip
dWX$PRCP[dWX$PRCP=='T']<-0.001
dWX$PRCP[dWX$PRCP=='']<-0
dWX$PRCP<-as.numeric(dWX$PRCP)
#fill errors with rollmean
  dRoll<-dWX
  dRoll<-dRoll %>% mutate(TAVG=roll_mean(TAVG,n=7,na.rm = TRUE,fill=0),
                          TMIN=roll_mean(TMIN,n=7,na.rm = TRUE,fill=0),
                          TMAX=roll_mean(TMAX,n=7,na.rm = TRUE,fill=0),
                          PRCP=roll_mean(PRCP,n=7,na.rm = TRUE,fill=0),
                          AWND=roll_median(AWND,n=14,na.rm = TRUE,fill=0),
                          AWNS=roll_mean(AWNS,n=7,na.rm = TRUE,fill=0),
                          RH  =roll_mean(  RH,n=7,na.rm = TRUE,fill=0),
                          AWBT=roll_mean(AWBT,n=10,na.rm = TRUE,fill=0),
                          ASP =roll_mean( ASP,n=7,na.rm = TRUE,fill=0))

  dWX$TAVG[is.na(dWX$TAVG)]<-dRoll$TAVG[is.na(dWX$TAVG)]
  dWX$TMIN[is.na(dWX$TMIN)]<-dRoll$TMIN[is.na(dWX$TMIN)]
  dWX$TMAX[is.na(dWX$TMAX)]<-dRoll$TMAX[is.na(dWX$TMAX)]
  dWX$PRCP[is.na(dWX$PRCP)]<-dRoll$PRCP[is.na(dWX$PRCP)]
  dWX$AWND[is.na(dWX$AWND)]<-dRoll$AWND[is.na(dWX$AWND)]
  dWX$AWNS[is.na(dWX$AWNS)]<-dRoll$AWNS[is.na(dWX$AWNS)]
  dWX$RH[is.na(dWX$RH)]<-dRoll$RH[is.na(dWX$RH)]
  dWX$AWBT[is.na(dWX$AWBT)]<-dRoll$AWBT[is.na(dWX$AWBT)]
  dWX$ASP[is.na(dWX$ASP)]<-dRoll$ASP[is.na(dWX$ASP)]
  #fill missing TAVG with avg of TminTmax
  dWX$TAVG[is.na(dWX$TAVG)]<-(dWX$TMIN[is.na(dWX$TAVG)]+dWX$TMAX[is.na(dWX$TAVG)])/2
#convert units
  WXm<-dWX 
  summary(WXm)
  WXm<-WXm %>% mutate(TAVG=(TAVG-32)*(5/9), #C
                      TMAX=(TMAX-32)*(5/9), #C
                      TMIN=(TMIN-32)*(5/9), #C
                      AWBT=(AWBT-32)*(5/9), #C
                      PRCP=PRCP*0.0254,     #m
                      AWNS=AWNS*0.44704,    #m/s
                      ASP=ASP*3386.39)      #pa
  summary(WXm)
#convert to 2 component wind
  RperD<-pi / 180
  WXm$AWND[WXm$AWND==360]<-0
  WXm<-WXm %>% mutate(Ugeo = -AWNS * sin(AWND * RperD),Vgeo = -AWNS * cos(AWND * RperD))
  WXm<-WXm %>% mutate(Ugeo = ifelse(abs(Ugeo)<0.1,0,Ugeo),Vgeo = ifelse(abs(Vgeo)<0.1,0,Vgeo))

dWXg<-gather(WXm,'var','value',-DateTime,-Date,-Time)

#format for LAKE
# N_Temp          1   air temperature, K |
# N_Pres          2   atmospheric pressure, Pa |
# N_LWdown        3   total solar radiation, W/m^2 |
# N_SWdown        4   atmospheric radiation, W/m^2) |
# N_Uspeed        5   x-component speed, m/s |
# N_Vspeed        6   y-component speed, m/s |  
# N_Hum           7   air humidity, kg/kg |
# N_Precip        8   precipitation, m/s |

#Ta, K
LakeWX<-WXm %>% mutate(T_K=round(TAVG+273.15,3)) %>% select(DateTime,T_K)
#pres
LakeWX$Pres<-round(WXm$ASP,0)
#rad from CERES
  #foxden  66.558800°, -164.456642°
  ncin<-nc_open(paste0(met_data,'CERES_SYN1deg-Day_Terra-Aqua-MODIS_Ed4.1_Subset_20000301-20200630_foxden.nc'))
  print(ncin) # metadata
  names(ncin$var)
  lon <- ncvar_get(ncin,"lon")
  nlon <- dim(lon)
  lat <- ncvar_get(ncin,"lat")
  nlat <- dim(lat)
  # get time
  time <- ncvar_get(ncin,"time")
  # time
  tunits <- ncatt_get(ncin,"time","units")
  nt <- dim(time)
  nt
  tunits
  ncDate<-as.POSIXct(time*24*60*60, origin='2000-03-01 00:00:00',tz='UCT')
  #to match datetime on WX
  ncDate2<-ncDate+(12*60*60)-60
  # get variable
  dname<-'adj_atmos_sw_down_all_surface_daily'
  d.array <- ncvar_get(ncin,dname)
  dlname <- ncatt_get(ncin,dname,"long_name")
  dunits <- ncatt_get(ncin,dname,"units")
  fillvalue <- ncatt_get(ncin,dname,"_FillValue")
  dname<-'adj_atmos_lw_down_all_surface_daily'
  d.array2 <- ncvar_get(ncin,dname)
  rad<-data.frame(DateTime=ncDate2,LW=round(d.array2,3),SW=round(d.array,3))
  LakeWX<-merge.data.frame(LakeWX,rad)

  #wind
LakeWX$U<-round(WXm$Ugeo,3)
LakeWX$V<-round(WXm$Vgeo,3)
#humidity, kg/kg
WXm<-WXm %>% mutate(pws = ifelse(TAVG>0,610.94*exp(1)^( (17.625*(TAVG)) / ((TAVG)+243.04) ),
                                 611.21*exp(1)^( (22.587*(TAVG))/((TAVG)+273.86) ) ),
                    pw = (pws*RH)/100,
                    hr = (0.62198 * (pw)) / (ASP - (pw )))
LakeWX$H<-round(WXm$hr,5)

#prec m/day to m/s
  LakeWX$P<-round(WXm$PRCP/(60*60*24),10)

#date
LakeWX<-LakeWX %>% mutate(Date=(as.Date(DateTime,format='%Y-%m-%s')))

#save for model
write.table(LakeWX %>% select (T_K,Pres,LW,SW,U,V,H,P),paste0(outDir,outFile,'.dat'),sep=',',dec='.',col.names = FALSE,row.names = FALSE)
write.table(LakeWX %>% select (DateTime,T_K,Pres,LW,SW,U,V,H,P),paste0(met_data,outFile,'Date.dat'),sep=',',dec='.',col.names = FALSE,row.names = FALSE)
LakeWX %>% summarise(min=min(DateTime),max=max(DateTime),nrows=nrow(.))
#                   min                 max nrows
# 1 2009-06-10 23:59:00 2013-06-10 23:59:00  1510

#plot model data
Lakeg<-gather(LakeWX,'var','value',-DateTime)
ggplot(Lakeg,aes(x=DateTime,y=value))+geom_line()+facet_wrap(~var,scales = 'free')
ggsave(paste0(met_data,outFile,'.pdf'))


# #365
# outFile<-'FoxDen_kotz365'
# outDir<-paste0('/Users/jac/Documents/GitHub/LAKE/data/',outFile,'/');dir.create(outDir)
# LakeWX<-LakeWX %>% slice(1:365)
# write.table(LakeWX %>% select (T_K,Pres,LW,SW,U,V,H,P),paste0(outDir,outFile,'.dat'),sep=',',dec='.',col.names = FALSE,row.names = FALSE)
# # write.table(LakeWX %>% select (T_K,Pres,LW,SW,U,V,H,P),paste0(outDir,'FoxDenKotz','.dat'),sep=',',dec='.',col.names = FALSE,row.names = FALSE)
# head(LakeWX);paste0(outDir,outFile,'.dat')
# 
# Lakeg<-dWXg<-gather(LakeWX,'var','value',-DateTime,-Date)
# ggplot(Lakeg,aes(x=DateTime,y=value))+geom_line()+facet_wrap(~var,scales = 'free')
# ggsave(paste0(outDir,outFile,'.pdf'))

# WT03 - Thunder
# WV01 - Fog, ice fog, or freezing fog (may include heavy fog)
# WT04 - Ice pellets, sleet, snow pellets, or small hail"
# PRCP - Precipitation
# WT05 - Hail (may include small hail)
# WV03 - Thunder
# WT06 - Glaze or rime
# WT07 - Dust, volcanic ash, blowing dust, blowing sand, or blowing obstruction
# WT08 - Smoke or haze
# SNWD - Snow depth
# WT09 - Blowing or drifting snow
# WDF2 - Direction of fastest 2-minute wind
# WDF5 - Direction of fastest 5-second wind
# WT10 - Tornado, waterspout, or funnel cloud"
# PGTM - Peak gust time
# WT11 - High or damaging winds
# TMAX - Maximum temperature
# WT13 - Mist
# WSF2 - Fastest 2-minute wind speed
# FMTM - Time of fastest mile or fastest 1-minute wind
# WSF5 - Fastest 5-second wind speed
# SNOW - Snowfall
# WT14 - Drizzle
# WT15 - Freezing drizzle
# WT16 - Rain (may include freezing rain, drizzle, and freezing drizzle)"
# TOBS - Temperature at the time of observation
# WT17 - Freezing rain
# WT18 - Snow, snow pellets, snow grains, or ice crystals
# WT19 - Unknown source of precipitation
# AWND - Average wind speed
# WT21 - Ground fog
# WT22 - Ice fog or freezing fog
# WV20 - Rain or snow shower
# WT01 - Fog, ice fog, or freezing fog (may include heavy fog)
# WESD - Water equivalent of snow on the ground
# WSFI - Highest instantaneous wind speed
# WT02 - Heavy fog or heaving freezing fog (not always distinguished from fog)
# TAVG - Average Temperature.
# TMIN - Minimum temperature
# MDPR - Multiday precipitation total (use with DAPR and DWPR, if available)
# TSUN - Total sunshine for the period