
#################################################################
##                 Pre-process Toolik met data                 ##
##                         Jason Clark                         ##
##                     jaclark2@alaska.edu                     ##
##                          2021-04-01                         ##
#################################################################

##---------------------------------------------------------------
##                import toolik data met data,                 --
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
met_data<-'/Volumes/GoogleDrive/.shortcut-targets-by-id/1GKqJHCp9fdSJGlCQb-ZAG0N_WcJZ91E4/Toolik/met data/'
#plots are output to met data folder

#output dir for model data
outDir<-'/Volumes/GoogleDrive/.shortcut-targets-by-id/1GKqJHCp9fdSJGlCQb-ZAG0N_WcJZ91E4/Toolik/model infiles/'
dir.create(outDir)
#output met
outFile<-'TooHr'

#lake temp data dates are 8/2012-8/2016
#after prcp gauge fixed, per metadata '2013-05-16' 
minDate<-'2013-05-16' 
# maxDate<-'2015-08-30'
maxDate<-'2016-09-19'

##................................................................
##                         Begin script                         ..
##................................................................

#1hour
toolik<-read.csv(paste0(met_data,'Toolik_1_hour_data.csv'),colClasses = c('character',rep('numeric',24)),na.strings = NaN)
toolik<-toolik %>% mutate(hour2=as.character(hour),hourl=nchar(hour2),hour3=ifelse(hourl==3,paste0('0',hour2),hour2))
toolik<-toolik %>% mutate(datetime2=paste(date,hour3,sep=' '))
toolik<-toolik %>% mutate(datetime=as.POSIXct(strptime(datetime2,format = '%m/%d/%Y %H%M'),tz='US/Alaska'),
                          Date=as.POSIXct(strptime(datetime,format = '%Y-%m-%d'),tz='US/Alaska'))
toolik %>% dplyr::select(date,hour, datetime2,datetime,air_temp_5m) %>% filter(datetime>minDate) %>% as_tibble()
str(toolik)
summary(toolik)

#3hour for snow data
toolik3hr<-read.csv(paste0(met_data,'Toolik_3_hour_data.csv'),colClasses = c('character',rep('numeric',20)),na.strings = NaN)
toolik3hr<-toolik3hr %>% mutate(hour2=as.character(hour),hourl=nchar(hour2),hour3=ifelse(hourl==3,paste0('0',hour2),hour2))
toolik3hr<-toolik3hr %>% mutate(datetime2=paste(date,hour3,sep=' '))
toolik3hr<-toolik3hr %>% mutate(datetime=as.POSIXct(strptime(datetime2,format = '%m/%d/%Y %H%M'),tz='US/Alaska'),
                          Date=as.POSIXct(strptime(datetime,format = '%Y-%m-%d'),tz='US/Alaska'))
toolik3hr %>% dplyr::select(date,hour, datetime2,datetime,lake_temp_c) %>% filter(datetime>minDate) %>% as_tibble()
str(toolik3hr)
summary(toolik3hr)

#spurious rain
toolik$rain[toolik$rain>50]<-0
#add 85mbar to P for sea level correction, per metadata
toolik$barometer_mbar<-toolik$barometer_mbar+85

#wind dir is missing from 2013/12/31 to 2015/6/2
  #fill with runif
  nWind<-length(toolik[is.na(toolik$wind_dir_5m),1])
  toolik %>% mutate(wind_dir_5m2=wind_dir_5m) %>% mutate(wind_dir_5m2=ifelse(is.na(wind_dir_5m2),1,0)) %>% pull(wind_dir_5m2) %>% sum()
  toolik<-toolik %>% mutate(wind_dir_5m2=wind_dir_5m) %>% mutate(wind_dir_5m2=ifelse(is.na(wind_dir_5m2),runif(nWind,0,359),wind_dir_5m2))
  #mean wind dir is similar
  summary(toolik$wind_dir_5m)
  summary(toolik$wind_dir_5m2)

#rename
tooHrly<-toolik %>% rename(
  TAVG=air_temp_3m,
  AWND=wind_dir_5m2,
  AWNS=wind_sp_5m,
  RH=rh_3m,
  SWIncRad=sw_up_avg,
  LWIncRad=lw_up_avg,
  ASP=barometer_mbar,
  PRCP=rain,
  PRCP2=rain2) %>%
  dplyr::select(datetime, TAVG,AWND, AWNS, RH, SWIncRad, LWIncRad,ASP,PRCP,PRCP2)

#snow from 3hr, also independent measure of lake depth and 2m temperature, 2m depth is not accurate per metadata.
tooSnowHr<-toolik3hr %>% rename(
  SnowDepth=snow_depth,
  lakeDepth=lake_depth,
  lakeTemp=lake_temp_c) %>%
  dplyr::select(datetime, SnowDepth,lakeDepth, lakeTemp)
summary(tooSnowHr)

#filter to dates of record +14days for rolling avg
tooHrly<-tooHrly %>% filter(datetime>as.POSIXct(minDate)-60*60*24*14,datetime<as.POSIXct(maxDate)+60*60*24*14)
tooSnowHr<-tooSnowHr %>% filter(datetime>as.POSIXct(minDate)-60*60*24*14,datetime<as.POSIXct(maxDate)+60*60*24*14)
summary(tooHrly)

#rollave for missing data filling
dRoll<-tooHrly
dRoll<-dRoll %>% mutate(TAVG=roll_mean(TAVG,n=7*24,na.rm = TRUE,fill=0),
                        PRCP=roll_mean(PRCP,n=7*24,na.rm = TRUE,fill=0),
                        PRCP2=roll_mean(PRCP2,n=7*24,na.rm = TRUE,fill=0),
                        AWND=roll_median(AWND,n=14*24,na.rm = TRUE,fill=0),
                        AWNS=roll_mean(AWNS,n=7*24,na.rm = TRUE,fill=0),
                        RH  =roll_mean(  RH,n=7*24,na.rm = TRUE,fill=0),
                        SWIncRad  =roll_mean(  SWIncRad,n=7*24,na.rm = TRUE,fill=0),
                        LWIncRad  =roll_mean(  LWIncRad,n=7*24,na.rm = TRUE,fill=0),
                        ASP =roll_mean( ASP,n=7*24,na.rm = TRUE,fill=0))

tooHrly$TAVG[is.na(tooHrly$TAVG)]<-dRoll$TAVG[is.na(tooHrly$TAVG)]
tooHrly$PRCP[is.na(tooHrly$PRCP)]<-dRoll$PRCP[is.na(tooHrly$PRCP)]
tooHrly$PRCP2[is.na(tooHrly$PRCP2)]<-dRoll$PRCP2[is.na(tooHrly$PRCP2)]
tooHrly$AWND[is.na(tooHrly$AWND)]<-dRoll$AWND[is.na(tooHrly$AWND)]
tooHrly$AWNS[is.na(tooHrly$AWNS)]<-dRoll$AWNS[is.na(tooHrly$AWNS)]
tooHrly$RH[is.na(tooHrly$RH)]<-dRoll$RH[is.na(tooHrly$RH)]
tooHrly$SWIncRad[is.na(tooHrly$SWIncRad)]<-dRoll$SWIncRad[is.na(tooHrly$SWIncRad)]
tooHrly$LWIncRad[is.na(tooHrly$LWIncRad)]<-dRoll$LWIncRad[is.na(tooHrly$LWIncRad)]
tooHrly$ASP[is.na(tooHrly$ASP)]<-dRoll$ASP[is.na(tooHrly$ASP)]

#rollave for missing data, 2nd pass
dRoll<-tooHrly
dRoll<-dRoll %>% mutate(TAVG=roll_mean(TAVG,n=7*24*4,na.rm = TRUE,fill=mean(TAVG,na.rm=TRUE)),
                        PRCP=roll_mean(PRCP,n=7*24*4,na.rm = TRUE,fill=mean(PRCP,na.rm=TRUE)),
                        PRCP2=roll_mean(PRCP2,n=7*24*4,na.rm = TRUE,fill=mean(PRCP2,na.rm=TRUE)),
                        AWND=roll_median(AWND,n=14*24*4,na.rm = TRUE,fill=mean(AWND,na.rm=TRUE)),
                        AWNS=roll_mean(AWNS,n=7*24*4,na.rm = TRUE,fill=mean(AWNS,na.rm=TRUE)),
                        RH  =roll_mean(  RH,n=7*24*4,na.rm = TRUE,fill=mean(RH,na.rm=TRUE)),
                        SWIncRad  =roll_mean(  SWIncRad,n=7*24*4,na.rm = TRUE,fill=mean(SWIncRad,na.rm=TRUE)),
                        LWIncRad  =roll_mean(  LWIncRad,n=7*24*4,na.rm = TRUE,fill=mean(LWIncRad,na.rm=TRUE)),
                        ASP =roll_mean( ASP,n=7*24*4,na.rm = TRUE,fill=mean(ASP,na.rm=TRUE)))

tooHrly$TAVG[is.na(tooHrly$TAVG)]<-dRoll$TAVG[is.na(tooHrly$TAVG)]
tooHrly$PRCP[is.na(tooHrly$PRCP)]<-dRoll$PRCP[is.na(tooHrly$PRCP)]
tooHrly$PRCP2[is.na(tooHrly$PRCP2)]<-dRoll$PRCP2[is.na(tooHrly$PRCP2)]
tooHrly$AWND[is.na(tooHrly$AWND)]<-dRoll$AWND[is.na(tooHrly$AWND)]
tooHrly$AWNS[is.na(tooHrly$AWNS)]<-dRoll$AWNS[is.na(tooHrly$AWNS)]
tooHrly$RH[is.na(tooHrly$RH)]<-dRoll$RH[is.na(tooHrly$RH)]
tooHrly$SWIncRad[is.na(tooHrly$SWIncRad)]<-dRoll$SWIncRad[is.na(tooHrly$SWIncRad)]
tooHrly$LWIncRad[is.na(tooHrly$LWIncRad)]<-dRoll$LWIncRad[is.na(tooHrly$LWIncRad)]
tooHrly$ASP[is.na(tooHrly$ASP)]<-dRoll$ASP[is.na(tooHrly$ASP)]

#filter to dates of record
tooHrly<-tooHrly %>% filter(datetime>minDate,datetime<maxDate)
summary(tooHrly)

#missing ASP for ~2months in 2015/6 filled with 2014/5
ASP2015<-tooHrly %>% filter(datetime>='2014-09-23',datetime<='2015-01-01') %>% pull(ASP)
tooHrly$ASP[tooHrly$datetime>='2015-09-23'&tooHrly$datetime<='2016-01-01']<-ASP2015

#spurious low ASP filled with mean
tooHrly$ASP[tooHrly$ASP<750]<-mean(tooHrly$ASP,na.rm = TRUE)

#missing PRCP filled with 0
tooHrly$PRCP[is.na(tooHrly$PRCP)]<-tooHrly$PRCP2[is.na(tooHrly$PRCP)]
tooHrly$PRCP[is.na(tooHrly$PRCP)]<-0

#bad RH filled with mean
tooHrly$RH[tooHrly$RH<10]<-mean(tooHrly$RH,na.rm = TRUE)

#bad LW, filled with mean
tooHrly$LWIncRad[tooHrly$LWIncRad<100]<-mean(tooHrly$LWIncRad,na.rm = TRUE)

# #plot filled data
# tG<-gather(tooHrly,'var','value', -datetime)
# ggplot(tG,aes(x=datetime,y=value))+geom_line()+facet_wrap(~var,scales = 'free')

#convert units for model, format

  #convert to 2 component wind
  RperD<-pi / 180
  tooHrly$AWND[tooHrly$AWND==360]<-0
  tooHrly<-tooHrly %>% mutate(Ugeo = -AWNS * sin(AWND * RperD),Vgeo = -AWNS * cos(AWND * RperD))
  tooHrly<-tooHrly %>% mutate(Ugeo = ifelse(abs(Ugeo)<0.1,0,Ugeo),Vgeo = ifelse(abs(Vgeo)<0.1,0,Vgeo))

  #AT
  #Ta, K
  LakeWX<-tooHrly %>% mutate(T_K=round(TAVG+273.15,3)) %>% dplyr::select(datetime,T_K)
  
  #ASP
  #mbar to pa
  tooHrly <- tooHrly %>% mutate(ASP=ASP*100)
  LakeWX$Pres<-round(tooHrly$ASP,0)
  
  #radiation
  #rad from station
  LakeWX$SWinstrument<-tooHrly$SWIncRad
  LakeWX$LWinstrument<-tooHrly$LWIncRad
  
  #wind
  LakeWX$U<-round(tooHrly$Ugeo,3)
  LakeWX$V<-round(tooHrly$Vgeo,3)
  
  #humidity
  #RH to kg/kg
  tooHrly<-tooHrly %>% mutate(pws = ifelse(TAVG>0,610.94*exp(1)^( (17.625*(TAVG)) / ((TAVG)+243.04) ),
                                   611.21*exp(1)^( (22.587*(TAVG))/((TAVG)+273.86) ) ),
                      pw = (pws*RH)/100,
                      hr = (0.62198 * (pw)) / (ASP - (pw )))
  LakeWX$H<-round(tooHrly$hr,5)
  
  #prcp 
  #mm/day to m/s
  # LakeWX$P<-round(tooHrly$PRCP/(60*60*24*1000),10)
  #mm/hr to m/s
  LakeWX$P<-round(tooHrly$PRCP/(60*60*1000),10)

#rad is missing until 2013-5-17
#clip for instrument rad
minDateInst<-'2013-05-17'
LakeWXins<-LakeWX %>% filter(datetime>minDateInst) %>% arrange(datetime)
summary(LakeWXins)

#save for model
write.table(LakeWXins %>% select (T_K,Pres,LWinstrument,SWinstrument,U,V,H,P),paste0(outDir,outFile,'.dat'),sep=',',dec='.',col.names = FALSE,row.names = FALSE)
write.table(LakeWXins %>% select (datetime,T_K,Pres,LWinstrument,SWinstrument,U,V,H,P),paste0(met_data,outFile,'Date.dat'),sep=',',dec='.',col.names = FALSE,row.names = FALSE)
LakeWXins %>% summarise(min=min(datetime),max=max(datetime),nrows=nrow(.))
#   min                   max      nrows
# 1 2013-05-17 2016-09-18 22:00:00 29300

#plot model data
Lakeg<-gather(LakeWXins,'var','value',-datetime)
ggplot(Lakeg,aes(x=datetime,y=value))+geom_line()+facet_wrap(~var,scales = 'free')
ggsave(paste0(met_data,outFile,'.pdf'))

#snow output
outFilesn<-paste('Snow_LakeDepth')
write.table(tooSnowHr,paste0(met_data,outFilesn,'.dat'),sep=',',dec='.',col.names = FALSE,row.names = FALSE)
Lakeg<-gather(tooSnowHr,'var','value',-datetime)
ggplot(Lakeg,aes(x=datetime,y=value))+geom_line()+facet_wrap(~var,scales = 'free')
ggsave(paste0(met_data,outFilesn,'.pdf'))


# toolik data columns
# air_temp_1m         Air Temp 1m                   1-hour data    C         1
# air_temp_3m         Air Temp 3m                   1-hour data    C         1
# air_temp_5m         Air Temp 5m                   1-hour data    C         1
# barometer_mbar      Barometric Pressure           1-hour data    mbar      1
# lw_dn_avg           Longwave looking down (1hr)   1-hour data    W/m2      3
# lw_up_avg           Longwave looking up (1hr)     1-hour data    W/m2      3
# rain2               Rainfall - tipping bucket     1-hour data    mm        3
# rh_1m               Relative Humidity 1m          1-hour data    %         0
# rh_3m               Relative Humidity 3m          1-hour data    %         0
# sw_dn_avg           Shortwave looking down (1hr)  1-hour data    W/m2      3
# sw_up_avg           Shortwave looking up (1hr)    1-hour data    W/m2      3
# snow_depth_1hr      Snow Depth (1hr)              1-hour data    cm        1
# mean_wind_direction_5m_sonicSonic Mean Wind Direction     1-hour data    deg       1
# mean_wind_speed_5m_sonicSonic Mean Wind Speed         1-hour data    m/s       3
# wind_dir_5m         Wind Direction 5m             1-hour data    deg       0
# wind_sp_1m          Wind Speed 1m                 1-hour data    m/s       3
# wind_sp_5m          Wind Speed 5m                 1-hour data    m/s       3
# rain                Year-round Precipitation      1-hour data    mm        3

#format, units for LAKE
# N_Temp          1   air temperature, K |
# N_Pres          2   atmospheric pressure, Pa |
# N_LWdown        3   total solar radiation, W/m^2 |
# N_SWdown        4   atmospheric radiation, W/m^2) |
# N_Uspeed        5   x-component speed, m/s |
# N_Vspeed        6   y-component speed, m/s |  
# N_Hum           7   air humidity, kg/kg |
# N_Precip        8   precipitation, m/s |
