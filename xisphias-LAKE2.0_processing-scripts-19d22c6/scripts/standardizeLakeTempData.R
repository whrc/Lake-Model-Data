
#################################################################
##                 Standardize lake temp data                  ##
##                      Jason Clark                            ##
##                  jaclark2@alaska.edu                        ##
##                          2021-05-06                         ##
#################################################################

##---------------------------------------------------------------
##                        for each lake                        --
##                   read data, format, save                   --
##---------------------------------------------------------------

library(dplyr)
library(tidyr)

##................................................................
##                            Inputs                            ..
##................................................................

tooOut<-'/Volumes/GoogleDrive/.shortcut-targets-by-id/1GKqJHCp9fdSJGlCQb-ZAG0N_WcJZ91E4/Toolik/lake temp/'
foxOut<-'/Volumes/GoogleDrive/.shortcut-targets-by-id/176DvgNPhsKW7TBVpTCtZG1D-yOhdnWiv/FoxDen/lake temp/'
atqOut<-'/Volumes/GoogleDrive/.shortcut-targets-by-id/1SEdqxXdOheha0mU1ulzQQvIJks1sLfl2/Atgasuk/lake temp/'

##................................................................
##                         Begin Script                         ..
##................................................................

#toolik
{
  measured<-read.csv('/Volumes/GoogleDrive/My Drive/beavers/LAKE/Toolik_Lake_Temp_2012_2016Allinterploate.csv')
  measured$DateTime<-as.POSIXct(measured$DateTime,format="%Y-%m-%d %H:%M:%S",tz="US/Alaska")
  mStd<-measured %>% select(DateTime,Depth,Temp_C)
  write.csv(mStd,file=paste0(tooOut,'lakeTemp.csv'),row.names = FALSE)
  saveRDS(mStd, file=paste0(tooOut,'lakeTemp.RDS'))
}

#foxden
{
  measured<-read.csv('/Volumes/GoogleDrive/.shortcut-targets-by-id/176DvgNPhsKW7TBVpTCtZG1D-yOhdnWiv/FoxDen/lake temp/Fox_Den_bed_temp_2009_2013.csv')
  measured$DateTime<-as.POSIXct(measured$time,format="%m/%d/%y %H:%M",tz="US/Alaska")
  measured$Temp_C<-measured$temp
  measured$Depth<-1.5
  mStd<-measured %>% select(DateTime,Depth,Temp_C)
  write.csv(mStd,file=paste0(foxOut,'lakeTemp.csv'),row.names = FALSE)
  saveRDS(mStd, file=paste0(foxOut,'lakeTemp.RDS'))
}
#atq
{
  # https://search.dataone.org/view/doi%3A10.18739%2FA21G0HV92
  m1<-read.csv('/Volumes/GoogleDrive/.shortcut-targets-by-id/1SEdqxXdOheha0mU1ulzQQvIJks1sLfl2/Atgasuk/lake temp/atqasuk-atq-200-2012-laketimeseries-calon.csv')
  m2<-read.csv('/Volumes/GoogleDrive/.shortcut-targets-by-id/1SEdqxXdOheha0mU1ulzQQvIJks1sLfl2/Atgasuk/lake temp/atqasuk-atq-200-2013-laketimeseries.csv',
               skip=2,header=FALSE)
  names(m2)<-c('Date',     'TimeAST',  'TempC30cm','AbsPres_kPa','TempC250cm','Pres_mBar','SensorDepth_m',rep(NULL,8))
  m2<-m2[,1:7]
  m3<-read.csv('/Volumes/GoogleDrive/.shortcut-targets-by-id/1SEdqxXdOheha0mU1ulzQQvIJks1sLfl2/Atgasuk/lake temp/atqasuk-atq-200-2014-laketimeseries-winter.csv',skip=1)
  m4<-read.csv('/Volumes/GoogleDrive/.shortcut-targets-by-id/1SEdqxXdOheha0mU1ulzQQvIJks1sLfl2/Atgasuk/lake temp/atqasuk-atq-200-2015-laketimeseries.csv',
               skip=2,header=FALSE)
  names(m4)<-c('Date',     'TimeAST',  'TempC30cm','TempC325cm','Pres_mBar','SensorDepth_m',rep(NULL,10))
  m4<-m4[,1:6]
  #2012
    m1$DateTime<-as.POSIXct(paste(m1$Date,m1$Time..ADT),format='%m/%d/%Y %H:%M:%S',tz="US/Alaska")
    m1g<-m1 %>% select(DateTime,'0.3'=TempC..30cm.,'2.5'=TempC..250cm.,SensorDepth_m=Water.Level..m) %>% gather('Depth', 'Temp_C',-DateTime,-SensorDepth_m)
    m1g$Depth<-as.numeric(m1g$Depth)
  #2013
    m2$DateTime<-as.POSIXct(paste(m2$Date,m2$TimeAST),format='%m/%d/%Y %H:%M:%S',tz="US/Alaska")
    m2g<-m2 %>% select(DateTime,'0.3'=TempC30cm,'2.5'=TempC250cm,SensorDepth_m=SensorDepth_m) %>% gather('Depth', 'Temp_C',-DateTime,-SensorDepth_m)
    m2g$Depth<-as.numeric(m2g$Depth)
  #2014
    m3$DateTime<-as.POSIXct(paste(m3$Date,m3$Time..AST),format='%m/%d/%Y %H:%M:%S',tz="US/Alaska")
    m3g<-m3 %>% select(DateTime,'0.3'=TempC..30cm,'3.25'=TempC..325cm,SensorDepth_m=Sensor.Depth..m) %>% gather('Depth', 'Temp_C',-DateTime,-SensorDepth_m)
    m3g$Depth<-as.numeric(m3g$Depth)
  #2015
    m4$DateTime<-as.POSIXct(paste(m4$Date,m4$TimeAST),format='%m/%d/%Y %H:%M:%S',tz="US/Alaska")
    m4g<-m4 %>% select(DateTime,'0.3'=TempC30cm,'3.25'=TempC325cm,SensorDepth_m=SensorDepth_m) %>% gather('Depth', 'Temp_C',-DateTime,-SensorDepth_m)
    m4g$Depth<-as.numeric(m4g$Depth)
  
  #combine
  measured<-rbind(m1g,m2g,m3g,m4g)
  measured$Depth<-factor(measured$Depth)
  mStd<-measured %>% select(DateTime,Depth,Temp_C)
  write.csv(mStd,file=paste0(atqOut,'lakeTemp.csv'),row.names = FALSE)
  saveRDS(mStd, file=paste0(atqOut,'lakeTemp.RDS'))
}
# library(ggplot2)
# ggplot(mStd, aes(x=DateTime,y=Temp_C))+geom_line()+facet_wrap(~Depth,ncol=1)
