#################################################################
##                 Pre-process Toolik discharge data           ##
##                         Jason Clark                         ##
##                     jaclark2@alaska.edu                     ##
##                          2021-04-01                         ##
#################################################################

##---------------------------------------------------------------
##                    import discharge data                    --
##                       format for LAKE                       --
##                      save LAKE met data                     --
##---------------------------------------------------------------
library(tidyverse)

##................................................................
##                            Inputs                            ..
##................................................................

discharge_data<-'/Volumes/GoogleDrive/.shortcut-targets-by-id/1GKqJHCp9fdSJGlCQb-ZAG0N_WcJZ91E4/Toolik/met data/'
outFile<-'Toolik_inflow.dat'
outDir<-'/Volumes/GoogleDrive/.shortcut-targets-by-id/1GKqJHCp9fdSJGlCQb-ZAG0N_WcJZ91E4/Toolik/model infiles/'
dir.create(outDir)

##................................................................
##                         Begin script                         ..
##................................................................

##Discharge data
d2013<-read.csv(paste0(discharge_data,'2013_Toolik_Inlet_Kling.csv'),skip=4,col.names = c('Date_Time','Water_Temp_C','Conductivity_uScm','Q_m3sec'),header=FALSE)
d2014<-read.csv(paste0(discharge_data,'2014_Toolik_Inlet_Kling.csv'),skip=4,col.names = c('Date_Time','Water_Temp_C','Conductivity_uScm','Q_m3sec'),header=FALSE)
d2015<-read.csv(paste0(discharge_data,'2015_Toolik_Inlet_Kling.csv'),skip=4,col.names = c('Date_Time','Water_Temp_C','Conductivity_uScm','Q_m3sec'),header=FALSE)
d2016<-read.csv(paste0(discharge_data,'2016_Toolik_Inlet_Kling.csv'),skip=4,col.names = c('Date_Time','Water_Temp_C','Conductivity_uScm','Q_m3sec'),header=FALSE)

#combine
dComb<-rbind(d2013,d2014,d2015,d2016)

#format
dComb$Water_Temp_C<-as.numeric(dComb$Water_Temp_C)
dComb$DateTime<-as.POSIXct(strptime(dComb$Date_Time,format="%Y-%m-%d %H:%M",tz='US/Alaska'))
dComb$Date<-as.POSIXct(strptime(dComb$DateTime,format="%Y-%m-%d",tz='US/Alaska'))
dCombDay<- dComb %>% group_by(Date) %>% summarise(Temp=mean(Water_Temp_C,na.rm = TRUE),Q=mean(Q_m3sec,na.rm = TRUE))

#plot
ggplot(dCombDay,aes(x=Date))+geom_line(aes(y=Temp))+geom_line(aes(y=Q),color=2)
ggsave(paste0(discharge_data,'ToolikDischarge.png'),width=8,height = 8,units = 'in',dpi=450)

#create data file for LAKE
inFlow<-realD %>% mutate(timestring=as.numeric(Dates2), discharge=Q, depth=0.5, width=2, temperature=Temp, salinity=0, riverflowX=1, riverflowY=0,DOC=-999, POC=-999, DIC=-999, CH_4=0,level=-999)
inFlow<-inFlow %>% select (timestring,discharge,depth,width,temperature,salinity,riverflowX,riverflowY,DOC,POC,DIC,CH_4,level)
write.table(inFlow,paste0(outDir,outFile),sep=',',dec='.',col.names = FALSE,row.names = FALSE)

#outflow is copy of inflow discharge
outFlow<-inFlow %>% select(timestring,discharge) %>% mutate(depth0=26,depth=26.5,width=2,level=-999)
write.table(outFlow,paste0(outDir,gsub('inflow','outflow',outFile)),sep=',',dec='.',col.names = FALSE,row.names = FALSE)
