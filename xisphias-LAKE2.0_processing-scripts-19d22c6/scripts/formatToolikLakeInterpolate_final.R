#################################################################
##                 Pre-process Toolik lake temp data           ##
##                         Jason Clark                         ##
##                     jaclark2@alaska.edu                     ##
##                          2021-04-01                         ##
#################################################################
library(tidyverse)
library(akima)

##----------------------------------------------------------------
##                            Inputs                            --
##----------------------------------------------------------------
data_dir<-'/Volumes/GoogleDrive/.shortcut-targets-by-id/1GKqJHCp9fdSJGlCQb-ZAG0N_WcJZ91E4/Toolik/lake temp/'

##...............................................................
##                      read data                              ..
##...............................................................
m2012w<-read.csv(paste0(data_dir,'2012_2013_winter_toolik_temperature.csv'),colClasses = rep('character',19))
m2013s<-read.csv(paste0(data_dir,'2013_summer_toolik_temperature.csv'     ),colClasses = rep('character',24))
m2013w<-read.csv(paste0(data_dir,'2013_2014_winter_toolik_temperature.csv'),colClasses = rep('character',17))
m2014s<-read.csv(paste0(data_dir,'2014_summer_toolik_temperature.csv'     ),colClasses = rep('character',24))
m2014w<-read.csv(paste0(data_dir,'2014_2015_winter_toolik_temperature.csv'),colClasses = rep('character',19))
m2015s<-read.csv(paste0(data_dir,'2015_summer_toolik_temperature.csv'     ),colClasses = rep('character',23))
m2015w<-read.csv(paste0(data_dir,'2015_2016_winter_toolik_temperature.csv'),colClasses = rep('character',20))
m2016s<-read.csv(paste0(data_dir,'2016_summer_toolik_temperature.csv'     ),colClasses = rep('character',24))

##...............................................................
##                      set year for data                      ..
##...............................................................

#2012-3
{
  m2012w$Year[m2012w$DOY<366]<-2012
  m2012w$Year[m2012w$DOY>=366]<-2013
  m2012w$DOY2<-as.numeric(m2012w$DOY)
  m2013s$DOY2<-as.numeric(m2013s$DOY)
  m2012w$DOY2[m2012w$DOY>=366]<-m2012w$DOY2[m2012w$DOY2>=366]-365
  m2013s$Year<-2013
  summary(m2013s)
  head(m2013s)
  m2013s2<-m2013s[m2013s$DOY>max(m2012w$DOY2[m2012w$Year==2013]),]
}
#2013-4
{
  m2013w$Year[m2013w$DOY<366]<-2013
  m2013w$Year[m2013w$DOY>=366]<-2014
  m2013w$DOY2<-as.numeric(m2013w$DOY)
  m2014s$DOY2<-as.numeric(m2014s$DOY)
  m2013w$DOY2[m2013w$DOY>=366]<-m2013w$DOY2[m2013w$DOY2>=366]-365
  m2014s$Year<-2014
  m2014s2<-m2014s[m2014s$DOY>max(m2013w$DOY2[m2013w$Year==2014]),]
}
#2014-5
{
  m2014w$Year[m2014w$DOY<366]<-2014
  m2014w$Year[m2014w$DOY>=366]<-2015
  m2014w$DOY2<-as.numeric(m2014w$DOY)
  m2015s$DOY2<-as.numeric(m2015s$DOY)
  m2014w$DOY2[m2014w$DOY>=366]<-m2014w$DOY2[m2014w$DOY2>=366]-365
  m2015s$Year<-2015
  m2015s2<-m2015s[m2015s$DOY>max(m2014w$DOY2[m2014w$Year==2015]),]
}
#2015-6
{
  m2015w$Year[m2015w$DOY<366]<-2015
  m2015w$Year[m2015w$DOY>=366]<-2016
  m2015w$DOY2<-as.numeric(m2015w$DOY)
  m2016s$DOY2<-as.numeric(m2016s$DOY)
  m2015w$DOY2[m2015w$DOY>=366]<-m2015w$DOY2[m2015w$DOY2>=366]-365
  m2016s$Year<-2016
  m2016s2<-m2016s[m2016s$DOY>max(m2015w$DOY2[m2015w$Year==2016]),]
}


##...............................................................
##                       format depths,                        ..
##                      convert datatime,                      ..
##                    aggregate to hourly,                     ..
##                    convert to unixtime,                     ..
##                 interpolate to 1m intervals                 ..
##        Interpolation is time/processor intensive!!!         ..
##...............................................................
#2012winter
{
  #calc depths
  m1g<-m2012w %>% pivot_longer(c(-Year,-DOY2,-DOY),names_to='Depth',values_to='Temp_C')
  m1g$Depth<-gsub('X','',m1g$Depth)
  m1g<-m1g %>% mutate(Depth=as.numeric(Depth),Temp_C=as.numeric(Temp_C)) %>% select(-DOY) %>%rename(DOY=DOY2)
  #convert to dateTime
  m1g <- m1g %>% mutate(Date=as.POSIXct(strptime(paste(Year,DOY),format = ' %Y %j'),tz='US/Alaska'))
  m1g <- m1g %>% mutate(DateTime=as.POSIXct((DOY*86400)-16*60*60, 
                     origin=strptime(paste0(Year,'-01-01 00:00:00'),tz='US/Alaska',format='%Y-%m-%d %H:%M:%S'),
                     tz="US/Alaska"))
  m1g <- m1g %>% mutate(DateTimeHour=as.POSIXct(strftime(DateTime,format='%Y-%m-%d %H',tz='US/Alaska'),format='%Y-%m-%d %H',tz='US/Alaska'))
  #hourly
  m1gHr<-m1g %>% group_by(DateTimeHour,Date,Depth) %>% summarise(Temp_C=mean(Temp_C,na.rm=TRUE))
  #unix time for akima
  m1gHr <- m1gHr %>% mutate(utime=as.numeric(strftime(DateTimeHour,format = '%s')))
  m2g <- m1gHr %>% drop_na() %>% arrange(DateTimeHour,Depth)

  #interp by chunk
  mInterp<-data.frame()
  a<-length(unique(m2g$utime))
  utimes<-unique(m2g$utime)
  aSeq<-c(seq(1,a,by=100),a)
  b=1
  for(b in 1:(length(aSeq)-1)){
    starttime<-Sys.time()
    m1sub<-m2g %>% filter(utime>=utimes[aSeq[b]],utime<utimes[aSeq[b+1]],)
    #/1000 is for akima, needs <4x mag between x,y
    test<-with(m1sub,interp(x=Depth,y=utime/1000,z=Temp_C,linear = FALSE,xo=c(floor(min(Depth)):ceiling(max(Depth))),yo=unique(utime)/1000,extrap = TRUE))
    print(paste(b,": "))
    testz<-data.frame(test$z)
    names(testz)<-test$y*1000
    testz$Depth<-test$x
    testzg<-testz %>% pivot_longer(-Depth, names_to='utime',values_to='Temp_C') %>% arrange(Depth,utime)
    mInterp<-rbind(mInterp,testzg)
    print(Sys.time()-starttime)
  }
  mInterp <- mInterp %>% mutate(DateTime=as.POSIXct(as.numeric(utime),origin='1970-01-01 00:00:00 UTC',tz='US/Alaska'))  
  #rm na
  m2012wDay<-mInterp %>% drop_na()
}
#2013summer
{
  #calc depths
  m1g<-m2013s2 %>% pivot_longer(c(-Year,-DOY2,-DOY),names_to='Depth',values_to='Temp_C')
  m1g$Depth<-gsub('X','',m1g$Depth)
  m1g<-m1g %>% mutate(Depth=as.numeric(Depth),Temp_C=as.numeric(Temp_C)) %>% select(-DOY) %>%rename(DOY=DOY2)
  #convert to dateTime
  m1g <- m1g %>% mutate(Date=as.POSIXct(strptime(paste(Year,DOY),format = ' %Y %j'),tz='US/Alaska'))
  m1g <- m1g %>% mutate(DateTime=as.POSIXct((DOY*86400)-16*60*60, 
                                            origin=strptime(paste0(Year,'-01-01 00:00:00'),tz='US/Alaska',format='%Y-%m-%d %H:%M:%S'),
                                            tz="US/Alaska"))
  m1g <- m1g %>% mutate(DateTimeHour=as.POSIXct(strftime(DateTime,format='%Y-%m-%d %H',tz='US/Alaska'),format='%Y-%m-%d %H',tz='US/Alaska'))
  #hourly
  m1gHr<-m1g %>% group_by(DateTimeHour,Date,Depth) %>% summarise(Temp_C=mean(Temp_C,na.rm=TRUE))
  #unix time for akima
  m1gHr <- m1gHr %>% mutate(utime=as.numeric(strftime(DateTimeHour,format = '%s')))
  m2g <- m1gHr %>% drop_na() %>% arrange(DateTimeHour,Depth)
  
  #interp by chunk
  mInterp<-data.frame()
  a<-length(unique(m2g$utime))
  utimes<-unique(m2g$utime)
  aSeq<-c(seq(1,a,by=100),a)
  b=1
  for(b in 1:(length(aSeq)-1)){
    starttime<-Sys.time()
    m1sub<-m2g %>% filter(utime>=utimes[aSeq[b]],utime<utimes[aSeq[b+1]],)
    #/1000 is for akima, needs <4x mag between x,y
    test<-with(m1sub,interp(x=Depth,y=utime/1000,z=Temp_C,linear = FALSE,xo=c(floor(min(Depth)):ceiling(max(Depth))),yo=unique(utime)/1000,extrap = TRUE))
    print(paste(b,": "))
    testz<-data.frame(test$z)
    names(testz)<-test$y*1000
    testz$Depth<-test$x
    testzg<-testz %>% pivot_longer(-Depth, names_to='utime',values_to='Temp_C') %>% arrange(Depth,utime)
    mInterp<-rbind(mInterp,testzg)
    print(Sys.time()-starttime)
  }
  mInterp <- mInterp %>% mutate(DateTime=as.POSIXct(as.numeric(utime),origin='1970-01-01 00:00:00 UTC',tz='US/Alaska'))  
  #rm na
  m2013sDay<-mInterp %>% drop_na()
}

#2013winter
{
  #calc depths
  m1g<-m2013w %>% pivot_longer(c(-Year,-DOY2,-DOY),names_to='Depth',values_to='Temp_C')
  m1g$Depth<-gsub('X','',m1g$Depth)
  m1g<-m1g %>% mutate(Depth=as.numeric(Depth),Temp_C=as.numeric(Temp_C)) %>% select(-DOY) %>%rename(DOY=DOY2)
  #convert to dateTime
  m1g <- m1g %>% mutate(Date=as.POSIXct(strptime(paste(Year,DOY),format = ' %Y %j'),tz='US/Alaska'))
  m1g <- m1g %>% mutate(DateTime=as.POSIXct((DOY*86400)-16*60*60, 
                                            origin=strptime(paste0(Year,'-01-01 00:00:00'),tz='US/Alaska',format='%Y-%m-%d %H:%M:%S'),
                                            tz="US/Alaska"))
  m1g <- m1g %>% mutate(DateTimeHour=as.POSIXct(strftime(DateTime,format='%Y-%m-%d %H',tz='US/Alaska'),format='%Y-%m-%d %H',tz='US/Alaska'))
  #hourly
  m1gHr<-m1g %>% group_by(DateTimeHour,Date,Depth) %>% summarise(Temp_C=mean(Temp_C,na.rm=TRUE))
  #unix time for akima
  m1gHr <- m1gHr %>% mutate(utime=as.numeric(strftime(DateTimeHour,format = '%s')))
  m2g <- m1gHr %>% drop_na() %>% arrange(DateTimeHour,Depth)
  
  #interp by chunk
  mInterp<-data.frame()
  a<-length(unique(m2g$utime))
  utimes<-unique(m2g$utime)
  aSeq<-c(seq(1,a,by=100),a)
  b=1
  for(b in 1:(length(aSeq)-1)){
    starttime<-Sys.time()
    m1sub<-m2g %>% filter(utime>=utimes[aSeq[b]],utime<utimes[aSeq[b+1]],)
    #/1000 is for akima, needs <4x mag between x,y
    test<-with(m1sub,interp(x=Depth,y=utime/1000,z=Temp_C,linear = FALSE,xo=c(floor(min(Depth)):ceiling(max(Depth))),yo=unique(utime)/1000,extrap = TRUE))
    print(paste(b,": "))
    testz<-data.frame(test$z)
    names(testz)<-test$y*1000
    testz$Depth<-test$x
    testzg<-testz %>% pivot_longer(-Depth, names_to='utime',values_to='Temp_C') %>% arrange(Depth,utime)
    mInterp<-rbind(mInterp,testzg)
    print(Sys.time()-starttime)
  }
  mInterp <- mInterp %>% mutate(DateTime=as.POSIXct(as.numeric(utime),origin='1970-01-01 00:00:00 UTC',tz='US/Alaska'))  
  #rm na
  m2013wDay<-mInterp %>% drop_na()
}
#2014summer
{
  #calc depths
  m1g<-m2014s2 %>% pivot_longer(c(-Year,-DOY2,-DOY),names_to='Depth',values_to='Temp_C')
  m1g$Depth<-gsub('X','',m1g$Depth)
  m1g<-m1g %>% mutate(Depth=as.numeric(Depth),Temp_C=as.numeric(Temp_C)) %>% select(-DOY) %>%rename(DOY=DOY2)
  #convert to dateTime
  m1g <- m1g %>% mutate(Date=as.POSIXct(strptime(paste(Year,DOY),format = ' %Y %j'),tz='US/Alaska'))
  m1g <- m1g %>% mutate(DateTime=as.POSIXct((DOY*86400)-16*60*60, 
                                            origin=strptime(paste0(Year,'-01-01 00:00:00'),tz='US/Alaska',format='%Y-%m-%d %H:%M:%S'),
                                            tz="US/Alaska"))
  m1g <- m1g %>% mutate(DateTimeHour=as.POSIXct(strftime(DateTime,format='%Y-%m-%d %H',tz='US/Alaska'),format='%Y-%m-%d %H',tz='US/Alaska'))
  #hourly
  m1gHr<-m1g %>% group_by(DateTimeHour,Date,Depth) %>% summarise(Temp_C=mean(Temp_C,na.rm=TRUE))
  #unix time for akima
  m1gHr <- m1gHr %>% mutate(utime=as.numeric(strftime(DateTimeHour,format = '%s')))
  m2g <- m1gHr %>% drop_na() %>% arrange(DateTimeHour,Depth)
  
  #interp by chunk
  mInterp<-data.frame()
  a<-length(unique(m2g$utime))
  utimes<-unique(m2g$utime)
  aSeq<-c(seq(1,a,by=100),a)
  b=1
  for(b in 1:(length(aSeq)-1)){
    starttime<-Sys.time()
    m1sub<-m2g %>% filter(utime>=utimes[aSeq[b]],utime<utimes[aSeq[b+1]],)
    #/1000 is for akima, needs <4x mag between x,y
    test<-with(m1sub,interp(x=Depth,y=utime/1000,z=Temp_C,linear = FALSE,xo=c(floor(min(Depth)):ceiling(max(Depth))),yo=unique(utime)/1000,extrap = TRUE))
    print(paste(b,": "))
    testz<-data.frame(test$z)
    names(testz)<-test$y*1000
    testz$Depth<-test$x
    testzg<-testz %>% pivot_longer(-Depth, names_to='utime',values_to='Temp_C') %>% arrange(Depth,utime)
    mInterp<-rbind(mInterp,testzg)
    print(Sys.time()-starttime)
  }
  mInterp <- mInterp %>% mutate(DateTime=as.POSIXct(as.numeric(utime),origin='1970-01-01 00:00:00 UTC',tz='US/Alaska'))  
  #rm na
  m2014sDay<-mInterp %>% drop_na()
  #fix errors
  m2014sDay %>% filter(Temp_C>(50))
  m2014sDay[m2014sDay$Temp_C>(50),]$Temp_C<-NA
  m2014sDay<-m2014sDay %>% drop_na()
  m2014sDay %>% filter(Temp_C<(-50))
  class(m2014sDay)
  m2014sDay<-as.data.frame(m2014sDay)
  m2014sDay<-m2014sDay %>% drop_na()
  m2014sDay[m2014sDay$Temp_C<(-50),]$Temp_C<-NA
  m2014sDay %>% filter(Temp_C<(1))
  m2014sDay<-m2014sDay %>% drop_na()
  m2014sDay[m2014sDay$Temp_C<(1),]$Temp_C<-NA
  m2014sDay<-m2014sDay %>% drop_na()
  m2014sDay[m2014sDay$Temp_C>(19),]$Temp_C<-NA
  m2014sDay<-m2014sDay %>% drop_na()
}

#2014winter
{
  #calc depths
  m1g<-m2014w %>% pivot_longer(c(-Year,-DOY2,-DOY),names_to='Depth',values_to='Temp_C')
  m1g$Depth<-gsub('X','',m1g$Depth)
  m1g<-m1g %>% mutate(Depth=as.numeric(Depth),Temp_C=as.numeric(Temp_C)) %>% select(-DOY) %>%rename(DOY=DOY2)
  #convert to dateTime
  m1g <- m1g %>% mutate(Date=as.POSIXct(strptime(paste(Year,DOY),format = ' %Y %j'),tz='US/Alaska'))
  m1g <- m1g %>% mutate(DateTime=as.POSIXct((DOY*86400)-16*60*60, 
                                            origin=strptime(paste0(Year,'-01-01 00:00:00'),tz='US/Alaska',format='%Y-%m-%d %H:%M:%S'),
                                            tz="US/Alaska"))
  m1g <- m1g %>% mutate(DateTimeHour=as.POSIXct(strftime(DateTime,format='%Y-%m-%d %H',tz='US/Alaska'),format='%Y-%m-%d %H',tz='US/Alaska'))
  #hourly
  m1gHr<-m1g %>% group_by(DateTimeHour,Date,Depth) %>% summarise(Temp_C=mean(Temp_C,na.rm=TRUE))
  #unix time for akima
  m1gHr <- m1gHr %>% mutate(utime=as.numeric(strftime(DateTimeHour,format = '%s')))
  m2g <- m1gHr %>% drop_na() %>% arrange(DateTimeHour,Depth)
  
  #interp by chunk
  mInterp<-data.frame()
  a<-length(unique(m2g$utime))
  utimes<-unique(m2g$utime)
  aSeq<-c(seq(1,a,by=100),a)
  b=1
  for(b in 1:(length(aSeq)-1)){
    starttime<-Sys.time()
    m1sub<-m2g %>% filter(utime>=utimes[aSeq[b]],utime<utimes[aSeq[b+1]],)
    #/1000 is for akima, needs <4x mag between x,y
    test<-with(m1sub,interp(x=Depth,y=utime/1000,z=Temp_C,linear = FALSE,xo=c(floor(min(Depth)):ceiling(max(Depth))),yo=unique(utime)/1000,extrap = TRUE))
    print(paste(b,": "))
    testz<-data.frame(test$z)
    names(testz)<-test$y*1000
    testz$Depth<-test$x
    testzg<-testz %>% pivot_longer(-Depth, names_to='utime',values_to='Temp_C') %>% arrange(Depth,utime)
    mInterp<-rbind(mInterp,testzg)
    print(Sys.time()-starttime)
  }
  mInterp <- mInterp %>% mutate(DateTime=as.POSIXct(as.numeric(utime),origin='1970-01-01 00:00:00 UTC',tz='US/Alaska'))  
  #rm na
  m2014wDay<-mInterp %>% drop_na()
  
  m2014wDay[m2014wDay$Temp_C<(0),]$Temp_C<-NA
  m2014wDay<-m2014wDay %>% drop_na()
}
#2015summer
{
  #calc depths
  m1g<-m2015s2 %>% pivot_longer(c(-Year,-DOY2,-DOY),names_to='Depth',values_to='Temp_C')
  m1g$Depth<-gsub('X','',m1g$Depth)
  m1g<-m1g %>% mutate(Depth=as.numeric(Depth),Temp_C=as.numeric(Temp_C)) %>% select(-DOY) %>%rename(DOY=DOY2)
  #convert to dateTime
  m1g <- m1g %>% mutate(Date=as.POSIXct(strptime(paste(Year,DOY),format = ' %Y %j'),tz='US/Alaska'))
  m1g <- m1g %>% mutate(DateTime=as.POSIXct((DOY*86400)-16*60*60, 
                                            origin=strptime(paste0(Year,'-01-01 00:00:00'),tz='US/Alaska',format='%Y-%m-%d %H:%M:%S'),
                                            tz="US/Alaska"))
  m1g <- m1g %>% mutate(DateTimeHour=as.POSIXct(strftime(DateTime,format='%Y-%m-%d %H',tz='US/Alaska'),format='%Y-%m-%d %H',tz='US/Alaska'))
  #hourly
  m1gHr<-m1g %>% group_by(DateTimeHour,Date,Depth) %>% summarise(Temp_C=mean(Temp_C,na.rm=TRUE))
  #unix time for akima
  m1gHr <- m1gHr %>% mutate(utime=as.numeric(strftime(DateTimeHour,format = '%s')))
  m2g <- m1gHr %>% drop_na() %>% arrange(DateTimeHour,Depth)
  
  #interp by chunk
  mInterp<-data.frame()
  a<-length(unique(m2g$utime))
  utimes<-unique(m2g$utime)
  aSeq<-c(seq(1,a,by=100),a)
  b=1
  for(b in 1:(length(aSeq)-1)){
    starttime<-Sys.time()
    m1sub<-m2g %>% filter(utime>=utimes[aSeq[b]],utime<utimes[aSeq[b+1]],)
    #/1000 is for akima, needs <4x mag between x,y
    test<-with(m1sub,interp(x=Depth,y=utime/1000,z=Temp_C,linear = FALSE,xo=c(floor(min(Depth)):ceiling(max(Depth))),yo=unique(utime)/1000,extrap = TRUE))
    print(paste(b,": "))
    testz<-data.frame(test$z)
    names(testz)<-test$y*1000
    testz$Depth<-test$x
    testzg<-testz %>% pivot_longer(-Depth, names_to='utime',values_to='Temp_C') %>% arrange(Depth,utime)
    mInterp<-rbind(mInterp,testzg)
    print(Sys.time()-starttime)
  }
  mInterp <- mInterp %>% mutate(DateTime=as.POSIXct(as.numeric(utime),origin='1970-01-01 00:00:00 UTC',tz='US/Alaska'))  
  #rm na
  m2015sDay<-mInterp %>% drop_na()
  #fix errors
  m2015sDay %>% filter(Temp_C<(0))
  m2015sDay[m2015sDay$Temp_C<(0),]$Temp_C<-NA
  m2015sDay<-m2015sDay %>% drop_na()
}


#2015winter
{
  #calc depths
  m1g<-m2015w %>% pivot_longer(c(-Year,-DOY2,-DOY),names_to='Depth',values_to='Temp_C')
  m1g$Depth<-gsub('X','',m1g$Depth)
  m1g<-m1g %>% mutate(Depth=as.numeric(Depth),Temp_C=as.numeric(Temp_C)) %>% select(-DOY) %>%rename(DOY=DOY2)
  #convert to dateTime
  m1g <- m1g %>% mutate(Date=as.POSIXct(strptime(paste(Year,DOY),format = ' %Y %j'),tz='US/Alaska'))
  m1g <- m1g %>% mutate(DateTime=as.POSIXct((DOY*86400)-16*60*60, 
                                            origin=strptime(paste0(Year,'-01-01 00:00:00'),tz='US/Alaska',format='%Y-%m-%d %H:%M:%S'),
                                            tz="US/Alaska"))
  m1g <- m1g %>% mutate(DateTimeHour=as.POSIXct(strftime(DateTime,format='%Y-%m-%d %H',tz='US/Alaska'),format='%Y-%m-%d %H',tz='US/Alaska'))
  #hourly
  m1gHr<-m1g %>% group_by(DateTimeHour,Date,Depth) %>% summarise(Temp_C=mean(Temp_C,na.rm=TRUE))
  #unix time for akima
  m1gHr <- m1gHr %>% mutate(utime=as.numeric(strftime(DateTimeHour,format = '%s')))
  m2g <- m1gHr %>% drop_na() %>% arrange(DateTimeHour,Depth)
  
  #interp by chunk
  mInterp<-data.frame()
  a<-length(unique(m2g$utime))
  utimes<-unique(m2g$utime)
  aSeq<-c(seq(1,a,by=100),a)
  b=1
  for(b in 1:(length(aSeq)-1)){
    starttime<-Sys.time()
    m1sub<-m2g %>% filter(utime>=utimes[aSeq[b]],utime<utimes[aSeq[b+1]],)
    #/1000 is for akima, needs <4x mag between x,y
    test<-with(m1sub,interp(x=Depth,y=utime/1000,z=Temp_C,linear = FALSE,xo=c(floor(min(Depth)):ceiling(max(Depth))),yo=unique(utime)/1000,extrap = TRUE))
    print(paste(b,": "))
    testz<-data.frame(test$z)
    names(testz)<-test$y*1000
    testz$Depth<-test$x
    testzg<-testz %>% pivot_longer(-Depth, names_to='utime',values_to='Temp_C') %>% arrange(Depth,utime)
    mInterp<-rbind(mInterp,testzg)
    print(Sys.time()-starttime)
  }
  mInterp <- mInterp %>% mutate(DateTime=as.POSIXct(as.numeric(utime),origin='1970-01-01 00:00:00 UTC',tz='US/Alaska'))  
  #rm na
  m2015wDay<-mInterp %>% drop_na()
  #remove bad temps
  m2015wDay$Temp_C[m2015wDay$Temp_C<(-1)]<-NA
  #rm na
  m2015wDay<-m2015wDay %>% drop_na()  
  
  m2015wDay[m2015wDay$Temp_C<(0),]$Temp_C<-NA
  m2015wDay<-m2015wDay %>% drop_na()
}
#2016summer
{
  #calc depths
  m1g<-m2016s2 %>% pivot_longer(c(-Year,-DOY2,-DOY),names_to='Depth',values_to='Temp_C')
  m1g$Depth<-gsub('X','',m1g$Depth)
  m1g<-m1g %>% mutate(Depth=as.numeric(Depth),Temp_C=as.numeric(Temp_C)) %>% select(-DOY) %>%rename(DOY=DOY2)
  #convert to dateTime
  m1g <- m1g %>% mutate(Date=as.POSIXct(strptime(paste(Year,DOY),format = ' %Y %j'),tz='US/Alaska'))
  m1g <- m1g %>% mutate(DateTime=as.POSIXct((DOY*86400)-16*60*60, 
                                            origin=strptime(paste0(Year,'-01-01 00:00:00'),tz='US/Alaska',format='%Y-%m-%d %H:%M:%S'),
                                            tz="US/Alaska"))
  m1g <- m1g %>% mutate(DateTimeHour=as.POSIXct(strftime(DateTime,format='%Y-%m-%d %H',tz='US/Alaska'),format='%Y-%m-%d %H',tz='US/Alaska'))
  #hourly
  m1gHr<-m1g %>% group_by(DateTimeHour,Date,Depth) %>% summarise(Temp_C=mean(Temp_C,na.rm=TRUE))
  #unix time for akima
  m1gHr <- m1gHr %>% mutate(utime=as.numeric(strftime(DateTimeHour,format = '%s')))
  m2g <- m1gHr %>% drop_na() %>% arrange(DateTimeHour,Depth)
  
  #interp by chunk
  mInterp<-data.frame()
  a<-length(unique(m2g$utime))
  utimes<-unique(m2g$utime)
  aSeq<-c(seq(1,a,by=100),a)
  b=1
  for(b in 1:(length(aSeq)-1)){
    starttime<-Sys.time()
    m1sub<-m2g %>% filter(utime>=utimes[aSeq[b]],utime<utimes[aSeq[b+1]],)
    #/1000 is for akima, needs <4x mag between x,y
    test<-with(m1sub,interp(x=Depth,y=utime/1000,z=Temp_C,linear = FALSE,xo=c(floor(min(Depth)):ceiling(max(Depth))),yo=unique(utime)/1000,extrap = TRUE))
    print(paste(b,": "))
    testz<-data.frame(test$z)
    names(testz)<-test$y*1000
    testz$Depth<-test$x
    testzg<-testz %>% pivot_longer(-Depth, names_to='utime',values_to='Temp_C') %>% arrange(Depth,utime)
    mInterp<-rbind(mInterp,testzg)
    print(Sys.time()-starttime)
  }
  mInterp <- mInterp %>% mutate(DateTime=as.POSIXct(as.numeric(utime),origin='1970-01-01 00:00:00 UTC',tz='US/Alaska'))  
  #rm na
  m2016sDay<-mInterp %>% drop_na()
}

##...............................................................
##               write interpolated data to file               ..
##...............................................................

write.csv(m2012wDay,file=paste0(data_dir,'Toolik_Lake_Temp_2012wInterpolate.csv'),row.names = FALSE)
write.csv(m2013sDay,file=paste0(data_dir,'Toolik_Lake_Temp_2013sInterpolate.csv'),row.names = FALSE)
write.csv(m2013wDay,file=paste0(data_dir,'Toolik_Lake_Temp_2013wInterpolate.csv'),row.names = FALSE)
write.csv(m2014sDay,file=paste0(data_dir,'Toolik_Lake_Temp_2014sInterpolate.csv'),row.names = FALSE)
write.csv(m2014wDay,file=paste0(data_dir,'Toolik_Lake_Temp_2014wInterpolate.csv'),row.names = FALSE)
write.csv(m2015sDay,file=paste0(data_dir,'Toolik_Lake_Temp_2015sInterpolate.csv'),row.names = FALSE)
write.csv(m2015wDay,file=paste0(data_dir,'Toolik_Lake_Temp_2015wInterpolate.csv'),row.names = FALSE)
write.csv(m2016sDay,file=paste0(data_dir,'Toolik_Lake_Temp_2016sInterpolate.csv'),row.names = FALSE)


##...............................................................
##                           Combine                           ..
##                       write to file                         ..
##...............................................................

m1combAll<-rbind(m2012wDay,
              m2013sDay,
              m2013wDay,
              m2014sDay,
              m2014wDay,
              m2015sDay,
              m2015wDay,
              m2016sDay)
write.csv(m1combAll,file=paste(data_dir,'Toolik_Lake_Temp_2012_2016Allinterpolate.csv'),row.names = FALSE)

