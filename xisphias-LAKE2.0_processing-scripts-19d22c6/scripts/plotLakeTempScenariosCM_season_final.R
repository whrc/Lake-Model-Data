##################################################################
##                Calculate errors for scenarios with seasons   ##
##          Combine errors for seasonal barplot                 ##
##                      Jason Clark                             ##
##                   jaclark2@alaska.edu                        ##
##                          2021-10-22                          ##
##################################################################
library(tidyverse)
library(gtable)
library(grid)
library(gridExtra)
library(ggh4x)
library(ggpubr)
library(scales)
##---------------------------------------------------------------
##               Function for calculating errors with seasons  --
##---------------------------------------------------------------
calcErrorScenarioSeason<-function(obsDataFile, modDataDir, modWxFile, Scenarios=c(''), scenDataDir,Depths=c(1,2),Season='Season_Ice'){
  starttime<-Sys.time()
  library(tidyverse)
  library(Metrics)
  measured<-read.csv(obsDataFile)
  modWater =read.table(paste0(modDataDir,'/time_series/water_temp  1  1.dat'),skip=6)
  modSoil  =read.table(paste0(modDataDir,'/time_series/soil_temp  1  1.dat'),skip=6)
  modLayers=read.table(paste0(modDataDir,'/time_series/layers  1  1.dat'),skip=19)
  modWX    =read.csv(modWxFile,header=FALSE)
  print('Data read')
  # prep data
  #measured
  measured$DateTime<-as.POSIXct(measured$DateTime,format="%Y-%m-%d %H:%M:%S",tz="US/Alaska")
  print(paste('Observed Water Depths:',paste(unique(measured$Depth),collapse = ', ')))
  if(length(grep('Atgasuk',obsDataFile))==1){
    #fix depth for atq, depth was recorded as 3.25 in datafile, but max depth of lake is 2.5m
    measured$Depth[measured$Depth==3.25]<-2.5
  }

  #wx
  #measured
  names(modWX)<-c('DateTime','T_K','Pres','LWinstrument','SWinstrument','U','V','H','P')
  # head(modWX$DateTime)
  modWX$DateTime<-as.POSIXct(modWX$DateTime,format="%Y-%m-%d %H:%M:%S",tz="US/Alaska")
  #calc season from daily AT
  modWX<-modWX %>% mutate(Date=as.POSIXct(format(DateTime,format='%Y-%m-%d'),tz='US/Alaska')) %>% 
    group_by(Date) %>% 
    mutate(meanAT=mean(T_K,na.rm = TRUE)) %>%
    mutate(Season_AT=factor(ifelse(meanAT>273.15,'Thawed','Frozen')))
  
  #prep model data
  names(modWater)[1:5]<-c('Year','Month','Day','Hour','integrationTime')
  names(modWater)[seq(6,dim(modWater)[2],by=2)]<-paste0('Depth',seq(1,(dim(modWater)[2]-5)/2,by=1))
  names(modWater)[seq(7,dim(modWater)[2],by=2)]<-paste0('Temp',seq(1,(dim(modWater)[2]-5)/2,by=1))
  modWater$DateTime<-with(modWater,as.POSIXct(strptime(paste0(Year,'-',Month,'-',Day,'H',Hour),format = '%Y-%m-%dH%H'),tz='US/Alaska'))
  print('model dates:');print(min(na.omit(modWater$DateTime)));print(max(na.omit(modWater$DateTime)))
  m1gD<-gather(modWater %>% select(-(starts_with('Temp'))),colName,Depth,-Year,-Month,-Day,-Hour,-integrationTime,-DateTime)
  m1gT<-gather(modWater %>% select(-(starts_with('Depth'))),colName,Temp,-Year,-Month,-Day,-Hour,-integrationTime,-DateTime)
  calcm<-cbind(m1gD %>% select(-colName) %>% mutate(Depth=factor(Depth)),Temp=m1gT %>% pull(Temp))
  #calc season from daily shallow modeled WT
  minDepth=min(as.numeric(as.character(calcm$Depth)))
  calcm_season<-calcm %>% 
    filter(Depth==minDepth) %>%
    mutate(Date=as.POSIXct(format(DateTime,format='%Y-%m-%d'),tz='US/Alaska')) %>%
    group_by(Date) %>%
    mutate(meanWT=mean(Temp,na.rm = TRUE)) %>%
    mutate(Season_WT=factor(ifelse(meanWT>0,'Thawed','Frozen')))

  #prep soil
  names(modSoil)[1:5]<-c('Year','Month','Day','Hour','integrationTime')
  names(modSoil)[seq(6,dim(modSoil)[2],by=2)]<-paste0('Depth',seq(1,(dim(modSoil)[2]-5)/2,by=1))
  names(modSoil)[seq(7,dim(modSoil)[2],by=2)]<-paste0('Temp',seq(1,(dim(modSoil)[2]-5)/2,by=1))
  modSoil$DateTime<-with(modSoil,as.POSIXct(strptime(paste0(Year,'-',Month,'-',Day,'H',Hour),format = '%Y-%m-%dH%H'),tz='US/Alaska'))
  modSoilgD<-gather(modSoil %>% select(-(starts_with('Temp'))),colName,Depth,-Year,-Month,-Day,-Hour,-integrationTime,-DateTime)
  modSoilgT<-gather(modSoil %>% select(-(starts_with('Depth'))),colName,Temp,-Year,-Month,-Day,-Hour,-integrationTime,-DateTime)
  Soilm<-cbind(modSoilgD %>% select(-colName) %>% mutate(Depth=factor(Depth)),Temp=modSoilgT %>% pull(Temp))
  # head(Soilm);levels(Soilm$Depth)
  
  #prep layers
  names(modLayers)<-c("year","month","day","hour","intTimeHr","waterlayer_m","Wmixedlayerthickness_m","Emixedlayerthickness_m","Smixedlayerthickness_m","Nmixedlayerthickness_m","Wlowerlayerthickness_m","Elowerlayerthickness_m","Slowerlayerthickness_m","Nlowerlayerthickness_m","icelayer_m","snowlayer_m","bottomice_m","reservoirvolume_m3","volumedeficit_accumulated_m3")
  modLayers$DateTime<-as.POSIXct(strptime(with(modLayers,paste0(year,'-',month,'-',day,' ',as.integer(hour),':00:00')),format='%Y-%m-%d %H:%M:%S'),tz='US/Alaska')
  modLayers$Date<-as.POSIXct(format(modLayers$DateTime,format='%Y-%m-%d'),tz='US/Alaska')
  # head(modLayers)
  #calc season from daily ice thickness
  modLayers_season<-modLayers %>% 
    mutate(Date=as.POSIXct(format(DateTime,format='%Y-%m-%d'),tz='US/Alaska')) %>%
    group_by(Date) %>%
    select(DateTime,Date,icelayer_m) %>%
    mutate(meanIce=mean(icelayer_m,na.rm = TRUE)) %>%
    mutate(Season_Ice=factor(ifelse(meanIce==0,'Thawed','Frozen')))
  
  cat(paste(paste('Ice thaw days:',modLayers_season %>% select(Date,meanIce,Season_Ice) %>% unique() %>% mutate(WT_bi=as.numeric(Season_Ice)-1) %>% pull(WT_bi) %>% sum()),
    paste('WT thaw days:',calcm_season %>% select(Date,meanWT,Season_WT) %>% unique() %>% mutate(WT_bi=as.numeric(Season_WT)-1) %>% pull(WT_bi) %>% sum()),
    paste('AT thaw days:',modWX %>% select(Date,meanAT,Season_AT) %>% unique() %>% mutate(WT_bi=as.numeric(Season_AT)-1) %>% pull(WT_bi) %>% sum()),
    sep='\n'
  ))
  
  #error catching, foxden is daily wx, but hourly output need to match by date and use datetime from model
  if(length(modWX$DateTime[modWX$DateTime %in% calcm_season$DateTime])>0){
    #no error, dates overlap merge will work
    seasonDat<-merge(
                    merge(modWX %>% select(DateTime,Date,Season_AT),
                          calcm_season %>% select(DateTime,Season_WT), 
                          all.x = TRUE),
                    modLayers_season %>% select(DateTime,Season_Ice),
                    all.x=TRUE) %>% tibble()
  }else{
    #error, no matching dates
    seasonDat<-merge(
                      merge(modWX %>% select(Date,Season_AT),
                            unique(calcm_season %>% select(Date,DateTime,Season_WT)), 
                            all.y = TRUE,by='Date'),
                      modLayers_season %>% select(DateTime,Date,Season_Ice),
                      all.x=TRUE,by=c('DateTime','Date')) %>% tibble()
  }
  
  print('Data prepped')
 
  a=1
  nScen=length(Scenarios)
  
  scenWater<-data.frame(Scenario=character(), Year=numeric(), Month=numeric(), Day=numeric(), Hour=numeric(), integrationTime=numeric(), DateTime=character(), Depth=numeric(), Temp=numeric()
                        )
  scenSoil<-data.frame(Scenario=character(), Year=numeric(), Month=numeric(), Day=numeric(), Hour=numeric(), integrationTime=numeric(), DateTime=character(), Depth=numeric(), Temp=numeric(),Date=character()
                        )
  scenLayers<-data.frame(Scenario=character(), Year=numeric(), Month=numeric(), Day=numeric(), Hour=numeric(), intTimeHr=numeric(), DateTime=character(), Date=character(),icelayer_m=numeric(), snowlayer_m=numeric()
                        )
  for(a in 1:nScen){
    print(Scenarios[a])
    modWaterScen =read.table(paste0(scenDataDir,'/',Scenarios[a],'/time_series/water_temp  1  1.dat'),skip=6)
    modSoilScen  =read.table(paste0(scenDataDir,'/',Scenarios[a],'/time_series/soil_temp  1  1.dat'),skip=6)
    modLayersScen=read.table(paste0(scenDataDir,'/',Scenarios[a],'/time_series/layers  1  1.dat'),skip=19)
    
    #prep model water temp
    names(modWaterScen)[1:5]<-c('Year','Month','Day','Hour','integrationTime')
    names(modWaterScen)[seq(6,dim(modWaterScen)[2],by=2)]<-paste0('Depth',seq(1,(dim(modWaterScen)[2]-5)/2,by=1))
    names(modWaterScen)[seq(7,dim(modWaterScen)[2],by=2)]<-paste0('Temp',seq(1,(dim(modWaterScen)[2]-5)/2,by=1))
    modWaterScen$DateTime<-with(modWaterScen,as.POSIXct(strptime(paste0(Year,'-',Month,'-',Day,'H',Hour),format = '%Y-%m-%dH%H'),tz='US/Alaska'))
    m1gD<-gather(modWaterScen %>% select(-(starts_with('Temp'))),colName,Depth,-Year,-Month,-Day,-Hour,-integrationTime,-DateTime)
    m1gT<-gather(modWaterScen %>% select(-(starts_with('Depth'))),colName,Temp,-Year,-Month,-Day,-Hour,-integrationTime,-DateTime)
    calcmScen<-cbind(m1gD %>% select(-colName) %>% mutate(Depth=factor(Depth)),Temp=m1gT %>% pull(Temp)) %>%
      filter(Depth%in%Depths)
    scenWater<-rbind(scenWater,data.frame(Scenario=Scenarios[a],calcmScen))
    print(paste('water levels: ',paste(levels(scenWater$Depth),collapse= ', ')))
    #prep soil
    names(modSoilScen)[1:5]<-c('Year','Month','Day','Hour','integrationTime')
    names(modSoilScen)[seq(6,dim(modSoilScen)[2],by=2)]<-paste0('Depth',seq(1,(dim(modSoilScen)[2]-5)/2,by=1))
    names(modSoilScen)[seq(7,dim(modSoilScen)[2],by=2)]<-paste0('Temp',seq(1,(dim(modSoilScen)[2]-5)/2,by=1))
    modSoilScen$DateTime<-with(modSoilScen,as.POSIXct(strptime(paste0(Year,'-',Month,'-',Day,'H',Hour),format = '%Y-%m-%dH%H'),tz='US/Alaska'))
    modSoilScengD<-gather(modSoilScen %>% select(-(starts_with('Temp'))),colName,Depth,-Year,-Month,-Day,-Hour,-integrationTime,-DateTime)
    modSoilScengT<-gather(modSoilScen %>% select(-(starts_with('Depth'))),colName,Temp,-Year,-Month,-Day,-Hour,-integrationTime,-DateTime)
    SoilmScen<-cbind(modSoilScengD %>% select(-colName) %>% mutate(Depth=factor(Depth)),Temp=modSoilScengT %>% pull(Temp))
    # head(Soilm);levels(Soilm$Depth)
    scenSoil<-rbind(scenSoil,data.frame(Scenario=Scenarios[a],SoilmScen))
    
    #prep layers
    names(modLayersScen)<-c("Year","Month","Day","Hour","intTimeHr","waterlayer_m","Wmixedlayerthickness_m","Emixedlayerthickness_m","Smixedlayerthickness_m","Nmixedlayerthickness_m","Wlowerlayerthickness_m","Elowerlayerthickness_m","Slowerlayerthickness_m","Nlowerlayerthickness_m","icelayer_m","snowlayer_m","bottomice_m","reservoirvolume_m3","volumedeficit_accumulated_m3")
    modLayersScen$DateTime<-as.POSIXct(strptime(with(modLayersScen,paste0(Year,'-',Month,'-',Day,' ',as.integer(Hour),':00:00')),format='%Y-%m-%d %H:%M:%S'),tz='US/Alaska')
    modLayersScen$Date<-as.POSIXct(format(modLayersScen$DateTime,format='%Y-%m-%d'),tz='US/Alaska')
    # head(modLayersScen)
    scenLayers<-rbind(scenLayers,data.frame(Scenario=Scenarios[a],modLayersScen %>% 
                                              select(Year, Month, Day, Hour, intTimeHr,DateTime, Date, icelayer_m, snowlayer_m)
    ))
  }
 
  print("Scenarios prepped")
  
  #calculate error from baseline
  {
    print('calculating errors...')
    #water temp
    thisDate3<-min(calcm$DateTime,na.rm = TRUE)
    thisDate4<-max(calcm$DateTime,na.rm = TRUE)
    #combine meausured, baseline
    dfHrly<-merge(data.frame(measured %>% filter(Depth %in% Depths) %>% filter(DateTime<=thisDate4,DateTime>thisDate3) %>% rename(temp_c=Temp_C)  %>% select(DateTime,Depth,observed=temp_c)),
                  data.frame( (calcm %>% filter(Depth %in% Depths) %>% select(DateTime,Depth,modeled=Temp) )) ,by=c('DateTime','Depth'),all=TRUE)
    str(dfHrly)  
    #combine scenarios
    newdat<-rbind(dfHrly %>% mutate(Scenario='Baseline',DateTime=as.character(DateTime),Temp=as.numeric(modeled)) %>% select(Scenario,DateTime,Depth,Temp) ,
                  scenWater %>% mutate(DateTime=as.character(DateTime)) %>% select(Scenario,DateTime,Depth,Temp)
    )
    # str(newdat)
    newdat$Temp<-as.numeric(newdat$Temp)
    newdat$Depth<-factor(newdat$Depth)
    # need to drop NA and remove duplicate records
    newdat_W <- newdat %>% drop_na() %>% unique() %>% pivot_wider(id_col=c(DateTime,Depth),names_from=Scenario,values_from=c(Temp),names_repair='unique')
    head(newdat_W$DateTime)
    head( seasonDat %>% select(-Date) %>% mutate(DateTime=as.character(DateTime)) %>% pull(DateTime))
    #add seasons
    newdat_W <- merge(newdat_W,
                      seasonDat %>% select(-Date) %>% mutate(DateTime=as.character(DateTime)),
                      all.x = TRUE,by='DateTime') %>% tibble()
    newDat_L<-newdat_W %>% pivot_longer(cols=starts_with('S_'),names_to = 'Scenario',values_to='Temp') %>% mutate(Scenario=factor(Scenario))
    
    (waterErr<-newDat_L %>% drop_na() %>% group_by(Scenario,Depth,.dots=setNames(nm=Season)) %>% summarise(MAE=mae(Baseline,Temp), #mean abs error
                                                                                 Bias=bias(Baseline,Temp), #mean ave error
                                                                                 MSE=mse(Baseline,Temp),   #mean sq error
                                                                                 PB=percent_bias(Baseline,Temp), #percent bias
                                                                                 RAE=rae(Baseline,Temp),      #relative abs error
                                                                                 RMSE=rmse(Baseline,Temp),
                                                                                 cumAE=sum(ae(Baseline,Temp)),
                                                                                 zScore_mean=mean((Temp-mean(Baseline))/sd(Baseline)),
                                                                                 zScore_med=median((Temp-mean(Baseline))/sd(Baseline))
    )
    )
    #soil temp, need to specify depth of interest
    #ice thickness
    newdat<-rbind(modLayers %>% mutate(Scenario='Baseline',DateTime=as.character(DateTime),Depth=as.numeric(icelayer_m)) %>% select(Scenario,DateTime,Depth) ,
                  scenLayers %>% mutate(DateTime=as.character(DateTime),Depth=as.numeric(icelayer_m)) %>% select(Scenario,DateTime,Depth)
    )
    newdat$Depth<-as.numeric(newdat$Depth)
    newdat_W <- newdat %>% drop_na() %>% pivot_wider(id_col=DateTime,names_from=Scenario,values_from=c(Depth),names_repair='unique')
    #add seasons
    newdat_W <- merge(newdat_W,
                      seasonDat %>% select(-Date) %>% mutate(DateTime=as.character(DateTime)),
                      all.x = TRUE,by='DateTime') %>% tibble()
    # summary(newdat_W)
    newDat_L<-newdat_W %>% pivot_longer(cols=starts_with('S_'),names_to = 'Scenario',values_to='Depth') %>% mutate(Scenario=factor(Scenario))
    (iceErr<-newDat_L %>% filter(Baseline>0) %>% group_by(Scenario,.dots=setNames(nm=Season))  %>% summarise(MAE=mae(Baseline,        Depth), #mean abs error
                                                                                   Bias=bias(Baseline,      Depth), #mean ave error
                                                                                   MSE=mse(Baseline,        Depth),   #mean sq error
                                                                                   PB=percent_bias(Baseline,Depth), #percent bias
                                                                                   RAE=rae(Baseline,        Depth),      #relative abs error
                                                                                   RMSE=rmse(Baseline,      Depth),
                                                                                   cumAE=sum(ae(Baseline,   Depth)),
                                                                                   zScore_mean=mean((Depth-mean(Baseline))/sd(Baseline)),
                                                                                   zScore_med=median((Depth-mean(Baseline))/sd(Baseline))
    )
    )
    #snow thickness
    newdat<-rbind(modLayers %>% mutate(Scenario='Baseline',DateTime=as.character(DateTime),Depth=as.numeric(snowlayer_m)) %>% select(Scenario,DateTime,Depth) ,
                  scenLayers %>% mutate(DateTime=as.character(DateTime),Depth=as.numeric(snowlayer_m)) %>% select(Scenario,DateTime,Depth)
    )
    newdat$Depth<-as.numeric(newdat$Depth)
    newdat_W <- newdat %>% drop_na() %>% pivot_wider(id_col=DateTime,names_from=Scenario,values_from=c(Depth),names_repair='unique')
    #add seasons
    newdat_W <- merge(newdat_W,
                      seasonDat %>% select(-Date) %>% mutate(DateTime=as.character(DateTime)),
                      all.x = TRUE,by='DateTime') %>% tibble()
    # summary(newdat_W)
    newDat_L<-newdat_W %>% pivot_longer(cols=starts_with('S_'),names_to = 'Scenario',values_to='Depth') %>% mutate(Scenario=factor(Scenario))
    # str(newDat_L)
    (snowErr<-newDat_L %>% filter(Baseline>0) %>% group_by(Scenario,.dots=setNames(nm=Season))  %>% summarise(MAE=mae(Baseline,        Depth), #mean abs error
                                                                                    Bias=bias(Baseline,      Depth), #mean ave error
                                                                                    MSE=mse(Baseline,        Depth),   #mean sq error
                                                                                    PB=percent_bias(Baseline,Depth), #percent bias
                                                                                    RAE=rae(Baseline,        Depth),      #relative abs error
                                                                                    RMSE=rmse(Baseline,      Depth),
                                                                                    cumAE=sum(ae(Baseline,   Depth)),
                                                                                    zScore_mean=mean((Depth-mean(Baseline))/sd(Baseline)),
                                                                                    zScore_med=median((Depth-mean(Baseline))/sd(Baseline))
    )
    )
  }
  return(list(waterErr,iceErr,snowErr))
  print(Sys.time()-starttime)
}

##---------------------------------------------------------------
##                        Toolik errors                        --
##---------------------------------------------------------------

runPre<-'T'
scenToolik<-c(paste('S',runPre,'P',80,'w',sep='_'),paste('S',runPre,'P',120,'w',sep='_'),paste('S',runPre,'P',50,sep='_'),paste('S',runPre,'P',80,sep='_'),paste('S',runPre,'P',90,sep='_'),paste('S',runPre,'P',110,sep='_'),paste('S',runPre,'P',120,sep='_'),paste('S',runPre,'P',200,sep='_'),paste('S',runPre,'SW',80,sep='_'),paste('S',runPre,'SW',90,sep='_'),paste('S',runPre,'SW',110,sep='_'),paste('S',runPre,'SW',120,sep='_'),paste('S',runPre,'T_K',-2,sep='_'),paste('S',runPre,'T_K',-1,sep='_'),paste('S',runPre,'T_K',1,sep='_'),paste('S',runPre,'T_K',2,sep='_')
)

TooErr<-calcErrorScenarioSeason(obsDataFile='/Volumes/GoogleDrive/My Drive/Lake/Toolik_UAF/lake temp/lakeTemp.csv',
                          modDataDir ='/Volumes/GoogleDrive/My Drive/Lake/Toolik_UAF/model output',
                          modWxFile  ='/Volumes/GoogleDrive/My Drive/Lake/Toolik_UAF/met data/TooHrDate.dat',
                          scenDataDir = '/Users/jac/Documents/GitHub/LAKE/results',
                          Scenarios=scenToolik,
                          Depths=c(1,2,3,5,10,19),
                          Season='Season_Ice')
##---------------------------------------------------------------
##                        Atqasuk errors                      --
##---------------------------------------------------------------

runPre<-'A'
scenAtq<-c(paste('S',runPre,'P',80,'w',sep='_'),paste('S',runPre,'P',120,'w',sep='_'),paste('S',runPre,'P',50,sep='_'),paste('S',runPre,'P',80,sep='_'),paste('S',runPre,'P',90,sep='_'),paste('S',runPre,'P',110,sep='_'),paste('S',runPre,'P',120,sep='_'),paste('S',runPre,'P',200,sep='_'),paste('S',runPre,'SW',80,sep='_'),paste('S',runPre,'SW',90,sep='_'),paste('S',runPre,'SW',110,sep='_'),paste('S',runPre,'SW',120,sep='_'),paste('S',runPre,'T_K',-2,sep='_'),paste('S',runPre,'T_K',-1,sep='_'),paste('S',runPre,'T_K',1,sep='_'),paste('S',runPre,'T_K',2,sep='_')
)

AtqErr<-calcErrorScenarioSeason(obsDataFile='/Volumes/GoogleDrive/My Drive/Lake/Atgasuk_UAF/lake temp/lakeTemp.csv',
                          modDataDir ='/Volumes/GoogleDrive/My Drive/Lake/Atgasuk_UAF/model output',
                          modWxFile  ='/Volumes/GoogleDrive/My Drive/Lake/Atgasuk_UAF/model_inputs/AtgasukDate.dat',
                          scenDataDir = '/Users/jac/Documents/GitHub/LAKE/results',
                          Scenarios=scenAtq,
                          Depths=c(0.3,2.5),
                          Season='Season_Ice')


##---------------------------------------------------------------
##                        Foxden errors                        --
##---------------------------------------------------------------
runPre<-'F'
scenFox<-c(paste('S',runPre,'P',80,'w',sep='_'),paste('S',runPre,'P',120,'w',sep='_'),paste('S',runPre,'P',50,sep='_'),paste('S',runPre,'P',80,sep='_'),paste('S',runPre,'P',90,sep='_'),paste('S',runPre,'P',110,sep='_'),paste('S',runPre,'P',120,sep='_'),paste('S',runPre,'P',200,sep='_'),paste('S',runPre,'SW',80,sep='_'),paste('S',runPre,'SW',90,sep='_'),paste('S',runPre,'SW',110,sep='_'),paste('S',runPre,'SW',120,sep='_'),paste('S',runPre,'T_K',-2,sep='_'),paste('S',runPre,'T_K',-1,sep='_'),paste('S',runPre,'T_K',1,sep='_'),paste('S',runPre,'T_K',2,sep='_')
)

FoxErr<-calcErrorScenarioSeason(obsDataFile='/Volumes/GoogleDrive/My Drive/Lake/FoxDen_UAF/lake temp/lakeTemp.csv',
                          modDataDir ='/Volumes/GoogleDrive/My Drive/Lake/FoxDen_UAF/model output',
                          modWxFile  ='/Volumes/GoogleDrive/My Drive/Lake/FoxDen_UAF/met data/FoxDenDayDate.dat',
                          scenDataDir = '/Users/jac/Documents/GitHub/LAKE/results',
                          Scenarios=scenFox,
                          Depths=1.5,
                          Season='Season_Ice')
##---------------------------------------------------------------
##                       save errors                           --
##---------------------------------------------------------------
saveRDS(list(TooErr,AtqErr,FoxErr),'error_season.RDS')

##----------------------------------------------------------------
##                    combine, format, label                    --
##----------------------------------------------------------------

waterError<-rbind(data.frame(Lake='Toolik', TooErr[[1]] %>% filter(Depth%in%c(1,3,5,10,19))),
                  data.frame(Lake='Atqasuk',AtqErr[[1]] %>% filter(Depth%in%c(2.5))),
                  data.frame(Lake='FoxDen', FoxErr[[1]] %>% filter(Depth%in%c(1.5))))
waterError_L<-waterError %>% pivot_longer(cols=c("MAE","Bias","MSE","PB","RAE","RMSE","cumAE","zScore_mean","zScore_med"),names_to = 'Metric',values_to='waterError')
waterError_L2 <- waterError_L %>% mutate(Scenario2=factor(gsub('S_[[:upper:]]_','',as.character(Scenario),perl=TRUE))) %>%
  select(-Scenario) %>% rename(Scenario=Scenario2)
waterError_L2$Scenario<-factor(waterError_L2$Scenario,
                               labels = c("P 110%","P 120%","P 120% w","P 200%","P 50%","P 80%","P 80% w","P 90%","SW 110%","SW 120%","SW 80%","SW 90%","TA -1C","TA -2C","TA +1C","TA +2C"))

##----------------------------------------------------------------
##              calculate mean error across depths              --
##----------------------------------------------------------------
waterError_Lmean<-waterError_L2 %>% group_by(Lake, Scenario, Metric,Season_Ice) %>% summarise(waterErrorMean=mean(waterError),waterErrorAbsMean=mean(abs(waterError)))

##---------------------------------------------------------------
##                   plot and save                             --
##---------------------------------------------------------------

waterError_Lmean %>% mutate(Scenario=factor(Scenario,levels = c("P 50%","P 80%","P 80% w","P 90%","P 110%","P 120%","P 120% w","P 200%","SW 80%","SW 90%","SW 110%","SW 120%","TA -1C","TA -2C","TA +1C","TA +2C"),
                                            labels = c("P 50%","P 80%","P 80% w","P 90%","P 110%","P 120%","P 120% w","P 200%","SW 80%","SW 90%","SW 110%","SW 120%","TA -1\u00B0C","TA -2\u00B0C","TA +1\u00B0C","TA +2\u00B0C"))) %>%
  filter(Metric=='zScore_mean') %>%
  ggplot(aes(x=waterErrorMean,y=Scenario,group=Lake))+
  geom_bar(aes(fill=Season_Ice,group=Season_Ice),stat='identity',position='dodge') + 
  scale_fill_manual(values=c('#4575b4','#d73027'),name='Season')+
  xlab('Z Score')+
  theme_pubr()+
  theme(text = element_text(family='sans',size=10),
        axis.text = element_text(family='sans',size=10,margin=0),
        axis.title.y = element_text(family='sans',size=10,margin = unit(c(0,0,0,1),'mm')),
        axis.text.y = element_text(family='sans',size=10,hjust=0))+
  theme(panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "lines"),
        panel.spacing = unit(0, "null"),
        # axis.ticks.length = unit(0, "null"),
        legend.margin = margin(0,0,-0.9,0, "lines"),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.position = 'none',
        strip.background=element_blank(),
        strip.text = element_text(family='sans',size=10)
  )+
  theme(legend.position = 'top',legend.title = element_text(family='sans',size=10),legend.text = element_text(family='sans',size=10))+
  facet_wrap(~Lake,ncol=3)

ggsave(filename = paste0('/Volumes/GoogleDrive/My Drive/beavers/LAKE/Figs/Sensitivity_Bar_Zseason.png'),
       width=8,height = 6)

ggsave(filename = paste0('/Volumes/GoogleDrive/My Drive/beavers/LAKE/Figs/Sensitivity_Bar_Zseason2.png'),
       width=3.25,height = 6.5)

