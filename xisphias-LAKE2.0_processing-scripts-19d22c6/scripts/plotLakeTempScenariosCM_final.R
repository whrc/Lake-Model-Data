
##################################################################
##                Calculate errors for scenarios                ##
##          Combine errors for color matrix plot/plots          ##
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
##               Function for calculating errors               --
##---------------------------------------------------------------
calcErrorScenario<-function(obsDataFile, modDataDir, modWxFile, Scenarios=c(''), scenDataDir,Depths=c(1,2)){
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
  modWX$DateTime<-as.POSIXct(modWX$DateTime,format="%Y-%m-%d %H:%M:%S",tz="US/Alaska")
  
  #prep model data
  names(modWater)[1:5]<-c('Year','Month','Day','Hour','integrationTime')
  names(modWater)[seq(6,dim(modWater)[2],by=2)]<-paste0('Depth',seq(1,(dim(modWater)[2]-5)/2,by=1))
  names(modWater)[seq(7,dim(modWater)[2],by=2)]<-paste0('Temp',seq(1,(dim(modWater)[2]-5)/2,by=1))
  modWater$DateTime<-with(modWater,as.POSIXct(strptime(paste0(Year,'-',Month,'-',Day,'H',Hour),format = '%Y-%m-%dH%H'),tz='US/Alaska'))
  print('model dates:');print(min(na.omit(modWater$DateTime)));print(max(na.omit(modWater$DateTime)))
  m1gD<-gather(modWater %>% select(-(starts_with('Temp'))),colName,Depth,-Year,-Month,-Day,-Hour,-integrationTime,-DateTime)
  m1gT<-gather(modWater %>% select(-(starts_with('Depth'))),colName,Temp,-Year,-Month,-Day,-Hour,-integrationTime,-DateTime)
  calcm<-cbind(m1gD %>% select(-colName) %>% mutate(Depth=factor(Depth)),Temp=m1gT %>% pull(Temp))
  
  #prep soil
  names(modSoil)[1:5]<-c('Year','Month','Day','Hour','integrationTime')
  names(modSoil)[seq(6,dim(modSoil)[2],by=2)]<-paste0('Depth',seq(1,(dim(modSoil)[2]-5)/2,by=1))
  names(modSoil)[seq(7,dim(modSoil)[2],by=2)]<-paste0('Temp',seq(1,(dim(modSoil)[2]-5)/2,by=1))
  modSoil$DateTime<-with(modSoil,as.POSIXct(strptime(paste0(Year,'-',Month,'-',Day,'H',Hour),format = '%Y-%m-%dH%H'),tz='US/Alaska'))
  modSoilgD<-gather(modSoil %>% select(-(starts_with('Temp'))),colName,Depth,-Year,-Month,-Day,-Hour,-integrationTime,-DateTime)
  modSoilgT<-gather(modSoil %>% select(-(starts_with('Depth'))),colName,Temp,-Year,-Month,-Day,-Hour,-integrationTime,-DateTime)
  Soilm<-cbind(modSoilgD %>% select(-colName) %>% mutate(Depth=factor(Depth)),Temp=modSoilgT %>% pull(Temp))
  
  #prep layers
  names(modLayers)<-c("year","month","day","hour","intTimeHr","waterlayer_m","Wmixedlayerthickness_m","Emixedlayerthickness_m","Smixedlayerthickness_m","Nmixedlayerthickness_m","Wlowerlayerthickness_m","Elowerlayerthickness_m","Slowerlayerthickness_m","Nlowerlayerthickness_m","icelayer_m","snowlayer_m","bottomice_m","reservoirvolume_m3","volumedeficit_accumulated_m3")
  modLayers$DateTime<-as.POSIXct(strptime(with(modLayers,paste0(year,'-',month,'-',day,' ',as.integer(hour),':00:00')),format='%Y-%m-%d %H:%M:%S'),tz='US/Alaska')
  modLayers$Date<-as.POSIXct(format(modLayers$DateTime,format='%Y-%m-%d'),tz='US/Alaska')
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
      # print('model dates:');print(min(na.omit(modWaterScen$DateTime)));print(max(na.omit(modWaterScen$DateTime)))
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
      scenSoil<-rbind(scenSoil,data.frame(Scenario=Scenarios[a],SoilmScen))
      
      #prep layers
      names(modLayersScen)<-c("Year","Month","Day","Hour","intTimeHr","waterlayer_m","Wmixedlayerthickness_m","Emixedlayerthickness_m","Smixedlayerthickness_m","Nmixedlayerthickness_m","Wlowerlayerthickness_m","Elowerlayerthickness_m","Slowerlayerthickness_m","Nlowerlayerthickness_m","icelayer_m","snowlayer_m","bottomice_m","reservoirvolume_m3","volumedeficit_accumulated_m3")
      modLayersScen$DateTime<-as.POSIXct(strptime(with(modLayersScen,paste0(Year,'-',Month,'-',Day,' ',as.integer(Hour),':00:00')),format='%Y-%m-%d %H:%M:%S'),tz='US/Alaska')
      modLayersScen$Date<-as.POSIXct(format(modLayersScen$DateTime,format='%Y-%m-%d'),tz='US/Alaska')
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
        newdat$Temp<-as.numeric(newdat$Temp)
        newdat$Depth<-factor(newdat$Depth)
        # need to drop NA and remove duplicate records
        newdat_W <- newdat %>% drop_na() %>% unique() %>% pivot_wider(id_col=c(DateTime,Depth),names_from=Scenario,values_from=c(Temp),names_repair='unique')
        newDat_L<-newdat_W %>% pivot_longer(cols=starts_with('S_'),names_to = 'Scenario',values_to='Temp') %>% mutate(Scenario=factor(Scenario))
        (waterErr<-newDat_L %>% drop_na() %>% group_by(Scenario,Depth) %>% summarise(MAE=mae(Baseline,Temp), #mean abs error
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
      newDat_L<-newdat_W %>% pivot_longer(cols=starts_with('S_'),names_to = 'Scenario',values_to='Depth') %>% mutate(Scenario=factor(Scenario))
      (iceErr<-newDat_L %>% filter(Baseline>0) %>% group_by(Scenario)  %>% summarise(MAE=mae(Baseline,        Depth), #mean abs error
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
      newDat_L<-newdat_W %>% pivot_longer(cols=starts_with('S_'),names_to = 'Scenario',values_to='Depth') %>% mutate(Scenario=factor(Scenario))
      (snowErr<-newDat_L %>% filter(Baseline>0) %>% group_by(Scenario)  %>% summarise(MAE=mae(Baseline,        Depth), #mean abs error
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

TooErr<-calcErrorScenario(obsDataFile='/Volumes/GoogleDrive/My Drive/Lake/Toolik_UAF/lake temp/lakeTemp.csv',
                        modDataDir ='/Volumes/GoogleDrive/My Drive/Lake/Toolik_UAF/model output',
                        modWxFile  ='/Volumes/GoogleDrive/My Drive/Lake/Toolik_UAF/met data/TooHrDate.dat',
                        scenDataDir = '/Users/jac/Documents/GitHub/LAKE/results',
                        Scenarios=scenToolik,
                        Depths=c(1,2,3,5,10,19))

##---------------------------------------------------------------
##                        Atqasuk errors                      --
##---------------------------------------------------------------


runPre<-'A'
scenAtq<-c(paste('S',runPre,'P',80,'w',sep='_'),paste('S',runPre,'P',120,'w',sep='_'),paste('S',runPre,'P',50,sep='_'),paste('S',runPre,'P',80,sep='_'),paste('S',runPre,'P',90,sep='_'),paste('S',runPre,'P',110,sep='_'),paste('S',runPre,'P',120,sep='_'),paste('S',runPre,'P',200,sep='_'),paste('S',runPre,'SW',80,sep='_'),paste('S',runPre,'SW',90,sep='_'),paste('S',runPre,'SW',110,sep='_'),paste('S',runPre,'SW',120,sep='_'),paste('S',runPre,'T_K',-2,sep='_'),paste('S',runPre,'T_K',-1,sep='_'),paste('S',runPre,'T_K',1,sep='_'),paste('S',runPre,'T_K',2,sep='_')
)

AtqErr<-calcErrorScenario(obsDataFile='/Volumes/GoogleDrive/My Drive/Lake/Atgasuk_UAF/lake temp/lakeTemp.csv',
                        modDataDir ='/Volumes/GoogleDrive/My Drive/Lake/Atgasuk_UAF/model output',
                        modWxFile  ='/Volumes/GoogleDrive/My Drive/Lake/Atgasuk_UAF/model_inputs/AtgasukDate.dat',
                        scenDataDir = '/Users/jac/Documents/GitHub/LAKE/results',
                        Scenarios=scenAtq,
                        Depths=c(0.3,2.5))


##---------------------------------------------------------------
##                        Foxden errors                        --
##---------------------------------------------------------------


runPre<-'F'
scenFox<-c(paste('S',runPre,'P',80,'w',sep='_'),paste('S',runPre,'P',120,'w',sep='_'),paste('S',runPre,'P',50,sep='_'),paste('S',runPre,'P',80,sep='_'),paste('S',runPre,'P',90,sep='_'),paste('S',runPre,'P',110,sep='_'),paste('S',runPre,'P',120,sep='_'),paste('S',runPre,'P',200,sep='_'),paste('S',runPre,'SW',80,sep='_'),paste('S',runPre,'SW',90,sep='_'),paste('S',runPre,'SW',110,sep='_'),paste('S',runPre,'SW',120,sep='_'),paste('S',runPre,'T_K',-2,sep='_'),paste('S',runPre,'T_K',-1,sep='_'),paste('S',runPre,'T_K',1,sep='_'),paste('S',runPre,'T_K',2,sep='_')
)
FoxErr<-calcErrorScenario(obsDataFile='/Volumes/GoogleDrive/My Drive/Lake/FoxDen_UAF/lake temp/lakeTemp.csv',
                        modDataDir ='/Volumes/GoogleDrive/My Drive/Lake/FoxDen_UAF/model output',
                        modWxFile  ='/Volumes/GoogleDrive/My Drive/Lake/FoxDen_UAF/met data/FoxDenDayDate.dat',
                        scenDataDir = '/Users/jac/Documents/GitHub/LAKE/results',
                        Scenarios=scenFox,
                        Depths=1.5)


##---------------------------------------------------------------
##                       save errors                           --
##---------------------------------------------------------------
saveRDS(list(TooErr,AtqErr,FoxErr),'error.RDS')


##---------------------------------------------------------------
##                         bind errors                         --
##---------------------------------------------------------------
  waterError<-rbind(data.frame(Lake='Toolik', TooErr[[1]] %>% filter(Depth%in%c(1,3,5,10,19))),
                  data.frame(Lake='Atqasuk',AtqErr[[1]] %>% filter(Depth%in%c(2.5))),
                  data.frame(Lake='FoxDen', FoxErr[[1]] %>% filter(Depth%in%c(1.5))))
  waterError_L<-waterError %>% pivot_longer(cols=c("MAE","Bias","MSE","PB","RAE","RMSE","cumAE","zScore_mean","zScore_med"),names_to = 'Metric',values_to='waterError')
  # spread watererror out by selecting depths for each lake and renaming error by depth
    waterError_LW<-waterError_L %>% pivot_wider(id_cols=c(Lake,Scenario,Metric),names_from = Depth,values_from=waterError)
  #ice
  iceError<-rbind(data.frame(Lake='Toolik',   TooErr[[2]]),
                    data.frame(Lake='Atqasuk',AtqErr[[2]]),
                    data.frame(Lake='FoxDen', FoxErr[[2]]))
  iceError_L<-iceError %>% pivot_longer(cols=c("MAE","Bias","MSE","PB","RAE","RMSE","cumAE","zScore_mean","zScore_med"),names_to = 'Metric',values_to='Ice')
  #snow  
  snowError<-rbind(data.frame(Lake='Toolik',   TooErr[[3]]),
                  data.frame(Lake='Atqasuk',  AtqErr[[3]]),
                  data.frame(Lake='FoxDen',   FoxErr[[3]]))
  snowError_L<-snowError %>% pivot_longer(cols=c("MAE","Bias","MSE","PB","RAE","RMSE","cumAE","zScore_mean","zScore_med"),names_to = 'Metric',values_to='Snow')
  
##----------------------------------------------------------------
##                         merge errors                         --
##----------------------------------------------------------------
  
  
  allError<-merge(merge(waterError_LW,iceError_L,by=c('Lake','Scenario','Metric'),all = TRUE),
                snowError_L,by=c('Lake','Scenario','Metric'),all = TRUE)
##---------------------------------------------------------------
##                new scenario + variable names                --
##---------------------------------------------------------------
  
  allError <- allError %>% mutate(Scenario2=factor(gsub('S_[[:upper:]]_','',as.character(Scenario),perl=TRUE))) %>%
                           select(-Scenario) %>% rename(Scenario=Scenario2)
  allError <- allError %>% mutate(Scenario=factor(Scenario,
                                 labels = c("P 110%","P 120%","P 120% w","P 200%","P 50%","P 80%","P 80% w","P 90%","SW 110%","SW 120%","SW 80%","SW 90%","TA -1C","TA -2C","TA +1C","TA +2C")))
  allError <- allError %>% mutate(Scenario=factor(Scenario,levels = c("P 50%","P 80%","P 80% w","P 90%","P 110%","P 120%","P 120% w","P 200%","SW 80%","SW 90%","SW 110%","SW 120%","TA -1C","TA -2C","TA +1C","TA +2C"),
                                              labels = c("P 50%","P 80%","P 80% w","P 90%","P 110%","P 120%","P 120% w","P 200%","SW 80%","SW 90%","SW 110%","SW 120%","TA -1\u00B0C","TA -2\u00B0C","TA +1\u00B0C","TA +2\u00B0C")))
  allError_L<-allError %>% pivot_longer(cols = c('19','10','5','3','2.5','1.5','1','Ice','Snow'),names_to='Type',values_to='Error') %>% 
                           mutate(Lake=factor(Lake),Metric=factor(Metric),Type=factor(Type))
  allError_L$Type<-factor(allError_L$Type,levels=c('1','1.5','2.5','3','5','10','19','Ice','Snow'),labels = c('1m','1.5m','2.5m','3m','5m','10m','19m','Ice','Snow'))
  saveRDS(allError_L,file='allError_L.RDS')

  
##---------------------------------------------------------------
##                        create plots                         --
##---------------------------------------------------------------
  
limits<-c(allError_L %>% filter(Metric=='zScore_mean') %>% pull(Error) %>% na.omit() %>% min(),
         allError_L %>% filter(Metric=='zScore_mean') %>% pull(Error) %>% na.omit() %>% max())
p1<-ggplot(allError_L %>% filter(Metric=='zScore_mean',Lake=='Atqasuk') %>% drop_na(),aes(y=Type,x=Scenario,group=Lake))+
  geom_tile(aes(fill = Error)) + 
  geom_text(aes(label = round(Error, 2)), size = 3) + 
  coord_fixed() + 
  scale_fill_fermenter(type = 'div',palette = 'RdYlBu',limits=limits,name='Z Score')+
  theme_bw()+
  theme(panel.border = element_blank())+
  theme(axis.text = element_text(family = 'sans',size=10))+
  theme(axis.ticks.x=element_blank())+
  theme(panel.grid =element_blank())+
  ylab('Atqasuk')

p2<-ggplot(allError_L %>% filter(Metric=='zScore_mean',Lake=='FoxDen') %>% drop_na(),aes(y=Type,x=Scenario,group=Lake))+
  geom_tile(aes(fill = Error)) + 
  geom_text(aes(label = round(Error, 2)), size = 3) + 
  coord_fixed() + 
  scale_fill_fermenter(type = 'div',palette = 'RdYlBu',limits=limits)+
  theme_bw()+
  theme(panel.border = element_blank())+
  theme(axis.text = element_text(family = 'sans',size=10))+
  theme(axis.ticks.x=element_blank())+
  theme(panel.grid =element_blank())+
  ylab('FoxDen')

p3<-ggplot(allError_L %>% filter(Metric=='zScore_mean',Lake=='Toolik') %>% drop_na(),aes(y=Type,x=Scenario,group=Lake))+
  geom_tile(aes(fill = Error)) + 
  geom_text(aes(label = round(Error, 2)), size = 3) + 
  coord_fixed() + 
  scale_fill_fermenter(type = 'div',palette = 'RdYlBu',limits=limits)+
  theme_bw()+
  theme(panel.border = element_blank())+
  theme(axis.text = element_text(family = 'sans',size=10))+
  theme(axis.text.x = element_text(angle = 90,hjust=0))+
  theme(panel.grid =element_blank())+
  ylab('Toolik')

##---------------------------------------------------------------
##                   arrange plots and save                    --
##---------------------------------------------------------------
  
  g1<-ggplotGrob(p1+theme(legend.position = 'top' ,legend.margin = margin(0,0,0,0),legend.box.margin = margin(0,0,-11,0),
                                                   axis.title.x = element_blank(),axis.text.x = element_blank())+
                   theme(plot.margin = margin(1,1,0,2)))
  g2<-ggplotGrob(p2+theme(legend.position = 'none',axis.title.x = element_blank(),axis.text.x = element_blank())+
                   theme(plot.margin = margin(-6,1,0,1)))
  g3<-ggplotGrob(p3+theme(legend.position = 'none')+
                   theme(plot.margin = margin(-4,1,1,1)))
  g<-rbind(g1,g2,g3)
  g$widths <- unit.pmax(g1$widths, g2$widths,g3$widths)
  gg<-arrangeGrob(g)
  grid.arrange(gg)
  ggsave(gg,filename = paste0('/Volumes/GoogleDrive/My Drive/beavers/LAKE/Figs/Sensitivity_CM_Z.png'),
         width=6.5)
