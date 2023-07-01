#################################################################
##                    Plot LAKE model temp data                ##
##                    compared to obs data and met data        ##
##                      Jason Clark                            ##
##                  jaclark2@alaska.edu                        ##
##                          2021-10-22                         ##
#################################################################


##----------------------------------------------------------------
##                 Function for plotting                        --
##----------------------------------------------------------------

plotLakeDepth<-function(obsDataFile, modDataDir, modWxFile, Depths=c(1,2,3),dateBreak='3 months',dateBreak_minor='1 month',plotPrefix='Test',outDir='/Test/',P_mult=60*60,sec_axis_coeff=10){
  #P_mult 60*60 for hourly, 60*60*24 for daily
  #sec_axis_coeff for scaling Precip axis and precip data, 50 for daily, 10 for monthly
  library(tidyverse)
  library(gtable)
  library(grid)
  library(gridExtra)
  library(ggh4x)
  library(ggtext)
  dir.create(outDir)
  setwd(outDir)
  measured<-read.csv(obsDataFile)
  modWater =read.table(paste0(modDataDir,'/time_series/water_temp  1  1.dat'),skip=6)
  modSoil  =read.table(paste0(modDataDir,'/time_series/soil_temp  1  1.dat'),skip=6)
  modLayers=read.table(paste0(modDataDir,'/time_series/layers  1  1.dat'),skip=19)
  modWX    =read.csv(modWxFile,header=FALSE)
  
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
    modWXday<-modWX %>% drop_na() %>% mutate(Date=as.POSIXct(format(DateTime,format='%Y-%m-%d'),tz='US/Alaska')) %>% group_by(Date) %>% summarise(P=sum(P*P_mult,na.rm = TRUE))
    #converted P from m/s to m/hr then summed to m/day
    modWXmonth<-modWX %>% drop_na() %>% mutate(Date=as.POSIXct(paste0(format(DateTime,format='%Y-%m'),'-28'),tz='US/Alaska')) %>% group_by(Date) %>% summarise(P=sum(P*P_mult,na.rm = TRUE))
    #converted P from m/s to m/hr then summed to m/month
    
    #prep model data
    names(modWater)[1:5]<-c('Year','Month','Day','Hour','integrationTime')
    names(modWater)[seq(6,dim(modWater)[2],by=2)]<-paste0('Depth',seq(1,(dim(modWater)[2]-5)/2,by=1))
    names(modWater)[seq(7,dim(modWater)[2],by=2)]<-paste0('Temp',seq(1,(dim(modWater)[2]-5)/2,by=1))
    modWater$DateTime<-with(modWater,as.POSIXct(strptime(paste0(Year,'-',Month,'-',Day,'H',Hour),format = '%Y-%m-%dH%H'),tz='US/Alaska'))
    print('model dates:');print(min(na.omit(modWater$DateTime)));print(max(na.omit(modWater$DateTime)))
    m1gD<-gather(modWater %>% select(-(starts_with('Temp'))),colName,Depth,-Year,-Month,-Day,-Hour,-integrationTime,-DateTime)
    m1gT<-gather(modWater %>% select(-(starts_with('Depth'))),colName,Temp,-Year,-Month,-Day,-Hour,-integrationTime,-DateTime)
    calcm<-cbind(m1gD %>% select(-colName) %>% mutate(Depth=factor(Depth)),Temp=m1gT %>% pull(Temp))
    print(paste('Model Water Depths:',paste(levels(calcm$Depth),collapse = ', ')))
    
    #prep soil
    names(modSoil)[1:5]<-c('Year','Month','Day','Hour','integrationTime')
    names(modSoil)[seq(6,dim(modSoil)[2],by=2)]<-paste0('Depth',seq(1,(dim(modSoil)[2]-5)/2,by=1))
    names(modSoil)[seq(7,dim(modSoil)[2],by=2)]<-paste0('Temp',seq(1,(dim(modSoil)[2]-5)/2,by=1))
    modSoil$DateTime<-with(modSoil,as.POSIXct(strptime(paste0(Year,'-',Month,'-',Day,'H',Hour),format = '%Y-%m-%dH%H'),tz='US/Alaska'))
    modSoilgD<-gather(modSoil %>% select(-(starts_with('Temp'))),colName,Depth,-Year,-Month,-Day,-Hour,-integrationTime,-DateTime)
    modSoilgT<-gather(modSoil %>% select(-(starts_with('Depth'))),colName,Temp,-Year,-Month,-Day,-Hour,-integrationTime,-DateTime)
    Soilm<-cbind(modSoilgD %>% select(-colName) %>% mutate(Depth=factor(Depth)),Temp=modSoilgT %>% pull(Temp))
    print(paste('Soil Depths:',paste(levels(Soilm$Depth),collapse = ', ')))
  
    #prep layers
    names(modLayers)<-c("year","month","day","hour","intTimeHr","waterlayer_m","Wmixedlayerthickness_m","Emixedlayerthickness_m","Smixedlayerthickness_m","Nmixedlayerthickness_m","Wlowerlayerthickness_m","Elowerlayerthickness_m","Slowerlayerthickness_m","Nlowerlayerthickness_m","icelayer_m","snowlayer_m","bottomice_m","reservoirvolume_m3","volumedeficit_accumulated_m3")
    modLayers$DateTime<-as.POSIXct(strptime(with(modLayers,paste0(year,'-',month,'-',day,' ',as.integer(hour),':00:00')),format='%Y-%m-%d %H:%M:%S'),tz='US/Alaska')
    modLayers$Date<-as.POSIXct(format(modLayers$DateTime,format='%Y-%m-%d'),tz='US/Alaska')
    
  #loop through depths and plot
  tempPlots<-vector('list')
  nDepth<-length(Depths)
  a=1
  for(a in 1:nDepth){
    thisDepth=Depths[a]
    thisDate3<-min(calcm$DateTime,na.rm = TRUE)
    thisDate4<-max(calcm$DateTime,na.rm = TRUE)
    
    dfHrly<-merge(data.frame(measured %>% filter(Depth==(thisDepth)) %>% filter(DateTime<=thisDate4,DateTime>thisDate3) %>% rename(temp_c=Temp_C)  %>% select(DateTime,observed=temp_c)),
                  data.frame( (calcm %>% filter(Depth==(thisDepth)) %>% select(DateTime,modeled=Temp) )) ,by='DateTime',all=TRUE)
    errHrly=mean(abs(dfHrly %>% mutate(err=observed-modeled) %>% pull(err) %>% na.omit()))
    dfHrlyg<-gather(dfHrly, 'type', 'temp_c', -DateTime)
    dfHrlyg$type<-factor(dfHrlyg$type)
    thisDate1<-min(dfHrlyg %>% filter(type=='observed') %>% pull(DateTime),na.rm=TRUE)
    thisDate2<-max(dfHrlyg %>% filter(type=='observed') %>% pull(DateTime),na.rm=TRUE)
    p<-ggplot(data=dfHrlyg,aes(x=DateTime,y=temp_c,color=type))+
      geom_abline(intercept = 0,slope = 0,color=4)+
      geom_line(alpha=0.9)+
      scale_color_manual(values=c('black','red'),name='')+
      xlab('Date')+
      scale_y_continuous(name='Temp. &deg;C',limits=c(0,20),
                         sec.axis = sec_axis(trans = ~.*1, name=paste0(thisDepth," m")))+
      theme_bw()+
      scale_x_datetime(date_breaks = dateBreak,date_minor_breaks = dateBreak_minor,date_labels = "%Y-%m",
                       limits = c(max(thisDate1,thisDate3),min(thisDate2,thisDate4)),
                       guide = "axis_minor" # this is added to the original code)+
      ) +
      theme(ggh4x.axis.ticks.length.minor = rel(1))+
      theme(axis.ticks.y.right = element_blank(),axis.text.y.right = element_blank())+
      theme(axis.text = element_text(size=10,family = 'sans'))+
      theme(axis.title = element_text(size=10,family = 'sans'))+
      theme(axis.title.y = element_markdown(size=10,family = 'sans'))+
      theme(legend.text = element_text(size=10,family = 'sans'))+
      theme(axis.text.x = element_text(size=10,family = 'sans'))+
      theme(panel.grid = element_blank())
    #p
    tempPlots[[a]]<-p
    ggsave(plot=p,filename=paste0(plotPrefix,'_mod-obs_',thisDepth,'m.png'),width=8,height = 3)
    print(paste0(plotPrefix,'_mod-obs_',thisDepth,'m.png'))
  }
  #plot water contour
  {
    calcm2<-calcm %>% arrange(DateTime,(-as.numeric(as.character(Depth))))
    calcG2<-with(calcm2, data.frame(x=DateTime,y=(-as.numeric(as.character(Depth))),z=Temp) )
    calcG2<-calcG2[calcG2$y!='-0.3',]
    p1<-ggplot(calcG2)+
      geom_raster(aes(x=x,y=y,z=z,fill=z),interpolate = FALSE, vjust=0)+
      geom_line(data=mlayersDay %>% mutate(icelayer_m=ifelse(icelayer_m==0,NA,icelayer_m)), aes(x=Date,y=icelayer_m*(-1)),color=1,show.legend = FALSE)+
      scale_fill_gradient2(low='blue',mid='white',high='red',aesthetics = c('fill'),limits=c(-20,20))+
      ylim(-25.1,0)+ylab('Depth, m')+xlab('')+guides(fill=guide_colorbar(title="Temperature, C"))+
      theme_bw()+
      scale_x_datetime(date_breaks = dateBreak,date_labels = "%Y-%m",limits = c(max(thisDate1,thisDate3),min(thisDate2,thisDate4)))+
      theme(axis.text.x = element_text(angle=30,hjust=1,vjust=1))
    # p1
    ggsave(plot=p1,filename=paste0(plotPrefix,'_contourWater.png'),width=8,height = 3)
    print(paste0(plotPrefix,'_contourWater.png'))
  }
  #plot soil contour
  {
    Soilm2<-Soilm %>% arrange(integrationTime,(-as.numeric(as.character(Depth))))
    SoilG<-with(Soilm2 %>% drop_na(), data.frame(x=DateTime,y=(-as.numeric(as.character(Depth))),z=Temp) )
    p2<-ggplot(SoilG, aes(x=x,y=y,z=z))+
      geom_raster(aes(fill=z),interpolate = TRUE,vjust=0)+
      geom_contour(aes(color=..level..),color=1,breaks=c(-100,0,100),show.legend = FALSE)+
      scale_fill_gradient2(low='blue',mid='white',high='red',aesthetics = c('fill'),limits=c(-5,10))+
      ylab('Depth, m')+guides(fill=guide_colorbar(title="Temperature, C"))+
      theme_bw()+
      scale_x_datetime(date_breaks = dateBreak,date_labels = "%Y-%m",limits = c(max(thisDate1,thisDate3),min(thisDate2,thisDate4)))+
      theme(axis.text.x = element_text(angle=30,hjust=1,vjust=1))
    # p2
    ggsave(plot=p2,filename=paste0(plotPrefix,'_contourSoil.png'),width=8,height = 3)
    print(paste0(plotPrefix,'_contourSoil.png'))
  }
  #plot snow, ice
  {
    layerg<-gather(modLayers %>% select(DateTime,icelayer_m,snowlayer_m),'var','value',-DateTime)
    layerg$var<-factor(layerg$var,levels=c("icelayer_m","snowlayer_m"),labels=c("Ice","Snow"))
    p3<-ggplot(layerg,aes(x=DateTime,y=value))+
      geom_bar(data=modWXmonth,aes(x=Date,y=P*sec_axis_coeff),stat='identity',width = 0.1,color=grey(0.1),fill=grey(0.1))+
      geom_line(aes(lty=var,color=var))+
      scale_linetype_manual(values=c(1,2),name='')+
      scale_color_manual(values=c(2,4),name='')+
      scale_y_continuous(
        name = "Depth, m",
        sec.axis = sec_axis( trans=~./sec_axis_coeff, name="Precip., m"),
        limits = c(0,1.5))+
      xlab('Date')+
      theme_bw()+
      scale_x_datetime(date_breaks = dateBreak,date_minor_breaks = dateBreak_minor,date_labels = "%Y-%m",
                       limits = c(max(thisDate1,thisDate3),min(thisDate2,thisDate4)),
                       guide = "axis_minor" # this is added to the original code)+
                       ) +
      theme(ggh4x.axis.ticks.length.minor = rel(1))+
      theme(axis.text = element_text(size=10,family = 'sans'))+
      theme(axis.title = element_text(size=10,family = 'sans'))+
      theme(legend.text = element_text(size=10,family = 'sans'))+
      theme(axis.text.x = element_text(size=10,family = 'sans'))+
      theme(legend.position = c(.99, 0.001),legend.justification = c('center'),legend.background = element_blank(),legend.box.background = element_blank())+
      theme(panel.grid = element_blank())
    #p3
    ggsave(plot=p3,filename=paste0(plotPrefix,'_iceSnow.png'),width=8,height = 3)
    print(paste0(plotPrefix,'_iceSnow.png'))
  }
  
  #combine temp and ice/snow/P plots
  if(length(tempPlots)==1){
    g1 <- ggplotGrob(
      tempPlots[[1]]+theme(axis.text.x = element_blank())+theme(axis.title.x=element_blank())+
        theme(legend.position = "top",legend.margin = margin(0,0,0,0),plot.margin = margin(1,-5,0,0),legend.box.margin = margin(-5,0,-11,0))+
        labs(tag='a')+
        theme(plot.tag.position = c(0,1),
              plot.tag = element_text(family = 'sans',size=14,
                                      margin=margin(6,15,0,2)))
    )
    g3 <-ggplotGrob(
      p3+theme(axis.text.x = element_text(family='sans',size=10,angle = 15,hjust=1,vjust=1))+
        theme(legend.position = "top",legend.margin = margin(0,0,0,0),plot.margin = margin(0,-5,0,11),legend.box.margin = margin(0,0,-11,0))+
        labs(tag='b')+
        theme(plot.tag.position = c(0,1),
              plot.tag = element_text(family = 'sans',size=14,
                                      margin=margin(12,15,0,2)))
    )
    
    g <- rbind(g1, g3, size = "first")
    g$widths <- unit.pmax(g1$widths, g3$widths)
    gg<-arrangeGrob(g)
    ggsave(plot=gg,filename=paste0(plotPrefix,'_tempIceSnow.png'),width=6.5,height = 2)
    print(paste0(plotPrefix,'_tempIceSnow.png'))
  }
  if(length(tempPlots)>1){
    g1 <- ggplotGrob(
                      tempPlots[[1]]+theme(axis.text.x = element_blank())+theme(axis.title.x=element_blank())+
                      theme(legend.position = "top",legend.margin = margin(0,0,0,0),plot.margin = margin(1,-5,0,0),legend.box.margin = margin(-5,0,-11,0))+
                      labs(tag='a')+
                      theme(plot.tag.position = c(0,1),
                      plot.tag = element_text(family = 'sans',size=14,
                                      margin=margin(6,15,0,2)))
      )
    g2 <- gtable_add_rows(
          ggplotGrob(
                      tempPlots[[2]]+theme(axis.text.x = element_blank())+theme(axis.title.x=element_blank())+
                      theme(legend.position = 'none',                              plot.margin = margin(1,-5,0,0))+
                      labs(tag='b')+
                      theme(plot.tag.position = c(0,1),
                      plot.tag = element_text(family = 'sans',size=14,
                                      margin=margin(2,15,1,2)))
                    ),
                    unit(c(0.01,0.01),"lines"))
    g3 <- ggplotGrob(
                    p3+theme(axis.text.x = element_text(family='sans',size=10,angle = 15,hjust=1,vjust=1))+
                      theme(legend.position = "top",legend.margin = margin(0,0,0,0),plot.margin = margin(0,-5,0,11),legend.box.margin = margin(0,0,-11,0))+
                    labs(tag='c')+
                    theme(plot.tag.position = c(0,1),
                    plot.tag = element_text(family = 'sans',size=14,
                                      margin=margin(12,15,0,2)))
      )
    
    g <- rbind(g1,g2, g3, size = "first")
    g$widths <- unit.pmax(g1$widths,g2$widths, g3$widths)

    gg<-arrangeGrob(g)

    ggsave(plot=gg,filename=paste0(plotPrefix,'_tempIceSnow.png'),width=6.5,height = 3)
    print(paste0(plotPrefix,'_tempIceSnow.png'))
  }
  
  #extra depths if present
  {
  if(length(tempPlots)==3){
    g1 <- ggplotGrob(
                      tempPlots[[3]]+
                      theme(axis.text.x = element_text(family='sans',size=10,angle = 15,hjust=1,vjust=1))+
                      theme(legend.position = "top",legend.margin = margin(0,0,0,0),plot.margin = margin(1,0,0,1),legend.box.margin = margin(-5,0,-11,0))
                    )
      
    g <- rbind(g1, size = "first")
      g$widths <- unit.pmax(g1$widths)

    }
    if(length(tempPlots)==4){
      g1 <- ggplotGrob(
        tempPlots[[3]]+theme(axis.text.x = element_blank())+theme(axis.title.x=element_blank())+
          theme(legend.position = "top",legend.margin = margin(0,0,0,0),plot.margin = margin(1,-5,0,0),legend.box.margin = margin(-5,0,-11,0))+
          labs(tag='a')+
          theme(plot.tag.position = c(0,1),
                plot.tag = element_text(family = 'sans',size=14,
                                        margin=margin(6,15,0,2)))
      )
      g2 <- gtable_add_rows(
              ggplotGrob(
                tempPlots[[4]]+theme(axis.text.x = element_text(family='sans',size=10,angle = 15,hjust=1,vjust=1))+
                  theme(legend.position = 'none',                        plot.margin = margin(1,-5,0,11))+
                  labs(tag='b')+
                  theme(plot.tag.position = c(0,1),
                        plot.tag = element_text(family = 'sans',size=14,
                                                margin=margin(2,15,1,2)))
              ),
              unit(c(0.01,0.01),"lines"))
      g <- rbind(g1,g2, size = "first")
      g$widths <- unit.pmax(g1$widths,g2$widths)
    }
    if(length(tempPlots)==5)
    {
      g1 <- ggplotGrob(
              tempPlots[[3]]+theme(axis.text.x = element_blank())+theme(axis.title.x=element_blank())+
              theme(legend.position = "top",legend.margin = margin(0,0,0,0),plot.margin = margin(1,-5,0,0),legend.box.margin = margin(-5,0,-11,0))+
              labs(tag='a')+
              theme(plot.tag.position = c(0,1),
                    plot.tag = element_text(family = 'sans',size=14,
                                            margin=margin(6,15,0,2)))
            )
      g2 <- gtable_add_rows(
              ggplotGrob(
              tempPlots[[4]]+theme(axis.text.x = element_blank())+theme(axis.title.x=element_blank())+
              theme(legend.position = 'none',                        plot.margin = margin(1,-5,0,0))+
              labs(tag='b')+
              theme(plot.tag.position = c(0,1),
                    plot.tag = element_text(family = 'sans',size=14,
                                            margin=margin(2,15,1,2)))
            ),
            unit(c(0.01,0.01),"lines"))
      g3 <- gtable_add_rows(
                ggplotGrob(
                tempPlots[[5]]+theme(axis.text.x = element_text(family='sans',size=10,angle = 15,hjust=1,vjust=1))+
                theme(legend.position = 'none',                        plot.margin = margin(1,-5,0,11))+
                labs(tag='c')+
                theme(plot.tag.position = c(0,1),
                      plot.tag = element_text(family = 'sans',size=14,
                                              margin=margin(2,15,1,2)))
            ),
            unit(c(0.01,0.01),"lines"))
      
      g <- rbind(g1,g2, g3, size = "first")
      g$widths <- unit.pmax(g1$widths,g2$widths, g3$widths)
    }
    if(length(tempPlots)>2){
      gg<-arrangeGrob(g)
      ggsave(plot=gg,filename=paste0(plotPrefix,'_tempExtra.png'),width=6.5,height = 3)
      print(paste0(plotPrefix,'_tempExtra.png'))
    }
  }
  
  print(paste0(plotPrefix,' plots done: ',outDir))
}

##----------------------------------------------------------------
##                            toolik                            --
##----------------------------------------------------------------
plotLakeDepth(obsDataFile='/Volumes/GoogleDrive/My Drive/Lake/Toolik_UAF/lake temp/lakeTemp.csv',
              modDataDir ='/Volumes/GoogleDrive/My Drive/Lake/Toolik_UAF/model output',
              modWxFile  ='/Volumes/GoogleDrive/My Drive/Lake/Toolik_UAF/met data/TooHrDate.dat',
              outDir='/Volumes/GoogleDrive/My Drive/beavers/LAKE/Figs/',
              plotPrefix='Too',
              Depths=c(1,3,5,10,19),
              dateBreak='4 months',
              dateBreak_minor='1 month',
              P_mult=60*60,
              sec_axis_coeff=10)
##----------------------------------------------------------------
##                            atqasuk                           --
##----------------------------------------------------------------
plotLakeDepth(obsDataFile='/Volumes/GoogleDrive/My Drive/Lake/Atgasuk_UAF/lake temp/lakeTemp.csv',
              modDataDir ='/Volumes/GoogleDrive/My Drive/Lake/Atgasuk_UAF/model output',
              modWxFile  ='/Volumes/GoogleDrive/My Drive/Lake/Atgasuk_UAF/model_inputs/AtgasukDate.dat',
              outDir='/Volumes/GoogleDrive/My Drive/beavers/LAKE/Figs/',
              plotPrefix='Atq',
              Depths=c(0.3,2.5),
              dateBreak='4 months',
              dateBreak_minor='1 month',
              P_mult=60*60,
              sec_axis_coeff=10)
##----------------------------------------------------------------
##                            foxden                            --
##----------------------------------------------------------------
#P_mult 60*60 for hourly, 60*60*24 for daily
#sec_axis_coeff for scaling Precip axis and precip data, 50 for daily, 10 for monthly
plotLakeDepth(obsDataFile='/Volumes/GoogleDrive/My Drive/Lake/FoxDen_UAF/lake temp/lakeTemp.csv',
              modDataDir ='/Volumes/GoogleDrive/My Drive/Lake/FoxDen_UAF/model output',
              modWxFile  ='/Volumes/GoogleDrive/My Drive/Lake/FoxDen_UAF/met data/FoxDenDayDate.dat',
              outDir='/Volumes/GoogleDrive/My Drive/beavers/LAKE/Figs/',
              plotPrefix='Fox',
              Depths=c(1.5),
              dateBreak='4 months',
              dateBreak_minor='1 month',
              P_mult=60*60*24,
              sec_axis_coeff=10)
