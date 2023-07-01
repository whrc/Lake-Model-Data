import h5py
import pandas as pd
import datetime
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal
import matplotlib.colors as colors

class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))
#####
#levels=np.linspace(-4,25,21)
#cntr1 = ax.contourf(np.flipud(wtemps),levels,norm=MidpointNormalize(midpoint=0),cmap="bwr")

def plot_water_temps_bwr(fdir):
    m=np.loadtxt(u"results/"+fdir+"/time_series/water_temp  1  1.dat",skiprows=6)
    [nx,ny]=np.shape(m)
    #print ('shape of the results:',nx,ny)
    nyy=int((ny-6)/2)
    #print ('number of water depth points:',nyy)

    wdepth=[]
    wtemps=[]
    for i in range(nyy+1):
        wdepth.append(m[0,5+2*i])
        wtemps.append(m[:,6+2*i])
    #print ('shape of the wtemps matrix:',np.shape(wtemps))
    #print ('water depths:',np.round(wdepth,2))

    #finding min and max for the colorbar lvl
    wtemps=np.array(wtemps)
    fstemps=wtemps.flatten()
    #levels=np.linspace(min(fstemps),max(fstemps),21)
    levels=np.linspace(-4,25,21)

    #plotting coutour plot of the water temperatures
    fig, ax = plt.subplots(figsize=(10, 6))
    #nn=np.normalize(np.flipud(wtemps))
    cntr1 = ax.contourf(np.flipud(wtemps),levels,norm=MidpointNormalize(midpoint=0),cmap="bwr")
    
    #cntr1=ax.pcolor(np.flipud(wtemps),norm=MidpointNormalize(midpoint=5),cmap="bwr")   # Set midpoint as 0
    plt.colorbar(cntr1, ax=ax) # To extend colorbar in the min values

    # setting ticks for the countour plot
    plt.title('Water temperature [$^oC$]', fontsize=16);
    plt.xlabel('Time [h]',fontsize=14)
    plt.ylabel('Depth [m]',fontsize=14)
    yticks = np.arange(len(wdepth)-1,-1,-1)
    ax.set_yticks(yticks)
    ax.set_yticklabels(wdepth, fontsize=14);
    
def plot_temps(fdir,ftemp):
    m=np.loadtxt(u"results/"+fdir+"/time_series/"+ftemp+"  1  1.dat",skiprows=6)
    [nx,ny]=np.shape(m)
    nyy=int((ny-6)/2)
    #print ('number of soil depth points:',nyy)

    sdepth=[]
    stemps=[]
    for i in range(nyy+1):
        sdepth.append(m[0,5+2*i])
        stemps.append(m[:,6+2*i])
    #print ('shape of the stemps matrix:',np.shape(stemps))

    #finding min and max for the colorbar lvl
    stemps=np.array(stemps)
    ftemps=stemps.flatten()
    levels=np.linspace(min(ftemps),max(ftemps),10)
    #print(levels)

    #plotting coutour plot of the water temperatures
    fig, ax = plt.subplots(figsize=(10, 8))
    cntr1 = ax.contourf(np.flipud(stemps),levels,cmap="jet")

    fig.colorbar(cntr1, ax=ax)

    # setting ticks for the countour plot
    ftitle='Soil temperature [$^oC$]'
    if ftemp=='water_temp':
        ftitle='Water temperature [$^oC$]'
    plt.title(ftitle, fontsize=16);
    plt.xlabel('Time [h]',fontsize=14)
    plt.ylabel('Depth [m]',fontsize=14)
    yticks = np.arange(len(sdepth)-1,-1,-1)
    ax.set_yticks(yticks)
    ax.set_yticklabels(sdepth, fontsize=14);

def plot_snow_ice(fdir):
    fig, ax = plt.subplots(figsize=(7, 5))
    m=np.loadtxt(u"results/"+fdir+"/time_series/layers  1  1.dat",skiprows=19)
    plt.plot(m[:,14],linewidth=2) 
    plt.plot(m[:,15],linewidth=2)
    plt.legend(['ice','snow'],fontsize=18);
    plt.title('snow and ice layers thickness',fontsize=14)
    plt.xlabel('Time [h]',fontsize=14)
    plt.ylabel('Depth [m]',fontsize=14)