"""
19/11/2018 - JFE
this scripts plots a map with the difference between observed AGB and the
AGB with LUH from 1850
"""

import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.geoaxes import GeoAxes
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import numpy as np
from mpl_toolkits.axes_grid1 import AxesGrid
from useful import *

#get areas and landmask
areas = get_areas()
_,landmask = get_predictors()

#load current AGB
med = xr.open_rasterio('/disk/scratch/local.2/jexbraya/AGB/Avitable_AGB_Map_0.25d.tif')
med.values[med.values == med.nodatavals[0]] = np.nan
med[0].values[~landmask] = np.nan

#load uncertainty to write some stats
unc = xr.open_rasterio('/disk/scratch/local.2/jexbraya/AGB/Avitable_AGB_Uncertainty_0.25d.tif')

#load agb with past land use
pot = xr.open_dataset('/disk/scratch/local.2/jexbraya/pantrop-AGB-LUH/output/AGB_hist.nc')

#print some statistics
lvls = ['mean','lower','upper']
for aa, agb in enumerate([med,med-unc,med+unc]):
    lvl = lvls[aa]
    agb.values[agb.values<0] = 0
    agb_tot = (agb[0]*areas).sum()*1e-13*.48
    print('AGB',lvl,np.round(agb_tot.values,1),' Pg C')

    pot_tot = (pot['AGB_' + lvl][:10].mean(axis=0)*areas).sum()*1e-13*.48
    print('Pot',lvl,np.round(pot_tot.values,1),' Pg C')

####
#### Plotting the map
####

#calculate dAGB to plot
dAGB = med[0].copy()
dAGB.values = pot.AGB_mean[:10].values.mean(axis=0) - med[0].values

#create a figure using the axesgrid to make the colorbar fit on the axis
projection = ccrs.PlateCarree()
axes_class = (GeoAxes,dict(map_projection=projection))

#create figure
fig = plt.figure('fig_3maps',figsize=(10,10))
fig.clf()

#create axes grid
axgr = AxesGrid(fig,111,nrows_ncols=(3,1),axes_class=axes_class,label_mode='',cbar_mode='each',cbar_pad = 0.25,cbar_size="3%",axes_pad=.5)

#plot
vmaxes = [200,200,50]
cmaps = ['viridis','viridis','YlOrRd']
titles = ['a) AGB$_{2000s}$','b) AGB$_{1850s}$','c) AGB$_{1850s}$ - AGB$_{2000s}$']
for mm,map2plot in enumerate([med[0],pot.AGB_mean[:10].mean(axis=0),dAGB]):
    (map2plot*.48).plot.imshow(ax=axgr[mm],cbar_ax=axgr.cbar_axes[mm],vmin=0,vmax=vmaxes[mm],extend='max',
                        interpolation='nearest',cbar_kwargs={'label':'Mg C ha$^{-1}$'},
                        cmap=cmaps[mm],xticks=np.arange(-120,161,40),yticks=np.arange(-60,41,20),
                        add_labels=False,ylim=(-60,40),xlim=(-120,160))

    #add grey mask for land regions outside the study, and black for the oceans
    axgr[mm].add_feature(cfeat.LAND,zorder=-1,facecolor='silver')
    axgr[mm].add_feature(cfeat.OCEAN,zorder=-1,facecolor='gray')

    #set labels
    axgr[mm].yaxis.set_major_formatter(LatitudeFormatter())
    axgr[mm].xaxis.set_major_formatter(LongitudeFormatter())

    axgr[mm].text(0.98,0.02,titles[mm],transform=axgr[mm].transAxes,ha='right',va='bottom',weight='bold')

fig.show()
fig.savefig('fig_3maps.png',bbox_inches='tight',dpi=300)
