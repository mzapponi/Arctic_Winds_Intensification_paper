import os
import numpy as np
import netCDF4 as nc
import math as mt
from pyproj import Transformer
import datetime
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.colors as mcolors
import cartopy.crs as crs
import cartopy.feature as cfeature
from collections import defaultdict
import geopy.distance
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from scipy.stats import pearsonr

seasons = 2
decades_ERA5 = 7
decades_CMIP6 = 15
dt_seasonmean_obs = np.zeros((seasons,decades_ERA5))
rw_seasonmean_obs = np.zeros((seasons,decades_ERA5))
dt_seasonmean_ens = np.zeros((seasons,decades_CMIP6))
rw_seasonmean_ens = np.zeros((seasons,decades_CMIP6))

dt_seasonmean_obs[0] = np.loadtxt('Figure4.txt', skiprows=2, delimiter=',', max_rows=1) #JFM
dt_seasonmean_obs[1] = np.loadtxt('Figure4.txt', skiprows=3, delimiter=',', max_rows=1) #JAS
rw_seasonmean_obs[0] = np.loadtxt('Figure4.txt', skiprows=5, delimiter=',', max_rows=1) #JFM
rw_seasonmean_obs[1] = np.loadtxt('Figure4.txt', skiprows=6, delimiter=',', max_rows=1) #JAS

dt_seasonmean_ens[0] = np.loadtxt('Figure4.txt', skiprows=9, delimiter=',', max_rows=1) #JFM
dt_seasonmean_ens[1] = np.loadtxt('Figure4.txt', skiprows=10, delimiter=',', max_rows=1) #JAS
rw_seasonmean_ens[0] = np.loadtxt('Figure4.txt', skiprows=12, delimiter=',', max_rows=1) #JFM
rw_seasonmean_ens[1] = np.loadtxt('Figure4.txt', skiprows=13, delimiter=',', max_rows=1) #JAS

color_hist = plt.cm.Blues(np.linspace(0,1,9))
color_fut = plt.cm.Oranges(np.linspace(0,1,10))
color_tot = np.vstack((color_hist[2:], color_fut[2:]))
mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', color_tot)

x = np.zeros((2))
y = [-1, 4]
plt.rcParams["axes.prop_cycle"] = plt.cycler("color", color_hist[2:])
winter= mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=10, label='Winter')
summer = mlines.Line2D([], [], color='k', marker='s', linestyle='None', markersize=10, label='Summer')
fig = plt.figure(figsize=(12,4))
markers = ['o', 's']
ax = plt.subplot(1,2,1)
ax.plot(x,y, color='k', linestyle='--')
for season in range(seasons):
    for decade in range(decades_ERA5):
        plt.plot(dt_seasonmean_obs[season,decade], rw_seasonmean_obs[season,decade]*100, marker=markers[season])
plt.ylim([60, 80])
plt.xlim([-5, 9])
plt.xticks([-4, -2, 0, 2, 4, 6, 8], fontsize=16)
ax.set_xlabel('T$_{850hPa}$-T$_{2m}$ (K)', fontsize=18)
plt.yticks([60, 70, 80], fontsize=14)
ax.set_ylabel('W$_{10m}$/W$_{850hPa}$ (%)', fontsize=18)
plt.title('ERA5', y=0.88, x=0.85, fontsize=20)
plt.legend(handles=[winter, summer], loc='lower left')
plt.rcParams["axes.prop_cycle"] = plt.cycler("color", color_tot)
plt.text(0.03, 0.99, 'a', ha='right', va='top', transform=plt.gca().transAxes, fontweight='bold')
ax = plt.subplot(1,2,2)
ax.plot(x,y, color='k', linestyle='--')
for season in range(seasons):
    for decade in range(decades_CMIP6):
        plt.plot(dt_seasonmean_ens[season,decade], rw_seasonmean_ens[season,decade]*100, marker=markers[season])
plt.ylim([60, 80])
plt.xlim([-5, 9])
plt.title('CMIP6', y=0.88, x=0.85, fontsize=20)
plt.xticks([-4, -2, 0, 2, 4, 6, 8], fontsize=16)
ax.set_xlabel('T$_{850hPa}$-T$_{2m}$ (K)', fontsize=18)
plt.yticks([60, 70, 80],[], fontsize=14)
plt.text(0.03, 0.99, 'b', ha='right', va='top', transform=plt.gca().transAxes, fontweight='bold')
plt.tight_layout()
plt.subplots_adjust(top=0.95, right=0.85)
bounds = np.arange(16)
norm = mcolors.BoundaryNorm(bounds, mymap.N, extend='both')
cbar_ax = fig.add_axes([0.9, 0.2, 0.02, 0.7])
cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=mymap), cax=cbar_ax, orientation='vertical', ticks=[0, 7, 15])
cbar.ax.tick_params(labelsize=14)
cbar.set_ticklabels([1950, 2020, 2100])
fig.savefig('Figure4.pdf')
