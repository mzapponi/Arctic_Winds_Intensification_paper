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
import geopy.distance
import matplotlib.path as mpath
import matplotlib.patches as mpatches


ncdata = nc.Dataset('Figure1.nc')
#Read ERA5 data
lon_era5 = ncdata['lon_era5'][:]
lat_era5 = ncdata['lat_era5'][:]
sic_era5 = ncdata['sic_trend_era5'][:,:]
t2m_era5 = ncdata['t2m_trend_era5'][:,:]
w10m_era5 = ncdata['w10m_trend_era5'][:,:]
#Read CMIP6 data
lon_cmip6 = ncdata['lon_CMIP6'][:]
lat_cmip6 = ncdata['lat_CMIP6'][:]
sic_cmip6_hist = ncdata['sic_trend_CMIP6'][0,:,:]
t2m_cmip6_hist = ncdata['t2m_trend_CMIP6'][0,:,:]
w10m_cmip6_hist = ncdata['w10m_trend_CMIP6'][0,:,:]
sic_cmip6_tot = ncdata['sic_trend_CMIP6'][1,:,:]
t2m_cmip6_tot = ncdata['t2m_trend_CMIP6'][1,:,:]
w10m_cmip6_tot = ncdata['w10m_trend_CMIP6'][1,:,:]
sic_mask_hist = ncdata['sic_signmask'][0,:,:]
t2m_mask_hist = ncdata['t2m_signmask'][0,:,:]
w10m_mask_hist = ncdata['w10m_signmask'][0,:,:]
sic_mask_tot = ncdata['sic_signmask'][1,:,:]
t2m_mask_tot = ncdata['t2m_signmask'][1,:,:]
w10m_mask_tot = ncdata['w10m_signmask'][1,:,:]
domain = ncdata['domain_mask'][:,:]


theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)

#Figure1#

#Sea ice concentration
levels = np.arange(-12.0, 12.1, 0.5)
vmin = levels[0]
vmax = levels[len(levels)-1]

fig=plt.figure(figsize=(15,5))
ax = plt.subplot(1, 3, 1, projection=crs.NorthPolarStereo())
ax.set_extent([0, 359, 60, 90], crs=crs.PlateCarree())
plt.contourf(lon_era5, lat_era5, sic_era5, levels=levels, transform=crs.PlateCarree(), cmap="RdBu", extend="both")
ax.add_feature(cfeature.LAND, facecolor='darkgray')
ax.coastlines()
ax.set_boundary(circle, transform=ax.transAxes)
plt.title('ERA5 (1950-2020)', x=0.95, y=0.0, fontsize=14)
ax = plt.subplot(1, 3, 2, projection=crs.NorthPolarStereo())
ax.set_extent([0, 359, 60, 90], crs=crs.PlateCarree())
plt.contourf(lon_cmip6, lat_cmip6, sic_cmip6_hist, levels=levels, transform=crs.PlateCarree(), cmap="RdBu", extend="both")
plt.contourf(lon_cmip6, lat_cmip6, sic_mask_hist, levels=[0.5, 1], transform=crs.PlateCarree(), colors='none', hatches=['..'])
ax.add_feature(cfeature.LAND, facecolor='darkgray')
ax.coastlines()
ax.set_boundary(circle, transform=ax.transAxes)
plt.title('CMIP6 (1950-2020)', x=0.95, y=0.0, fontsize=14)
ax = plt.subplot(1, 3, 3, projection=crs.NorthPolarStereo())
ax.set_extent([0, 359, 60, 90], crs=crs.PlateCarree())
plt.contourf(lon_cmip6, lat_cmip6, sic_cmip6_tot, levels=levels, transform=crs.PlateCarree(), cmap="RdBu", extend="both")
plt.contourf(lon_cmip6, lat_cmip6, sic_mask_tot, levels=[0.5, 1], transform=crs.PlateCarree(), colors='none', hatches=['..'])
ax.add_feature(cfeature.LAND, facecolor='darkgray')
ax.coastlines()
ax.set_boundary(circle, transform=ax.transAxes)
plt.title('CMIP6 (1950-2100)', x=0.95, y=0.0, fontsize=14)
plt.tight_layout()
plt.subplots_adjust(top=0.95, bottom=0.05, left=0.02, right=0.85)
cbar_ax = fig.add_axes([0.9, 0.2, 0.012, 0.5])
m = plt.cm.ScalarMappable(cmap="RdBu")
m.set_clim(vmin, vmax)
cbar = fig.colorbar(m, cax=cbar_ax)
cbar.set_label('% dec$^{-1}$', fontsize=12)
fig.text(0.04, 0.8, 'a', fontsize=12, fontweight='bold')
fig.text(0.32, 0.8, 'b', fontsize=12, fontweight='bold')
fig.text(0.6, 0.8, 'c', fontsize=12, fontweight='bold')
fig.savefig('Figure1_a-c.pdf')

#2-meter temperature#
levels = np.linspace(-2.0, 2.0, num=81)
vmin = levels[0]
vmax = levels[len(levels)-1]

fig=plt.figure(figsize=(15,5))
ax = plt.subplot(1, 3, 1, projection=crs.NorthPolarStereo())
ax.set_extent([0, 360, 60, 90], crs=crs.PlateCarree())
plt.contourf(lon_era5, lat_era5, t2m_era5, levels=levels, transform=crs.PlateCarree(), cmap="RdBu_r", extend="both")
ax.add_feature(cfeature.LAND, facecolor='darkgray')
ax.coastlines()
ax.set_boundary(circle, transform=ax.transAxes)
plt.title('ERA5 (1950-2020)', x=0.95, y=0.0, fontsize=14)
ax = plt.subplot(1, 3, 2, projection=crs.NorthPolarStereo())
ax.set_extent([0, 360, 60, 90], crs=crs.PlateCarree())
plt.contourf(lon_cmip6, lat_cmip6, t2m_cmip6_hist, levels=levels, transform=crs.PlateCarree(), cmap="RdBu_r", extend="both")
plt.contourf(lon_cmip6, lat_cmip6, t2m_mask_hist, levels=[0.5, 1], transform=crs.PlateCarree(), colors='none', hatches=['..'])
plt.contour(lon_cmip6, lat_cmip6+0.625, domain, levels=[0.5,1,1.5], transform=crs.PlateCarree(), colors='cyan', linewidths=2, zorder=3)
ax.add_feature(cfeature.LAND, facecolor='darkgray')
ax.coastlines()
ax.set_boundary(circle, transform=ax.transAxes)
gl = ax.gridlines(draw_labels=True, crs=crs.PlateCarree(), linewidth=1, y_inline=False, color='gray', zorder=2)
gl.ylines = False
plt.title('CMIP6 (1950-2020)', x=0.95, y=0.0, fontsize=14)
ax = plt.subplot(1, 3, 3, projection=crs.NorthPolarStereo())
ax.set_extent([0, 360, 60, 90], crs=crs.PlateCarree())
plt.contourf(lon_cmip6, lat_cmip6, t2m_cmip6_tot, levels=levels, transform=crs.PlateCarree(), cmap="RdBu_r", extend="both")
plt.contourf(lon_cmip6, lat_cmip6, t2m_mask_tot, levels=[0.5, 1], transform=crs.PlateCarree(), colors='none', hatches=['..'])
ax.add_feature(cfeature.LAND, facecolor='darkgray')
ax.coastlines()
ax.set_boundary(circle, transform=ax.transAxes)
plt.title('CMIP6 (1950-2100)', x=0.95, y=0.0, fontsize=14)
plt.tight_layout()
plt.subplots_adjust(top=0.95, bottom=0.05, left=0.02, right=0.85)
cbar_ax = fig.add_axes([0.9, 0.2, 0.012, 0.5])
m = plt.cm.ScalarMappable(cmap="RdBu_r")
m.set_clim(vmin, vmax)
cbar = fig.colorbar(m, cax=cbar_ax)
cbar.set_label('K dec$^{-1}$', fontsize=12)
fig.text(0.04, 0.8, 'd', fontsize=12, fontweight='bold')
fig.text(0.32, 0.8, 'e', fontsize=12, fontweight='bold')
fig.text(0.6, 0.8, 'f', fontsize=12, fontweight='bold')
fig.savefig('Figure1_d-f.pdf')

#10-meter wind speed 
levels = np.linspace(-0.25, 0.25, num=41)
vmin = levels[0]
vmax = levels[len(levels)-1]

fig=plt.figure(figsize=(15,5))
ax = plt.subplot(1, 3, 1, projection=crs.NorthPolarStereo())
ax.set_extent([0, 359, 60, 90], crs=crs.PlateCarree())
plt.contourf(lon_era5, lat_era5, w10m_era5, levels=levels, transform=crs.PlateCarree(), cmap="RdBu_r", extend="both")
ax.add_feature(cfeature.LAND, facecolor='darkgray')
ax.coastlines()
ax.set_boundary(circle, transform=ax.transAxes)
plt.title('ERA5 (1950-2020)', x=0.95, y=0.0, fontsize=14)
ax = plt.subplot(1, 3, 2, projection=crs.NorthPolarStereo())
ax.set_extent([0, 359, 60, 90], crs=crs.PlateCarree())
plt.contourf(lon_cmip6, lat_cmip6, w10m_cmip6_hist, levels=levels, transform=crs.PlateCarree(), cmap="RdBu_r", extend="both")
plt.contourf(lon_cmip6, lat_cmip6, w10m_mask_hist, levels=[0.5, 1], transform=crs.PlateCarree(), colors='none', hatches=['..'])
ax.add_feature(cfeature.LAND, facecolor='darkgray')
ax.coastlines()
ax.set_boundary(circle, transform=ax.transAxes)
plt.title('CMIP6 (1950-2020)', x=0.95, y=0.0, fontsize=14)
ax = plt.subplot(1, 3, 3, projection=crs.NorthPolarStereo())
ax.set_extent([0, 359, 60, 90], crs=crs.PlateCarree())
plt.contourf(lon_cmip6, lat_cmip6, w10m_cmip6_tot, levels=levels, transform=crs.PlateCarree(), cmap="RdBu_r", extend="both")
plt.contourf(lon_cmip6, lat_cmip6, w10m_mask_tot, levels=[0.5, 1], transform=crs.PlateCarree(), colors='none', hatches=['..'])
ax.add_feature(cfeature.LAND, facecolor='darkgray')
ax.coastlines()
ax.set_boundary(circle, transform=ax.transAxes)
plt.title('CMIP6 (1950-2100)', x=0.95, y=0.0, fontsize=14)
plt.tight_layout()
plt.subplots_adjust(top=0.95, bottom=0.05, left=0.02, right=0.85)
cbar_ax = fig.add_axes([0.9, 0.2, 0.012, 0.5])
m = plt.cm.ScalarMappable(cmap="RdBu_r")
m.set_clim(vmin, vmax)
cbar = fig.colorbar(m, cax=cbar_ax)
cbar.set_label('m/s dec$^{-1}$', fontsize=12)
fig.text(0.04, 0.8, 'g', fontsize=12, fontweight='bold')
fig.text(0.32, 0.8, 'h', fontsize=12, fontweight='bold')
fig.text(0.6, 0.8, 'i', fontsize=12, fontweight='bold')
fig.savefig('Figure1_g-i.pdf')