import os
import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset
import math as mt
from pyproj import Transformer
import geopy.distance
import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as crs
import cartopy.feature as cfeature
from collections import defaultdict
import matplotlib.path as mpath
import matplotlib.patches as patches

models = ['ACCESS-CM2', 'CanESM5', 'CMCC-CM2-SR5', 'CMCC-ESM2', 'CESM2-WACCM', 'CESM2', 'CNRM-CM6-1', 'CNRM-ESM2-1', 'FGOALS-g3', 'EC-Earth3', 'INM-CM5-0', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MIROC6', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL', 'MRI-ESM2-0']
awicm3_members = ['r1', 'r2', 'r3', 'r4', 'r5']

ncdata = nc.Dataset('Suppl_Figure2_3_4.nc')
#Read ERA5 data
lon_era5 = ncdata['ERA5_lon'][:]
lat_era5 = ncdata['ERA5_lat'][:]
sic_era5 = ncdata['ERA5_trends'][0,:,:]
t2m_era5 = ncdata['ERA5_trends'][1,:,:]
w10m_era5 = ncdata['ERA5_trends'][2,:,:]
uwind_era5 = ncdata['ERA5_uv_trends'][0,:,:]
vwind_era5 = ncdata['ERA5_uv_trends'][1,:,:]

#Read CMIP6 multi-model mean data
lon_cmip6 = ncdata['ENSEMBLE_lon'][:]
lat_cmip6 = ncdata['ENSEMBLE_lat'][:]
sic_cmip6 = ncdata['ENSEMBLE_trends'][0,:,:]
t2m_cmip6 = ncdata['ENSEMBLE_trends'][1,:,:]
w10m_cmip6 = ncdata['ENSEMBLE_trends'][2,:,:]
uwind_cmip6 = ncdata['ENSEMBLE_uv_trends'][0,:,:]
vwind_cmip6 = ncdata['ENSEMBLE_uv_trends'][1,:,:]

lon_model = defaultdict(list)
lat_model = defaultdict(list)
sic_model = defaultdict(list)
t2m_model = defaultdict(list)
w10m_model = defaultdict(list)
uwind_model = defaultdict(list)
vwind_model = defaultdict(list)

#Read cmip6 models data
for model in models:
    lon_model[model] = ncdata[f'{model}_lon'][:]
    lat_model[model] = ncdata[f'{model}_lat'][:]
    sic_model[model] = ncdata[f'{model}_trends'][0,:,:]
    t2m_model[model] = ncdata[f'{model}_trends'][1,:,:]
    w10m_model[model] = ncdata[f'{model}_trends'][2,:,:]
    if model in ['ACCESS-CM2', 'CanESM5', 'CMCC-CM2-SR5', 'CMCC-ESM2', 'CNRM-CM6-1', 'CNRM-ESM2-1', 'EC-Earth3', 'INM-CM5-0', 'IPSL-CM6A-LR', 'MIROC6', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'MRI-ESM2-0', 'UKESM1-0-LL']:
        uwind_model[model] = ncdata[f'{model}_uv_trends'][0,:,:]
        vwind_model[model] = ncdata[f'{model}_uv_trends'][1,:,:]
    
#Read awicm3 ensemble members data
for member in awicm3_members:
    model = 'AWI-CM-1-1-MR_%s' %member 
    lon_model[model] = ncdata[f'{model}_lon'][:]
    lat_model[model] = ncdata[f'{model}_lat'][:]
    sic_model[model] = ncdata[f'{model}_trends'][0,:,:]
    t2m_model[model] = ncdata[f'{model}_trends'][1,:,:]
    w10m_model[model] = ncdata[f'{model}_trends'][2,:,:]
    uwind_model[model] = ncdata[f'{model}_uv_trends'][0,:,:]
    vwind_model[model] = ncdata[f'{model}_uv_trends'][1,:,:]

#Plot
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)

#Suppl_Figure2

levels = [-12.5, -12.0, -11.5, -11.0, -10.5, -10.0, -9.5, -9.0, -8.5, -8.0, -7.5, -7.0, -6.5, -6.0, -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5]
vmin = levels[0]
vmax = levels[len(levels)-1]

fig=plt.figure(figsize=(20,25))
ax = plt.subplot(5, 5, 1, projection=crs.NorthPolarStereo())
ax.set_extent([0, 359, 60, 90], crs=crs.PlateCarree())
plt.contourf(lon_era5, lat_era5, sic_era5, levels=levels, transform=crs.PlateCarree(), cmap="RdBu", extend="both")
ax.add_feature(cfeature.LAND, facecolor='darkgray')
ax.coastlines()
ax.set_boundary(circle, transform=ax.transAxes)
plt.title('ERA5', x=0.85, y=-0.05, fontsize=11)
ax = plt.subplot(5, 5, 2, projection=crs.NorthPolarStereo())
ax.set_extent([0, 359, 60, 90], crs=crs.PlateCarree())
plt.contourf(lon_cmip6, lat_cmip6, sic_cmip6, levels=levels, transform=crs.PlateCarree(), cmap="RdBu", extend="both")
plt.title('CMIP6', x=0.85, y=-0.05, fontsize=11)
ax.add_feature(cfeature.LAND, facecolor='darkgray')
ax.coastlines()
ax.set_boundary(circle, transform=ax.transAxes)
axes=[3,4,6,7,8,9,11,12,13,14,16,17,18,19,21,22,23,24]
i=0
for model in models:
    ax = plt.subplot(5, 5, axes[i], projection=crs.NorthPolarStereo())
    ax.set_extent([0, 359, 60, 90], crs=crs.PlateCarree())
    plt.contourf(lon_model[model], lat_model[model], sic_model[model], levels=levels, transform=crs.PlateCarree(), cmap="RdBu", extend="both")
    ax.add_feature(cfeature.LAND, facecolor='darkgray')
    ax.coastlines()
    ax.set_boundary(circle, transform=ax.transAxes)
    plt.title('%s' %model, x=0.85, y=-0.05, fontsize=11)
    i+=1
m=0
awicm3_axes=[5,10,15,20,25]
for member in awicm3_members:
    model = 'AWI-CM-1-1-MR_%s' %member
    ax = plt.subplot(5, 5, awicm3_axes[m], projection=crs.NorthPolarStereo())
    ax.set_extent([0, 359, 60, 90], crs=crs.PlateCarree())
    plt.contourf(lon_model[model], lat_model[model], sic_model[model], levels=levels, transform=crs.PlateCarree(), cmap="RdBu", extend="both")
    ax.add_feature(cfeature.LAND, facecolor='darkgray')
    ax.coastlines()
    ax.set_boundary(circle, transform=ax.transAxes)
    plt.title('%s' %model, x=0.85, y=-0.05, fontsize=11)
    m+=1
rect = patches.Rectangle((0.8, 0.065), 0.2, 0.925, transform=fig.transFigure, linewidth=2, edgecolor='none', facecolor='gray', alpha=0.2, zorder=-1)
fig.add_artist(rect)
plt.tight_layout()
plt.subplots_adjust(bottom=0.05)
cbar_ax = fig.add_axes([0.2, 0.03, 0.6, 0.02])
m = plt.cm.ScalarMappable(cmap="RdBu")
m.set_clim(vmin, vmax)
cbar = fig.colorbar(m, cax=cbar_ax, orientation='horizontal')
cbar_ax.set_xticks([-12, -8, -4, 0, 4, 8, 12])
cbar_ax.set_xticklabels([-12, -8, -4, 0, 4, 8, 12], fontsize=20)
cbar.set_label('Sea ice extent trend over the period 1950-2020 [% dec$^{-1}$]', fontsize=20)
fig.savefig('Suppl_Figure2.pdf')

#Suppl_Figure3

levels = [-2.5, -2.375, -2.25, -2.125, -2.0, -1.875, -1.75, -1.625, -1.5, -1.375, -1.25, -1.125, -1.0, -0.875, -0.75, -0.625, -0.5, -0.375, -0.25, -0.125, 0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 2.0, 2.125, 2.25, 2.375, 2.5]
vmin = levels[0]
vmax = levels[len(levels)-1]

fig=plt.figure(figsize=(20,25))
ax = plt.subplot(5, 5, 1, projection=crs.NorthPolarStereo())
ax.set_extent([0, 359, 60, 90], crs=crs.PlateCarree())
plt.contourf(lon_era5, lat_era5, t2m_era5, levels=levels, transform=crs.PlateCarree(), cmap="RdBu_r", extend="both")
ax.add_feature(cfeature.LAND, facecolor='darkgray')
ax.coastlines()
ax.set_boundary(circle, transform=ax.transAxes)
plt.title('ERA5', x=0.85, y=-0.05, fontsize=11)
ax = plt.subplot(5, 5, 2, projection=crs.NorthPolarStereo())
ax.set_extent([0, 359, 60, 90], crs=crs.PlateCarree())
plt.contourf(lon_cmip6, lat_cmip6, t2m_cmip6, levels=levels, transform=crs.PlateCarree(), cmap="RdBu_r", extend="both")
plt.title('CMIP6', x=0.85, y=-0.05, fontsize=11)
ax.add_feature(cfeature.LAND, facecolor='darkgray')
ax.coastlines()
ax.set_boundary(circle, transform=ax.transAxes)
axes=[3,4,6,7,8,9,11,12,13,14,16,17,18,19,21,22,23,24]
i=0
for model in models:
    ax = plt.subplot(5, 5, axes[i], projection=crs.NorthPolarStereo())
    ax.set_extent([0, 359, 60, 90], crs=crs.PlateCarree())
    plt.contourf(lon_model[model], lat_model[model], t2m_model[model], levels=levels, transform=crs.PlateCarree(), cmap="RdBu_r", extend="both")
    ax.add_feature(cfeature.LAND, facecolor='darkgray')
    ax.coastlines()
    ax.set_boundary(circle, transform=ax.transAxes)
    plt.title('%s' %model, x=0.85, y=-0.05, fontsize=11)
    i+=1
m=0
awi_axes=[5,10,15,20,25]
for member in awicm3_members:
    model = 'AWI-CM-1-1-MR_%s' %member
    ax = plt.subplot(5, 5, awi_axes[m], projection=crs.NorthPolarStereo())
    ax.set_extent([0, 359, 60, 90], crs=crs.PlateCarree())
    plt.contourf(lon_model[model], lat_model[model], t2m_model[model], levels=levels, transform=crs.PlateCarree(), cmap="RdBu_r", extend="both")
    ax.add_feature(cfeature.LAND, facecolor='darkgray')
    ax.coastlines()
    ax.set_boundary(circle, transform=ax.transAxes)
    plt.title('%s' %model, x=0.85, y=-0.05, fontsize=11)
    m+=1
rect = patches.Rectangle((0.8, 0.065), 0.2, 0.925, transform=fig.transFigure, linewidth=2, edgecolor='none', facecolor='gray', alpha=0.2, zorder=-1)
fig.add_artist(rect)
plt.tight_layout()
plt.subplots_adjust(bottom=0.05)
cbar_ax = fig.add_axes([0.2, 0.03, 0.6, 0.02])
m = plt.cm.ScalarMappable(cmap="RdBu_r")
m.set_clim(vmin, vmax)
cbar = fig.colorbar(m, cax=cbar_ax, orientation='horizontal')
cbar_ax.set_xticks([-2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5])
cbar_ax.set_xticklabels([-2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5], fontsize=18)
cbar.set_label('2m temperature trend over the period 1950-2020 [K dec$^{-1}$]', fontsize=20)
fig.savefig('Suppl_Figure3.pdf')

#Suppl_Figure4

levels=[-0.350, -0.325, -0.300, -0.275, -0.250, -0.225, -0.200, -0.175, -0.150, -0.125, -0.100, -0.075, -0.050, -0.025, 0.000, 0.025, 0.050, 0.075, 0.100, 0.125, 0.150, 0.175, 0.200, 0.225, 0.250, 0.275, 0.300, 0.325, 0.350]
vmin = levels[0]
vmax = levels[len(levels)-1]

fig=plt.figure(figsize=(20,25))
ax = plt.subplot(5, 5, 1, projection=crs.NorthPolarStereo())
ax.set_extent([0, 359, 60, 90], crs=crs.PlateCarree())
plt.contourf(lon_era5, lat_era5, w10m_era5, levels=levels, transform=crs.PlateCarree(), cmap="RdBu_r", extend="both")
q=plt.quiver(lon_era5, lat_era5, uwind_era5, vwind_era5, pivot='tail', transform=crs.PlateCarree(), regrid_shape=30, units='inches', scale=0.5)
plt.title('ERA5', x=0.85, y=-0.05, fontsize=11)
ax.add_feature(cfeature.LAND, facecolor='darkgray')
ax.coastlines()
ax.set_boundary(circle, transform=ax.transAxes)
ax = plt.subplot(5, 5, 2, projection=crs.NorthPolarStereo())
ax.set_extent([0, 359, 60, 90], crs=crs.PlateCarree())
plt.contourf(lon_cmip6, lat_cmip6, w10m_cmip6, levels=levels, transform=crs.PlateCarree(), cmap="RdBu_r", extend="both")
q=plt.quiver(lon_cmip6, lat_cmip6, uwind_cmip6, vwind_cmip6, pivot='tail', transform=crs.PlateCarree(), regrid_shape=30, units='inches', scale=0.5)
plt.title('CMIP6', x=0.85, y=-0.05, fontsize=11)
ax.add_feature(cfeature.LAND, facecolor='darkgray')
ax.coastlines()
ax.set_boundary(circle, transform=ax.transAxes)
axes=[3,4,6,7,8,9,11,12,13,14,16,17,18,19,21,22,23,24]
i=0
for model in models:
    if model in ['CESM2', 'CESM2-WACCM', 'FGOALS-g3', 'NorESM2-LM']:
        ax = plt.subplot(5, 5, axes[i], projection=crs.NorthPolarStereo())
        ax.set_extent([0, 359, 60, 90], crs=crs.PlateCarree())
        plt.contourf(lon_model[model], lat_model[model], w10m_model[model], levels=levels, transform=crs.PlateCarree(), cmap="RdBu_r", extend="both")
        ax.add_feature(cfeature.LAND, facecolor='darkgray')
        ax.coastlines()
        ax.set_boundary(circle, transform=ax.transAxes)
        plt.title('%s' %model, x=0.85, y=-0.05, fontsize=11)
        i+=1
    else:        
        ax = plt.subplot(5, 5, axes[i], projection=crs.NorthPolarStereo())
        ax.set_extent([0, 359, 60, 90], crs=crs.PlateCarree())
        plt.contourf(lon_model[model], lat_model[model], w10m_model[model], levels=levels, transform=crs.PlateCarree(), cmap="RdBu_r", extend="both")
        q=plt.quiver(lon_model[model], lat_model[model], uwind_model[model], vwind_model[model], pivot='tail', transform=crs.PlateCarree(), regrid_shape=30, units='inches', scale=0.5)
        ax.add_feature(cfeature.LAND, facecolor='darkgray')
        ax.coastlines()
        ax.set_boundary(circle, transform=ax.transAxes)
        plt.title('%s' %model, x=0.85, y=-0.05, fontsize=11)
        i+=1
m=0
awi_axes=[5,10,15,20,25]
for member in awicm3_members:
    model = 'AWI-CM-1-1-MR_%s' %member
    ax = plt.subplot(5, 5, awi_axes[m], projection=crs.NorthPolarStereo())
    ax.set_extent([0, 359, 60, 90], crs=crs.PlateCarree())
    plt.contourf(lon_model[model], lat_model[model], w10m_model[model], levels=levels, transform=crs.PlateCarree(), cmap="RdBu_r", extend="both")
    q=plt.quiver(lon_model[model], lat_model[model], uwind_model[model], vwind_model[model], pivot='tail', transform=crs.PlateCarree(), regrid_shape=30, units='inches', scale=0.5)
    ax.add_feature(cfeature.LAND, facecolor='darkgray')
    ax.coastlines()
    ax.set_boundary(circle, transform=ax.transAxes)
    plt.title('%s' %model, x=0.85, y=-0.05, fontsize=11)
    m+=1
rect = patches.Rectangle((0.8, 0.065), 0.2, 0.925, transform=fig.transFigure, linewidth=2, edgecolor='none', facecolor='gray', alpha=0.2, zorder=-1)
fig.add_artist(rect)
plt.tight_layout()
plt.subplots_adjust(bottom=0.05)
cbar_ax = fig.add_axes([0.2, 0.03, 0.6, 0.02])
m = plt.cm.ScalarMappable(cmap="RdBu_r")
m.set_clim(vmin, vmax)
cbar = fig.colorbar(m, cax=cbar_ax, orientation='horizontal')
cbar_ax.set_xticks([-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3])
cbar_ax.set_xticklabels([-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3], fontsize=18)
cbar.set_label('10m wind speed trend over the period 1950-2020 [m/s dec$^{-1}$]', fontsize=20)
plt.quiverkey(q, 0.85,0.03, 0.1, '0.1 m/s dec$^{-1}$', coordinates = 'figure', fontproperties={'size' : 12})
fig.savefig('Suppl_Figure4.pdf')