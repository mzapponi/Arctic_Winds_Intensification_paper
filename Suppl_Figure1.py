import os
import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset
import math as mt
from pyproj import Transformer
import datetime
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.colors as mcolors
import matplotlib.path as mpath
import matplotlib.patches as mpatches

years_era5 = 71
years_tot = 151
years_79 = 43

#Read ERA5 data
t2m_era5_global = np.zeros((years_era5))
sic_era5 = np.zeros((years_era5))
t2m_era5 = np.zeros((years_era5))
w10m_era5 = np.zeros((years_era5))
t850_era5 = np.zeros((years_era5))
w850_era5 = np.zeros((years_era5))
dt_era5 = np.zeros((years_era5))
wr_era5 = np.zeros((years_era5))

t2m_era5_global = np.loadtxt('Suppl_Figure1.txt', skiprows=2, delimiter=',', max_rows=1)
sic_era5 = np.loadtxt('Suppl_Figure1.txt', skiprows=4, delimiter=',', max_rows=1)
t2m_era5 = np.loadtxt('Suppl_Figure1.txt', skiprows=6, delimiter=',', max_rows=1)
w10m_era5 = np.loadtxt('Suppl_Figure1.txt', skiprows=8, delimiter=',', max_rows=1)
t850_era5 = np.loadtxt('Suppl_Figure1.txt', skiprows=10, delimiter=',', max_rows=1)
w850_era5 = np.loadtxt('Suppl_Figure1.txt', skiprows=12, delimiter=',', max_rows=1)
dt_era5 = np.loadtxt('Suppl_Figure1.txt', skiprows=14, delimiter=',', max_rows=1)
wr_era5 = np.loadtxt('Suppl_Figure1.txt', skiprows=16, delimiter=',', max_rows=1)

#Read CMIP6 multi-model mean data
sic_cmip6 = np.zeros((years_tot))
t2m_cmip6 = np.zeros((years_tot))
w10m_cmip6 = np.zeros((years_tot))
t850_cmip6 = np.zeros((years_tot))
w850_cmip6 = np.zeros((years_tot))
dt_cmip6 = np.zeros((years_tot))
wr_cmip6 = np.zeros((years_tot))

sic_cmip6 = np.loadtxt('Suppl_Figure1.txt', skiprows=19, delimiter=',', max_rows=1)
t2m_cmip6 = np.loadtxt('Suppl_Figure1.txt', skiprows=21, delimiter=',', max_rows=1)
w10m_cmip6 = np.loadtxt('Suppl_Figure1.txt', skiprows=23, delimiter=',', max_rows=1)
t850_cmip6 = np.loadtxt('Suppl_Figure1.txt', skiprows=25, delimiter=',', max_rows=1)
w850_cmip6 = np.loadtxt('Suppl_Figure1.txt', skiprows=27, delimiter=',', max_rows=1)
dt_cmip6 = np.loadtxt('Suppl_Figure1.txt', skiprows=29, delimiter=',', max_rows=1)
wr_cmip6 = np.loadtxt('Suppl_Figure1.txt', skiprows=31, delimiter=',', max_rows=1)

x = np.arange(years_era5)
x_tot = np.arange(years_tot)

t2m_m_era5 = np.zeros(1)
t2m_c_era5 = np.zeros(1)
w10m_m_era5 = np.zeros(1)
w10m_c_era5 = np.zeros(1)
t850_m_era5 = np.zeros(1)
t850_c_era5 = np.zeros(1)
w850_m_era5 = np.zeros(1)
w850_c_era5 = np.zeros(1)
sic_m_era5 = np.zeros(1)
sic_c_era5 = np.zeros(1)
dt_m_era5 = np.zeros(1)
dt_c_era5 = np.zeros(1)
wr_m_era5 = np.zeros(1)
wr_c_era5 = np.zeros(1)
    
fit = np.polyfit(x,t2m_era5,1, full=True)
t2m_m_era5 = fit [0][0]
t2m_c_era5 = fit [0][1]
fit = np.polyfit(x,t850_era5,1, full=True)
t850_m_era5 = fit [0][0]
t850_c_era5 = fit [0][1]
fit = np.polyfit(x,w10m_era5,1, full=True)
w10m_m_era5 = fit [0][0]
w10m_c_era5 = fit [0][1]
fit = np.polyfit(x,w850_era5,1, full=True)
w850_m_era5 = fit [0][0]
w850_c_era5  = fit [0][1]
fit = np.polyfit(x,sic_era5,1, full=True)
sic_m_era5 = fit [0][0]
sic_c_era5 = fit [0][1]
fit = np.polyfit(x,dt_era5,1, full=True)
dt_m_era5 = fit [0][0]
dt_c_era5 = fit [0][1]
fit = np.polyfit(x,wr_era5,1, full=True)
wr_m_era5 = fit [0][0]
wr_c_era5 = fit [0][1]


x_79=np.arange(years_era5-28)
w10m_m_era5_79 = np.zeros(1)
w10m_c_era5_79 = np.zeros(1)
fit = np.polyfit(x_79,w10m_era5[28:],1, full=True)
w10m_m_era5_79 = fit [0][0]
w10m_c_era5_79 = fit [0][1]

m_globera5 = np.zeros((1))
c_globera5 = np.zeros((1))
fit = np.polyfit(x, t2m_era5_global, deg=1, full=True)
m_globera5 = fit[0][0]
c_globera5 = fit[0][1]


t2m_m_cmip6 = np.zeros(1)
t2m_c_cmip6 = np.zeros(1)
w10m_m_cmip6 = np.zeros(1)
w10m_c_cmip6 = np.zeros(1)
t850_m_cmip6 = np.zeros(1)
t850_c_cmip6 = np.zeros(1)
w850_m_cmip6 = np.zeros(1)
w850_c_cmip6 = np.zeros(1)
sic_m_cmip6 = np.zeros(1)
sic_c_cmip6 = np.zeros(1)
dt_m_cmip6 = np.zeros(1)
dt_c_cmip6 = np.zeros(1)
wr_m_cmip6 = np.zeros(1)
wr_c_cmip6 = np.zeros(1)
    
fit = np.polyfit(x_tot,t2m_cmip6,1, full=True)
t2m_m_cmip6 = fit [0][0]
t2m_c_cmip6 = fit [0][1]
fit = np.polyfit(x_tot,t850_cmip6,1, full=True)
t850_m_cmip6 = fit [0][0]
t850_c_cmip6 = fit [0][1]
fit = np.polyfit(x_tot,w10m_cmip6,1, full=True)
w10m_m_cmip6 = fit [0][0]
w10m_c_cmip6 = fit [0][1]
fit = np.polyfit(x_tot,w850_cmip6,1, full=True)
w850_m_cmip6 = fit [0][0]
w850_c_cmip6 = fit [0][1]
fit = np.polyfit(x_tot,sic_cmip6,1, full=True)
sic_m_cmip6 = fit [0][0]
sic_c_cmip6 = fit [0][1]
fit = np.polyfit(x_tot,dt_cmip6,1, full=True)
dt_m_cmip6 = fit [0][0]
dt_c_cmip6 = fit [0][1]
fit = np.polyfit(x_tot,wr_cmip6,1, full=True)
wr_m_cmip6 = fit [0][0]
wr_c_cmip6 = fit [0][1]

#Plot

color_hist = plt.cm.Blues(np.linspace(0,1,9))
color_fut = plt.cm.Oranges(np.linspace(0,1,10))
x_79 = np.arange(29,years_era5)
fig = plt.figure(figsize=(10,7))
ax = plt.subplot(4,2,1)
ax.plot(x, sic_era5, color='k', label='ERA5')
ax.plot(x, (x*sic_m_era5)+sic_c_era5, color='r', alpha=0.7)
ax.plot(x_tot, sic_cmip6, color=color_hist[5], label='CMIP6')
ax.plot(x_tot, (x_tot*sic_m_cmip6)+sic_c_cmip6, color='r', alpha=0.7)
ax.set_ylabel('SIC [%]', fontsize=10)
plt.xticks([0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150],[], fontsize=9)
ax.set_xlim([0,150])
ax.set_ylim([60,105])
plt.yticks([60,70,80,90,100],[60,70,80,90,100], fontsize=9)
plt.legend(loc='lower left')
plt.text(0.99, 0.99, 'a', ha='right', va='top', transform=plt.gca().transAxes, fontweight='bold')
ax = plt.subplot(4,2,2)
ax.plot(x, t2m_era5_global, color='k')
ax.plot(x, (x*m_globera5)+c_globera5, color='r', alpha=0.7)
ax.set_ylabel('Global average T$_{2m}$ [K]', fontsize=10)
plt.xticks([0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150],[], fontsize=9)
ax.set_xlim([0,150])
plt.text(0.99, 0.99, 'b', ha='right', va='top', transform=plt.gca().transAxes, fontweight='bold')
ax = plt.subplot(4,2,3)
ax.plot(x, t2m_era5, color='k')
ax.plot(x, (x*t2m_m_era5)+t2m_c_era5, color='r', alpha=0.7)
ax.plot(x_tot, t2m_cmip6, color=color_hist[5])
ax.plot(x_tot, (x_tot*t2m_m_cmip6)+t2m_c_cmip6, color='r', alpha=0.7)
ax.set_ylabel('T$_{2m}$ [K]', fontsize=10)
plt.xticks([0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150],[], fontsize=9)
ax.set_xlim([0,150])
plt.text(0.03, 0.99, 'c', ha='right', va='top', transform=plt.gca().transAxes, fontweight='bold')
ax = plt.subplot(4,2,4)
ax.plot(x, w10m_era5, color='k')
ax.plot(x, (x*w10m_m_era5)+w10m_c_era5, color='r', alpha=0.7)
ax.plot(x_79, (x_79*w10m_m_era5_79)+w10m_c_era5_79, linestyle='--', color='r', alpha=0.7)
ax.plot(x_tot, w10m_cmip6, color=color_hist[5])
ax.plot(x_tot, (x_tot*w10m_m_cmip6)+w10m_c_cmip6, color='r', alpha=0.7)
ax.set_ylabel('W$_{10m}$ [m/s]', fontsize=10)
plt.xticks([0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150],[], fontsize=9)
ax.set_xlim([0,150])
plt.text(0.99, 0.99, 'd', ha='right', va='top', transform=plt.gca().transAxes, fontweight='bold')
ax = plt.subplot(4,2,5)
ax.plot(x, t850_era5, color='k')
ax.plot(x, (x*t850_m_era5)+t850_c_era5, color='r', alpha=0.7)
ax.plot(x_tot, t850_cmip6, color=color_hist[5])
ax.plot(x_tot, (x_tot*t850_m_cmip6)+t850_c_cmip6, color='r', alpha=0.7)
ax.set_ylabel('T$_{850hPa}$ [K]', fontsize=10)
plt.xticks([0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150],[], fontsize=9)
ax.set_xlim([0,150])
plt.text(0.03, 0.99, 'e', ha='right', va='top', transform=plt.gca().transAxes, fontweight='bold')
ax = plt.subplot(4,2,6)
ax.plot(x, w850_era5, color='k')
ax.plot(x, (x*w850_m_era5)+w850_c_era5, color='r', alpha=0.7)
ax.plot(x_tot, w850_cmip6, color=color_hist[5])
ax.plot(x_tot, (x_tot*w850_m_cmip6)+w850_c_cmip6, color='r', alpha=0.7)
ax.set_ylabel('W$_{850hPa}$ [m/s]', fontsize=10)
plt.xticks([0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150],[], fontsize=9)
ax.set_xlim([0,150])
plt.text(0.99, 0.99, 'f', ha='right', va='top', transform=plt.gca().transAxes, fontweight='bold')
ax = plt.subplot(4,2,7)
ax.plot(x, dt_era5, color='k')
ax.plot(x, (x*dt_m_era5)+dt_c_era5, color='r', alpha=0.7)
ax.plot(x_tot, dt_cmip6, color=color_hist[5])
ax.plot(x_tot, (x_tot*dt_m_cmip6)+dt_c_cmip6, color='r', alpha=0.7)
ax.set_ylabel('T$_{850hPa}$-T$_{2m}$ [K]', fontsize=9)
plt.xticks([0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150],[1950,1960,1970,1980,1990,2000,2010,2020,2030,2040,2050,2060,2070,2080,2090,2100], rotation=45, fontsize=10)
ax.set_xlabel('Year')
ax.set_xlim([0,150])
plt.text(0.99, 0.99, 'g', ha='right', va='top', transform=plt.gca().transAxes, fontweight='bold')
ax = plt.subplot(4,2,8)
ax.plot(x, wr_era5, color='k')
ax.plot(x, (x*wr_m_era5)+wr_c_era5, color='r', alpha=0.7)
ax.plot(x_tot, wr_cmip6, color=color_hist[5])
ax.plot(x_tot, (x_tot*wr_m_cmip6)+wr_c_cmip6, color='r', alpha=0.7)
ax.set_ylabel('W$_{10m}$/W$_{850hPa}$', fontsize=9)
plt.xticks([0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150],[1950,1960,1970,1980,1990,2000,2010,2020,2030,2040,2050,2060,2070,2080,2090,2100], rotation=45, fontsize=10)
plt.text(0.03, 0.99, 'h', ha='right', va='top', transform=plt.gca().transAxes, fontweight='bold')
ax.set_xlabel('Year')
ax.set_xlim([0,150])
plt.tight_layout()
fig.savefig('Suppl_Figure1.pdf')