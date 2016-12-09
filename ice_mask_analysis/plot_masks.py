#!/usr/bin/env python
"""
   DESCRIPTION
      Plot different ice masks
      Produces both a multipage PDF and PNGs of each page

   AUTHOR
      Leo van Kampenhout

   DATE
      09-DEC-2016
"""

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.basemap import Basemap

import numpy as np
import sys,os
from netCDF4 import Dataset, default_fillvals

A4_size=(8.27, 11.69) # size in inches

def main():
   pp = PdfPages('masks.pdf')
  
   varname='ICE_MASK'

   filename = 'cesm_4km/CESM_ICE_1DEG.nc'
   with Dataset(filename, mode='r') as fid:
      lons = fid.variables['lon'][:]
      lats = fid.variables['lat'][:]
      mask4km = fid.variables['ICE_MASK'][:,:]

   filename = 'cesm_5km/CESM_ICE_1DEG.nc'
   with Dataset(filename, mode='r') as fid:
      mask5km = fid.variables['ICE_MASK'][:,:]

   filename = 'racmo23/RACMO23_ICE_1DEG.nc'
   with Dataset(filename, mode='r') as fid:
      mask_rac23 = fid.variables['GrIS_mask'][:,:] * 100.

   filename = 'racmo24/RACMO24_ICE_1DEG.nc'
   with Dataset(filename, mode='r') as fid:
      mask_rac24 = fid.variables['Promicemask'][:,:] * 100.

   lons2d, lats2d = np.meshgrid(lons, lats)

   #fig = plt.figure(figsize=A4_size)
   fig = plt.figure(figsize=A4_size)
   #fig, axarr = plt.subplots(2, 2)
   
   clevs=np.arange(0,110,10)
   #cs = m.contourf(x,y,mask4km,clevs,cmap='magma')

   ax = plt.subplot(221)
   m = setup_map_greenland(ax)
   x, y = m(lons2d, lats2d) # compute map proj coordinates.
   #cs = m.pcolormesh(x,y,mask4km,cmap='magma_r')
   cs = m.contourf(x,y,mask4km,cmap='magma_r',levels=clevs)
   ax.set_title('CESM 2.0 (off 4 km)')
   cbar = m.colorbar(cs,location='right',pad="5%")
   cbar.set_label("Percent ice sheet")

   ax = plt.subplot(222)
   m = setup_map_greenland(ax)
   cs = m.contourf(x,y,mask5km,cmap='magma_r',levels=clevs)
   ax.set_title('CESM 1.6 (off 5 km)')
   cbar = m.colorbar(cs,location='right',pad="5%")
   cbar.set_label("Percent ice sheet")

   ax = plt.subplot(223)
   m = setup_map_greenland(ax)
   cs = m.contourf(x,y,mask_rac23,cmap='magma_r',levels=clevs)
   ax.set_title('RACMO 2.3')
   cbar = m.colorbar(cs,location='right',pad="5%")
   cbar.set_label("Percent ice sheet")

   ax = plt.subplot(224)
   m = setup_map_greenland(ax)
   cs = m.contourf(x,y,mask_rac24,cmap='magma_r',levels=clevs)
   ax.set_title('RACMO 2.4')
   cbar = m.colorbar(cs,location='right',pad="5%")
   cbar.set_label("Percent ice sheet")

   plt.tight_layout()
   pp.savefig(fig,dpi=300)
   plt.savefig('all_masks_1deg.png',dpi=300)
   
   fig = plt.figure(figsize=A4_size)
   ax = plt.subplot(111)
   m = setup_map_greenland(ax)
   clevs=np.arange(-100,110,10)
   #cs = m.pcolormesh(x,y,mask5km-mask_rac23,cmap='bwr',levels=clevs)
   cs = m.contourf(x,y,mask5km-mask_rac23,cmap='bwr',levels=clevs)
   ax.set_title('CESM 1.6 - RACMO2.3')
   cbar = m.colorbar(cs,location='right',pad="5%")
   cbar.set_label("Percent ice sheet")

   pp.savefig(fig,dpi=300)
   #plt.savefig('masks.pdf',dpi=300)
   plt.savefig('cesm16_sub_racmo23.png',dpi=300)
   pp.close()


def setup_map_greenland(ax):
   """
   Attach basemap for Greenland map plotting
   """
   m = Basemap(projection='stere',lon_0=-40,lat_0=60., \
            llcrnrlat=59.,urcrnrlat=82.,\
            llcrnrlon=-55,urcrnrlon=10,\
            rsphere=6371200.,resolution='l',area_thresh=10000,ax=ax)
   m.drawcoastlines()
   return m


if(__name__)=="__main__":
   main()
