#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 11:06:44 2017

@author: leo
"""

GLC_NEC = 10 # maximum number of elevation classes present in input file
COLUNIT_GLCMEC = 4 # for landunit types and column types, land ice = 7 (older CLM) land ice = 4 (newer CLM)

import numpy as np
from netCDF4 import Dataset, default_fillvals
import matplotlib
import matplotlib.pyplot as plt
import os
import os.path

from string import Template
from DataReaderCesmVector import DataReaderCesmVector



# =========================
# 2-D Grid information
# =========================
ice_cover = '../cesm2beta_run144/ice_cover144.nc'
vector_info = 'vector_info183.nc'

with Dataset(vector_info,'r')  as fid:
   nlon = len(fid.dimensions['lon'])
   print('Longitude = '+str(nlon))
   nlat = len(fid.dimensions['lat'])
   print('Latitude = '+str(nlat))
   lats_      = fid.variables['lat'][:]               # latitude array
   lons_      = fid.variables['lon'][:]               # longitude array


#filetype="ymonmean"
filetype="yearmean"

S = Template('/glade/scratch/lvank/archive/avg2/${case}/${period}/vector/${varname}_${case}_'+filetype+'.nc')

cases = []
#cases += [('b.e20.BHIST.f09_g17.20thC.183_01_cpl3h.003','1870','1873')]
#cases += [('b.e20.BHIST.f09_g17.20thC.183_01_gs300_nmelt10','1870','1887')]
#cases += [('b.e20.BHIST.f09_g17.20thC.183_01_gs204_nmelt10','1870','1878')]
#cases += [('b.e20.BHIST.f09_g17.20thC.183_01_gs300_nmelt10_reset','1873','1878')]
#cases += [('b.e20.BHIST.f09_g17.20thC.183_01_ramp54_204_nmelt10_reset','1870','1875')]
cases += [('b.e20.BHIST.f09_g17.20thC.183_01_gs300_nmelt10_rfgrain_iwc5_glc1m','1870','1884')]
cases += [('b.e20.BHIST.f09_g17.20thC.183_01_gs300_nmelt10_reset','1893','1897')]

xlist = []
#xlist+=[('SNO_GS','microns','grain size',1)]
xlist+=[('H2OSNO','m','snow depth',1)]

ice_cover_including_ANT = 'pct_glc_mec_icesheet.nc'

for case,ys,ye in cases:
   period  = ys +'-' + ye

   for x in xlist:
      print(case,x)

      varname, units, vardesc, fac = x
      filename = S.substitute(varname=varname, case=case, period=period)
   
      datareader1 = DataReaderCesmVector(varname, filename, vector_info, ice_cover, fac)
      datareader1.setGlcFracSurfdat(ice_cover_including_ANT) # Antarctica missing in coupler ice frac

      sno_gs = datareader1.getGridded2d()
   
      # =====================================
      # Write NetCDF file with netcdf4-python
      # =====================================
      name = "gridded2d_%s_%s_%s_%s.nc" % (varname,case,period,filetype)
   
      # Open a new NetCDF file to write the data to. For format, you can choose from
      # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
      ncfile = Dataset(name, 'w', format='NETCDF4')
      ncfile.description = 'Gridded field from CLM vector output'
   
      # Create dimensions
      ncfile.createDimension('longitude', nlon)
      ncfile.createDimension('latitude', nlat)
      ncfile.createDimension('time', None)
   
      # Define the coordinate var
      lons   = ncfile.createVariable('longitude', 'f4', ('longitude',))
      lats   = ncfile.createVariable('latitude', 'f4', ('latitude',))
      times    = ncfile.createVariable('time', 'f8', ('time',))
   
      # Assign units attributes to coordinate var data
      lons.units   = "degrees_east"
      lats.units   = "degrees_north"
      #times.units    = "days since 1-01-01 00:00:00"
      times.units = datareader1.time_units
   
      #levs.units   = "MEC level number"
   
      # Write data to coordinate var
      lons[:]    = lons_
      lats[:]    = lats_
      #times[:]   = times_
      times[:] = datareader1.time
      #levs[:]    = range(0,GLC_NEC)
   
   	# ----------
   	# WRITE DATA
   	# ----------
   
      print(sno_gs.shape) #(12, 192, 288, 10)
      #sno_gs = sno_gs.transpose((0,1,2)) # permute columns
      #print(sno_gs.shape) #(12, 192, 288, 10)
   
      # Create output variable of correct dimensions
      var = ncfile.createVariable(varname,'f4',('time','latitude','longitude',),fill_value=default_fillvals['f4'])
      var.units = "microns"
      var.long_name = vardesc 

#      var[:,:,:] = default_fillvals['f4'] # Initialise with missing value everywhere (will be replaced later)
   	
      var[:] = sno_gs[:]
   
   #   for imonth in range(0,12):
      ncfile.close()
   
   

