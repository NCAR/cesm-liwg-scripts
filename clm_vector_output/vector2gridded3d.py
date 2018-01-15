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
ice_cover = '/Users/leo/workspace/data/cesm/ice_cover144.nc'
vector_info = '/Users/leo/workspace/data/cesm/vector_info183.nc'

with Dataset(vector_info,'r')  as fid:
   nlon = len(fid.dimensions['lon'])
   print('Longitude = '+str(nlon))
   nlat = len(fid.dimensions['lat'])
   print('Latitude = '+str(nlat))
   lats_      = fid.variables['lat'][:]               # latitude array
   lons_      = fid.variables['lon'][:]               # longitude array


#filetype="ymonmean"
filetype="ymonmean"

#S = Template('/glade/scratch/lvank/archive/avg2/${case}/${period}/vector/${varname}_${case}_'+filetype+'.nc')
S = Template('/Users/leo/workspace/data/cesm/${case}/${period}/vector/${varname}_${case}_'+filetype+'.nc')

cases = []
#cases += [('b.e20.BHIST.f09_g17.20thC.183_01_gs300_nmelt10_reset','1893','1897')]
#cases += [('b.e20.BHIST.f09_g17.20thC.183_01_gs300_nmelt10_reset','1893','1897')]
cases += [('b.e20.BHIST.f09_g17.20thC.190_ramp204_reset.002','1980','2005')]


#varlist = ['TSA']
varlist = ['RAIN_REPARTITIONED']
#varlist = 'EFLX_LH_TOT FGR FIRA FIRE FSA FSH FSM FSR QRUNOFF QSNOFRZ QSNOMELT QSOIL RAIN_REPARTITIONED SNOW_REPARTITIONED TG U10'.split()

#! Stream 2: monthly, vector
#hist_fincl2 = 'QRUNOFF', 'QSOIL', 'QSNOMELT', 'QSNOFRZ', 'H2OSNO', 'FSR', 'FSA', 'FIRE', 'FIRA', 'EFLX_LH_TOT', 'FSH', 'FGR', 'FSM', 'TG ', 'TSA', 'U10', 'TSOI', 'H2OSNO', 'SNOW_DEPTH', 'SNOWDP', 'SNOWLIQ', 'SNOWICE', 'FSNO', 'FSNO_EFF', 'SNOTTOPL', 'QICE_OLD', 'QICE', 'QICE_FRZ', 'QICE_MELT', 'SNO_TK', 'SNO_ABS', 'SNO_Z', 'SNO_T', 'SNO_LIQH2O','SNO_ICE', 'SNO_BW','SNO_EXISTENCE','SNO_GS', 'SNO_MELT', 'SNO_FRZ', 'RAIN_REPARTITIONED', 'SNOW_REPARTITIONED'

for case,ys,ye in cases:
   period  = ys +'-' + ye

   for varname in varlist:
      print(case,varname)

      filename = S.substitute(varname=varname, case=case, period=period)
      #filename = '/Users/leo/workspace/data/cesm/b.e20.BHIST.f09_g17.20thC.190_ramp204_reset.002/1980-2005/vector/rain_scaled.nc'

      # TODO: remove this
      #filename = '/glade/p/work/lvank/vector/topo_col_bg190.nc'
      #varname = 'TOPO_COL'
   
      datareader1 = DataReaderCesmVector(varname, filename, vector_info, ice_cover, 1)
      print('number of vectors = ',datareader1.nvec)
      var3d = datareader1.getGridded3d()
   
      # =====================================
      # Write NetCDF file with netcdf4-python
      # =====================================
      name = "gridded3d_%s_%s_%s_%s.nc" % (varname,case,period,filetype)
   
      # Open a new NetCDF file to write the data to. For format, you can choose from
      # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
      ncfile = Dataset(name, 'w', format='NETCDF4')
      ncfile.description = 'Gridded field from CLM vector output'
   
      # Create dimensions
      ncfile.createDimension('longitude', nlon)
      ncfile.createDimension('latitude', nlat)
      ncfile.createDimension('time', None)
      ncfile.createDimension('lev',GLC_NEC)
   
      # Define the coordinate var
      lons   = ncfile.createVariable('longitude', 'f4', ('longitude',))
      lats   = ncfile.createVariable('latitude', 'f4', ('latitude',))
      times    = ncfile.createVariable('time', 'f8', ('time',))
      levs   = ncfile.createVariable('lev', 'f4', ('lev',))
   
      # Assign units attributes to coordinate var data
      lons.units   = "degrees_east"
      lons.axis = "Y"
      lats.units   = "degrees_north"
      lats.axis = "X"
      #times.units    = "days since 1-01-01 00:00:00"
      times.units = datareader1.time_units
   
      levs.units   = "MEC level number"
   
      # Write data to coordinate var
      lons[:]    = lons_
      lats[:]    = lats_
      #times[:]   = times_
      times[:] = datareader1.time
      levs[:]    = range(0,GLC_NEC)
   
   	# ----------
   	# WRITE DATA
   	# ----------
      print(var3d.shape) #(12, 192, 288, 10)
      var3d = var3d.transpose((0,3,1,2)) # permute columns
      print(var3d.shape) #(12, 10, 192, 288)
   
      # Create output variable of correct dimensions
      var            = ncfile.createVariable(varname,'f8',('time','lev','latitude','longitude',))
      var.units      = datareader1.units
      var.long_name  = datareader1.long_name
      var[:,:,:,:]   = default_fillvals['f8'] # Initialise with missing value everywhere (will be replaced later)
   	
      var[:] = var3d[:]
   
   #   for imonth in range(0,12):
      ncfile.close()
   
   

