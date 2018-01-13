#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
20 Oct 2017

@author: Leo van Kampenhout
"""
import numpy as np
import xarray as xr

case = 'b.e20.BHIST.f09_g17.20thC.190_ramp204_reset.002'
#period = '1980-2005'
period = '1961-1990'
filetype="ymonmean_racmo"

topo_vec    = 'topo_racmo.nc'
#var_vec     = 'tsa_racmo.nc'

topo_target = '/glade/p/work/lvank/racmo/racmo23p2_GRN_monthly/elev.nc'
#topo_target = '/Users/leo/workspace/data/racmo/racmo23p2/elev.nc'

#varlist = 'EFLX_LH_TOT FGR FIRA FIRE FSA FSH FSM FSR QRUNOFF QSNOFRZ QSNOMELT QSOIL U10'.split()
#varlist += 'RAIN_REPARTITIONED SNOW_REPARTITIONED'.split()
#varlist += 'TSA TG'.split()
varlist = 'QICE'.split()
print(varlist)

for varname in varlist: 
   var_vec = "gridded3d_%s_%s_%s_%s.nc" % (varname,case,period,filetype)
   print(var_vec)

   d1 = xr.open_dataset(topo_vec)
   print(d1.dims)
   d2 = xr.open_dataset(var_vec)
   print(d2.dims)
   d3 = xr.open_dataset(topo_target)
   print(d3.dims)
   d3.coords['time'] = d2.coords['time']
   print(d3.dims)
   
   h = d1['TOPO_COL'].values[0,:,:,:]
   print(h.shape)
   
   v = d2[varname].values[:]
   print(v.shape)
   print(v.dtype)
   (ntime, nlev, nrlon, nrlat) = v.shape
   
   elev = d3['Elevation'].values
   
   # mask missing values
   hh=np.ma.masked_greater(h,1e30)
   vv=np.ma.masked_greater(v,1e30)
   
   td = np.ma.zeros((ntime,nrlon,nrlat)) # result
   td[:] = np.ma.masked
   print(varname,td.shape,td.dtype)
   
   import time
   t1 = time.time()
   
   # Note: this loop takes about 70 sec... Swapping the order of loop indices does not help...
   for itime in np.arange(ntime):
      for ilon in np.arange(nrlon):
         for ilat in np.arange(nrlat):
            if (not np.any(vv.mask[itime,:,ilon,ilat])):
               td[itime, ilon,ilat] = np.interp(elev[ilon,ilat],hh[:,ilon,ilat],vv[itime,:,ilon,ilat])
   
   t2 = time.time()
   print('elapsed [s]',t2-t1)
   
   d4 = xr.Dataset()
   d4[varname] = (('time','rlat','rlon'), td)
   print(d4)
   
   d4.to_netcdf(varname+'_downscaled.nc')
   

