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

import matplotlib.pyplot as plt
import matplotlib



class DataReaderCesmVector(object):
   def __init__(self, varname, filename, vector_info, ice_cover, fac):
      """
      Read column-output CESM variable from file and store for later use.
      
      varname:       variable name
      filename:      data file
      vector_info:   auxiliary data
      fac:           conversion factor that is immediately applied to the variable
      """

      print(filename)
      with Dataset(filename,'r') as fid:
         self.time = fid.variables['time'][:]
         self.time_units = fid.variables['time'].units
         self.data_cesm = fid.variables[varname][:]
         #print(fid.variables)
         self.var_type = fid.variables[varname].dimensions[-1] # 'col' or 'pft'
         self.long_name = fid.variables[varname].long_name
         self.units = fid.variables[varname].units
      
      print(self.time)
      print(self.time_units)

      # WORKAROUND shift data by one month
      #self.data_cesm = self.data_cesm[[1,2,3,4,5,6,7,8,9,10,11,0]]
      
      print(varname, self.var_type)
      
      self.data_cesm *= fac   
   
      # Hack for SNO_GS / layered data: use top layer only
      if (varname[0:4] == "SNO_"):
         self.data_cesm = self.data_cesm[:,0,:]
      
      self.ntime, self.nvec = self.data_cesm.shape
   
               
      """
      vector data descriptors
      """
      filename = vector_info
      with Dataset(filename,'r') as fid:
         self.lats = fid.variables['lat'][:]
         self.lons = fid.variables['lon'][:]           
         if (self.var_type == 'column'):
            self.ixy = fid.variables['cols1d_ixy'][:]
            self.jxy = fid.variables['cols1d_jxy'][:]
            self.lunit   = fid.variables['cols1d_itype_lunit'][:]   # col landunit type (vegetated,urban,lake,wetland,glacier or glacier_mec)
            self.coltype   = fid.variables['cols1d_itype_col'][:] 
         elif (self.var_type == 'pft'):
            self.ixy = fid.variables['pfts1d_ixy'][:]
            self.jxy = fid.variables['pfts1d_jxy'][:]
            self.lunit   = fid.variables['pfts1d_itype_lunit'][:]   # col landunit type (vegetated,urban,lake,wetland,glacier or glacier_mec)
            self.coltype   = fid.variables['pfts1d_itype_col'][:] 
         else:
            raise Exception('Unknown variable type: '+self.var_type)   
      
      
      self.nlats = len(self.lats)
      self.nlons = len(self.lons)
     
      self.setGlcFracCPL(ice_cover)


   def setGlcFracCPL(self,filename):
      """
      set glacier fraction from coupler history file

      variable: 
         x2lavg_Sg_ice_covered00 etc.
      """     
      self.mec_frac = np.zeros((GLC_NEC+1,self.nlats,self.nlons)) # One extra for tundra class
      with Dataset(filename,'r') as fid:
         #print(fid.variables)
         for i in range(0,11):
            self.mec_frac[i,:,:] = fid.variables['x2lavg_Sg_ice_covered%02d' % i][:]
      self.mec_frac = np.ma.masked_greater(self.mec_frac,2)


   def setGlcFracSurfdat(self, filename):
      """
      set glacier fraction from a surfdat file

      variable:
         double PCT_GLC_MEC_ICESHEET(nglcec, lsmlat, lsmlon)
      """
      self.mec_frac = np.zeros((GLC_NEC+1,self.nlats,self.nlons)) # One extra for tundra class
      with Dataset(filename,'r') as fid:
         #print(fid.variables)
         for i in range(1,11):
            self.mec_frac[i,:,:] = fid.variables['PCT_GLC_MEC_ICESHEET'][i-1,:,:]


   def getGridded3d(self):
      """
      Returns vector data as gridded (lat/lon) numpy array with levels (MEC).
      """
      lon2d,lat2d = np.meshgrid(self.lons,self.lats)
      var_out = np.ma.zeros((self.ntime,self.nlats,self.nlons,GLC_NEC))

      # Mask out all points without GLC_MEC
      var_out[:] = np.ma.masked
      
      for lev in range(GLC_NEC):
         mask = (self.coltype==(400+lev+1))
         idx, = np.where(mask)
         ix = self.ixy[idx]-1
         iy = self.jxy[idx]-1
         
         #print(tskin[:,idx].shape, var_out[:,iy,ix,lev].shape)
         var_out[:,iy,ix,lev] = self.data_cesm[:,idx]
   
      # Mask out points with missing value
      var_out = np.ma.masked_greater(var_out, 1e34)

      # report number of non-missing points
      print('number of non-zero points:', var_out.count() / self.ntime)
      return var_out
   

   def getGridded2d(self):
      """
      Returns vector data as gridded (lat/lon) numpy array.
      No levels, so weighted by ice_cover percentage.
      """
      var3d = self.getGridded3d() # dimensions (ntime, nlats, nlons, NGLC_NEC)

      #self.mec_frac = np.zeros((GLC_NEC+1,self.nlats,self.nlons)) # One extra for tundra class

      var_out = np.ma.zeros((self.ntime,self.nlats,self.nlons))

      frac = np.ma.filled(self.mec_frac, 0.0)
      for lev in range(GLC_NEC):
         #var_out += self.mec_frac[lev+1,:,:] * var3d[:,:,:,lev]
         var_out += frac[lev+1,:,:] * np.ma.filled(var3d[:,:,:,lev], 0.0)

      var_out /= np.sum(frac[1:,:,:], axis=0) # normalize for total fraction ( /= 1.0 when tundra present)

      # Mask out all points without GLC_MEC
      var_out = np.ma.masked_less(var_out, 1e-4)
      return var_out
           

   def get_profiles(self, imonth, lat_bounds, lon_bounds):
      """
      Return monthly profile integrated over an area which given by the bounds
      """

      var_out = np.zeros((self.ntime,self.nlats,self.nlons,GLC_NEC))    
      
      # :ctype_landice = 3 ;
      # :ctype_landice_multiple_elevation_classes = "4*100+m, m=1,glcnec" ;
      for lev in range(GLC_NEC):
         mask = (self.coltype==(400+lev+1))
         idx, = np.where(mask)
         ix = self.ixy[idx]-1
         iy = self.jxy[idx]-1
         
         #print(tskin[:,idx].shape, var_out[:,iy,ix,lev].shape)
         var_out[:,iy,ix,lev] = self.data_cesm[:,idx]
            
      # Mask out all points without GLC_MEC
      var_out = np.ma.masked_less(var_out, 1e-4)
   
      # Mask out points with missing value
      var_out = np.ma.masked_greater(var_out, 1e34)
   
   
      # some region around Greenland   
      r1=range(158,186)
      r2=range(220,288)
      ind = np.ix_(r1,r2)
          
      # ---------------------------------------------------
      # Plotting TG on a single level
      # ---------------------------------------------------
      lev = 0 # fifth ICE level (lev=4) ranges (1000,1300)m
   
      #levels = np.linspace(260,290,11)
      levels = np.linspace(265,275,11)
      cmap=plt.get_cmap()
      norm=matplotlib.colors.BoundaryNorm(levels,ncolors=cmap.N,clip=False)
      
      Z = var_out[imonth,:,:,lev]
      
      # Mask away virtual cells (mec_frac = 0.0)
      # ICE LEVELS START AT 1; ZERO IS BARE LAND
      #print(mec_frac.shape)
      Z = np.ma.masked_where(self.mec_frac[lev+1] < 1e-3, Z) 
      
      if (False):
         plt.figure()   
         plt.pcolormesh(lon2d[ind],lat2d[ind],Z[ind],cmap=cmap,norm=norm) #,levels=levels,extend='both')
         plt.title('TG lev = %d'%lev)
         plt.colorbar()
      
      """
      Elevation class heights
      
      shr/glc_elevclass_mod.F90
          case(10)
             topomax = [0._r8,   200._r8,   400._r8,   700._r8,  1000._r8,  1300._r8,  &
                                1600._r8,  2000._r8,  2500._r8,  3000._r8, 10000._r8]
             
      """
      topomax = np.asarray([0.,   200.,   400.,   700.,  1000.,  1300., 1600.,  2000.,  2500.,  3000., 10000])
      bin_center = (topomax[1::] + topomax[:-1])*0.5
                            
      b1 = np.logical_and(lon_bounds[0] <= lon2d, lon2d <= lon_bounds[1])
      b2 = np.logical_and(lat_bounds[0] <= lat2d, lat2d <= lat_bounds[1])
      B = np.logical_and(b1,b2)
      Bind = np.where(B) # indices
      Bmask = np.where(B,1,0) # mask (0,1)
      
      if (False):
         plt.figure()
         plt.pcolormesh(lon2d[ind],lat2d[ind],(Z*Bmask)[ind],edgecolor='k', lw=0.01,antialiased=False)
         plt.title('Western greenland mask lev = %d'%lev)
         plt.colorbar()   
         print('number of cells in region', Z[Bind].count())
      
      # CALCULATE MEAN TEMPERATURE AT EACH ELEVATION CLASS IN THIS REGION
      zmeans = np.zeros(GLC_NEC)
      zstd = np.zeros(GLC_NEC)
      for lev in range(GLC_NEC):
         Z = var_out[imonth,:,:,lev]
      
         Z = np.ma.masked_where(self.mec_frac[lev+1] < 1e-3, Z)
         zmeans[lev] = Z[Bind].mean()
         zstd[lev] = Z[Bind].std()
         #print('lev=',lev,'number of cells in region', Z[Bind].count())
   
      return zmeans, bin_center,zstd


   def divideByGriddedField(self,gfield):
      """ 
      To calculate albedo, we need to divide FSR by a regular lat/lon field. 
      This routine implements the division by a gridded field. 
      """
      print(self.nlats, self.nlons)
      print(gfield.shape)
      if(self.nlats != gfield.shape[0] or self.nlons != gfield.shape[1]):
         raise ValueError('grid dimensions do not match!')
      
      #print(self.ixy[0:10]-1)
      #print(self.jxy[0:10]-1)      


      #coords = list(zip(self.ixy-1 ,self.jxy-1))
      #print(coords)
      print(gfield[0,0])
      print(gfield[100,100])

      #print(gfield[self.jxy-1,self.ixy-1][1000:1010])
      fsds = gfield[self.jxy-1,self.ixy-1]
      print(self.data_cesm.shape)
      print(fsds.shape)
      self.data_cesm /= fsds
         #print(tskin[:,idx].shape, var_out[:,iy,ix,lev].shape)
         #var_out[:,iy,ix,lev] = self.data_cesm[:,idx]
