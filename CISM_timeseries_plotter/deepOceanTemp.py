#!/usr/bin/env python
#
# Compute depth and area weighted ocean temperature
#
# Author: 
# M Lofverstrom
# NCAR, Aug, 2016
#
###########################

import time
from netCDF4 import *
from numpy import *

###########################

class deepOceanTemp(object):
    def __init__(self,dataFile,Field='TEMP',minDepth=None):

        ## Dataset
        self.data = Dataset(dataFile,'r')
        self.var  = self.data.variables[Field][0,:,:,:]

        self.maxIdx = self.data.variables['KMT'][:]
        self.area   = self.data.variables['TAREA'][:]
        self.mask   = self.data.variables['REGION_MASK'][:]
        self.zt     = self.data.variables['z_t'][:]
        self.dz     = self.data.variables['dz'][:]
        if minDepth != None: self.minDepth = minDepth
        else: self.minDepth = 500.*100. # 500 m in cm
        
        ###

        """
        self.lat = self.data.variables['TLAT'][:]
        self.lon = self.data.variables['TLONG'][:]

        self.latMax = 30.
        self.latMin = -30.
        
        self.lonMax = 360.
        self.lonMin = 0.
        """

        self.sanityBound = 1.e2 # discard values greater than this number (fillvalues)



    def minZidx(self,depth=None):
        if depth == None: depth = self.minDepth
        return abs(self.zt-depth).argmin()


    def getIntegratedValue(self,depth=None):
        
        idxMin = self.minZidx(depth=depth) ## first index in vertical integration
        idxMax = shape(self.var)[0]        ## last  index in vertical integration

        var = zeros((shape(self.area)))*0. ## flattened variable (z collapse)

        ## Calculate area where variable is defined at start index
        totArea = sum(sum(where(self.var[idxMin,:,:]<self.sanityBound,\
                                    self.area,0.),axis=-1),axis=-1)

        ## Sum over vertical indices in column
        z = zeros((shape(self.area)))*0.
        for iz in range(idxMin,idxMax):
            z[:,:] += where(self.var[iz,:,:]<self.sanityBound,self.dz[iz,None,None],1.)
        
            var[:,:] += where(self.var[iz,:,:]<self.sanityBound,self.var[iz,:,:]* \
                                  self.dz[iz,None,None],0.)


        ## Divide by column depth from start index
        var = ma.masked_array(var/z)

        ## Sum over weighted area elements
        xx = nansum(nansum(var*self.area/totArea,axis=-1),axis=-1)

        print 'Integrated temperature:', xx

        return xx


###########################

if __name__ == "__main__":

    sim = 'bg.c15a06fv1_PI_synch_pseudoPlastSlid_cont36'
    dataDir = '/glade/u/home/marcusl/archive/%s/ocn/hist'%(sim)

    startYear = 120
    endYear   = 120


    for year in range(startYear,endYear+1):
        yr = '%04i'%year
        
        for month in range(1,13):
#        for month in range(1,2):
            mon = '%02i'%month

            dataFile = '%s/%s.pop.h.%s-%s.nc'%(dataDir,sim,yr,mon)

            dd = deepOceanTemp(dataFile)

            t = time.time()
            xx = dd.getIntegratedValue(depth=0.*1000.*100.)
            print 'Elapsed time:',time.time()-t


###########################
## === end of script === ##
###########################



