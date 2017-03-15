#!/usr/bin/env python

from netCDF4 import Dataset
from numpy import *

from plotObj import *

from MidpointNormalize import *

import matplotlib.pyplot as plt
import matplotlib.colors as colors

import glob

from progress import ProgressMeter

###########################


class dataServer:
    def __init__(self,Files):

        self.Files = Files

        try:
            self.lat = Dataset(self.Files[0],'r').variables['lat_aux_grid'][:]
            self.z   = Dataset(self.Files[0],'r').variables['moc_z'][:]
        except: 
            pass

        self.iz = self.findZidx(500.)

        latS = 28.
        latN = 80.

        if latS != None:
            self.latSidx = self.findLatIdx(latS)
        else: self.latSidx = 0

        if latN != None:
            self.latNidx = self.findLatIdx(latN)
        else: self.latNidx = len(self.lat)


    def extractTimeSeries(self,Field='MOC',latS=None,latN=None):

        if latS != None: sIdx = self.findLatIdx(latS)
        else: sIdx = self.latSidx 

        if latN != None: nIdx = self.findLatIdx(latN)            
        else: nIdx = self.latNidx

        ii = 0
        mocMax = 0.
        for File in self.Files:
            moc = amax(Dataset(File,'r').variables['MOC'][0,1,0,self.iz:,sIdx:nIdx+1])
            mocMax += moc
            ii += 1 
        return mocMax/ii


    def extractTimeSeriesFromSingleFile(self,Field='MOC'):
        sIdx = self.latSidx
        nIdx = self.latNidx

        dd = Dataset(self.Files,'r')
        time = dd.variables['time']
        tt = linspace(0,len(time)-12,len(time)/12)

        var = []
        for it in tt:
            it = int(it)
            xx = 0.
            for im in range(0,12):
                xx += amax(dd.variables['MOC'][it+im,1,0,self.iz:,sIdx:nIdx+1])/12.
            var.append(xx)
        return var


    def findLatIdx(self,ilat):
        return abs(self.lat-ilat).argmin()

    def findZidx(self,iz):
        return abs(self.z-iz*100.).argmin()


###########################

def catTimeSeries(x1,x2):
    ## Concatenate two time series
    return concatenate(x1,x2)


###########################

def plotVar(var,time):
    pp = plotObj()
    pp.linewidths = 1.5
    pp.color = '-k'
    pp.linePlot(var,time)
    pp.showFigure()



def plotVar2(var,time,start,end):

    pp = plotObj()
    pp.linewidths = 1.5

    idx = cumsum(end)

    ii = 0
    for it in range(len(time)):
        
        if it < idx[ii]:
            pp.color = '-g'
            pp.linePlot(var,time)

    pp.showFigure()




def getVarAndTimeJG(sims,start,end):
    latS = 28.
    latN = 80.

#    latS = -34.
#    latN = -34.

    var = []
    time = []
    ii = 0
    for sim in sims:
        dataPath = '/glade/scratch/jfyke/%s/run'%(sim)

        meter = ProgressMeter(total=end[ii]+1-start[ii])
        for yr in range(start[ii],end[ii]+1):
            year = '%04i'%yr

            Files = sorted(glob.glob('%s/%s.pop.h.%s-??.nc'%(dataPath,sim,year)))

            dd = dataServer(Files)
            var.append(dd.extractTimeSeries(latS=latS,latN=latN))
            meter.update(1)
            
        ii += 1
        

    var = squeeze(ma.array(var).filled().astype("f"))    
#    var = ma.array(var).filled().astype("f")
    time = arange(len(var))

    print var

    return var,time




###########################

def JG():
    sims = [
#        'JG_iteration_4',
#        'BG_iteration_4',
#        'JG_iteration_5',
#        'BG_iteration_5',
#        'JG_iteration_6',
#        'BG_iteration_6',
#        'JG_iteration_7',
        'BG_iteration_7',
        ]

    start = [
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        ]


    end = [
#        21,
#        40,
#        34,
#        40,
#        44,
#        40,
#        150,
        20,
           ]

    var,time = getVarAndTimeJG(sims,start,end)
    plotVar(var,time)
#    plotVar2(var,time,start,end)


    """
    pp = plotObj()
    pp.linewidths = 1.5

    idx = cumsum(end)

    time = arange(sum(end))
    var = zeros((len(time)))


    print idx
#    print len(time)
#    print time

    ii = 0
    for it in range(len(time)):

        if it <= idx[ii]:

            pass
            
#            print it,pp.color

        else:

            print 'here'

            mask = zeros(len(time),dtype=int)
            mask[idx[ii]:idx[ii+1]] = 1

            print idx[ii], idx[ii+1]

            print mask


            if ii % 2 == 0: 
                pp.color = '-g'

                var = ma.masked_array(var,mask=mask)
                pp.linePlot(var,time)

            else: 
                pp.color = '-r'

                var = ma.masked_array(var,mask=mask)
                pp.linePlot(var,time)

#            pp.linePlot(var[idx[ii]:idx[ii+1]],time)


            ii += 1

        
    


    pp.showFigure()
    """

###########################

if __name__ == "__main__":

    JG()

###########################
## === end of script === ##
###########################

