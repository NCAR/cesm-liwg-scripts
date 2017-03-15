#!/usr/bin/env python

import sys
import gzip
import glob
import os
import re
import numpy as np
import pandas as pd
from scipy.stats import linregress
import matplotlib.pyplot as plt
import deepOceanTemp
import netCDF4 as nc

import plotly.plotly as py              # for sending things to plotly
import plotly.graph_objs as go
import plotly.tools as tls              # for mpl, config, etc.
py.sign_in('JeremyFyke', 'u28du98cke')
plot_online=0 #if 0, then prints static plots locally.

plot_GLC_timeseries=1
plot_ocean_temperature_timeseries=0
plot_AMOC_timeseries=0
plot_calving_flux=1

suffix=''
suffix2='_restoring'
user_space='jfyke'

output_dir='plots'+suffix
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
def BuildFileList(GlobString):
    log_dir='/glade/scratch/'+user_space
    loglist=[]
    for i in np.arange(1,10):
	list_add=sorted(glob.glob(os.path.join(log_dir,'JG_iteration_'+str(i)+suffix+'/run/'+GlobString)),key=os.path.getmtime)
	loglist=loglist+list_add
	list_add=sorted(glob.glob(os.path.join(log_dir,'JG_iteration_'+str(i)+suffix2+'/run/'+GlobString)),key=os.path.getmtime)
	loglist=loglist+list_add	
	list_add=sorted(glob.glob(os.path.join(log_dir,'BG_iteration_'+str(i)+suffix+'/run/'+GlobString)),key=os.path.getmtime)
	loglist=loglist+list_add
    return loglist

def ExtractTSDataFromLog(TS,line,ValuePosition):
    #If the string associated with the Series 'name' attribute is matched...
    if re.search(re.escape(TS.name),line):
       #...then extract numerical value from string position...
       ExtractedValue=np.float(line.split()[ValuePosition])
       #...and append to time series object.
       TS=TS.set_value(len(TS),ExtractedValue)
    return TS

def TScolor(RunTypeColor,src_name):
    if re.search('BG',src_name):
	RunTypeColor.append(0)
    elif re.search('JG',src_name):
	RunTypeColor.append(1)
    else:
	print 'REGEX not working'
    return RunTypeColor

def runningMean(x, N):
    y = np.zeros((len(x),))
    for ctr in range(len(x)):
         y[ctr] = np.sum(x[ctr-6:(ctr-6+N)])
    return y/N

#General function to plot output time series.  This plot also contains
#a warning about zero-length time series, which is likely an indication
#of either bad paths/filenames, or a change to the formatting of CISM2
#output.
def WritePlot(Series,SeriesName,time,pointlabels,pointcolor):
    if Series.empty:
        print ''
        print 'Warning: '+SeriesName+' time series contains no data.'
        print SeriesName+'.png not written.'
    else:
        print 'Plotting: '+SeriesName+' timeseries.'
	RM=runningMean(Series,12)
	trace1=go.Scatter(x=np.array(time),y=np.array(Series),
	                  mode="markers",marker=dict(size='4',color=pointcolor,colorscale='Jet'),
			  name=SeriesName,text=pointlabels)
	trace2=go.Scatter( x=time[11:-12],y=RM[11:-12],name = SeriesName+'_smoothed',
	                  line=dict(color=('grey'),width=2) )
	figure=go.Figure(data=go.Data([trace1,trace2]),
	                 layout=go.Layout(hovermode='closest',title='', xaxis={'title':'time'}, yaxis={'title':SeriesName}))
	if plot_online:
	    py.iplot(figure, filename=SeriesName)  
	else:
	    py.image.save_as(figure,'plots'+suffix+'/'+SeriesName+'.png')
	
#First: define separate, empty Pandas Series for each time series.
#The 'name' attribute is used both as a regex to find the relevant
#output in the log files, and also to label the y-axis of the 
#resulting plots.
CISM2time=pd.Series(name='Diagnostic output, time',dtype=float)
Area=pd.Series(name='Total ice area (km^2)',dtype=float)
Volume=pd.Series(name='Total ice volume (km^3)',dtype=float)
Energy=pd.Series(name='Total ice energy (J)',dtype=float)
MeanThickness=pd.Series(name='Mean thickness (m)',dtype=float)
MeanTemperature=pd.Series(name='Mean temperature (C)',dtype=float)
MeanSMB=pd.Series(name='Mean accum/ablat (m/yr)',dtype=float)
MeanBasalMelt=pd.Series(name='Mean basal melt (m/yr)',dtype=float)
MaxThickness=pd.Series(name='Max thickness (m)',dtype=float)
MaxTemperature=pd.Series(name='Max temperature',dtype=float)
MinTemperature=pd.Series(name='Min temperature',dtype=float)
MaxSurfaceSpeed=pd.Series(name='Max sfc spd (m/yr)',dtype=float)
MaxBasalSpeed=pd.Series(name='Max base spd (m/yr)',dtype=float)
MeanOceanTemp=pd.Series(name='Ocean Temperature (C)',dtype=float)
CalvingFlux=pd.Series(name='Calving Flux (Sv)',dtype=float)
AMOC=pd.Series(name='AMOC (Sv)',dtype=float)

f = open('filelist.txt', 'wb')
f.write(str(BuildFileList('glc.log.*.gz')))

if plot_GLC_timeseries:
    RunTypeColor=[]
    FileNames=BuildFileList('glc.log.*.gz')
    for src_name in FileNames: #All tarred (successful) output files
	print 'Reading log file: '+src_name
	base=os.path.basename(src_name)
	with gzip.open(src_name,'rb') as f:
	  for line in f:
              #Set run type color for datapoint.  Used to plot JG/BG runs.
              if re.search(re.escape(CISM2time.name),line):
                  RunTypeColor=TScolor(RunTypeColor,src_name)
              CISM2time=ExtractTSDataFromLog(CISM2time,line,5)
              Area=ExtractTSDataFromLog(Area,line,5)
              Volume=ExtractTSDataFromLog(Volume,line,5)
              Energy=ExtractTSDataFromLog(Energy,line,5)
              MeanThickness=ExtractTSDataFromLog(MeanThickness,line,4)
              MeanTemperature=ExtractTSDataFromLog(MeanTemperature,line,4)
              MeanSMB=ExtractTSDataFromLog(MeanSMB,line,4)
              MeanBasalMelt=ExtractTSDataFromLog(MeanBasalMelt,line,5)
              MaxThickness=ExtractTSDataFromLog(MaxThickness,line,6)      
              MaxTemperature=ExtractTSDataFromLog(MaxTemperature,line,6)
              MinTemperature=ExtractTSDataFromLog(MinTemperature,line,6)  
              MaxSurfaceSpeed=ExtractTSDataFromLog(MaxSurfaceSpeed,line,7)
              MaxBasalSpeed=ExtractTSDataFromLog(MaxBasalSpeed,line,7)
    #Third: generate time series plots from log files
    WritePlot(Area,'Area',CISM2time,FileNames,RunTypeColor)
    WritePlot(Volume,'Volume',CISM2time,FileNames,RunTypeColor)
    WritePlot(Energy,'Energy',CISM2time,FileNames,RunTypeColor)
    WritePlot(MeanThickness,'MeanThickness',CISM2time,FileNames,RunTypeColor)
    WritePlot(MeanTemperature,'MeanTemperature',CISM2time,FileNames,RunTypeColor)
    WritePlot(MeanSMB,'MeanSMB',CISM2time,FileNames,RunTypeColor)
    WritePlot(MeanBasalMelt,'MeanBasalMelt',CISM2time,FileNames,RunTypeColor)
    WritePlot(MaxThickness,'MaxThickness',CISM2time,FileNames,RunTypeColor)
    WritePlot(MaxTemperature,'MaxTemperature',CISM2time,FileNames,RunTypeColor)
    WritePlot(MinTemperature,'MinTemperature',CISM2time,FileNames,RunTypeColor)
    WritePlot(MaxSurfaceSpeed,'MaxSurfaceSpeed',CISM2time,FileNames,RunTypeColor)
    WritePlot(MaxBasalSpeed,'MaxBasalSpeed',CISM2time,FileNames,RunTypeColor)

if plot_ocean_temperature_timeseries:
    RunTypeColor=[]
    FileNames=BuildFileList('*pop.h.0*')
    for src_name in FileNames: #All history files
	print 'Reading history file: '+src_name
	dd=deepOceanTemp.deepOceanTemp(src_name,Field='TEMP',minDepth=0.)
	MeanOceanTemp=MeanOceanTemp.set_value(len(MeanOceanTemp),dd.getIntegratedValue())
	RunTypeColor=TScolor(RunTypeColor,src_name)
    WritePlot(MeanOceanTemp,'Ocean_Temp',np.arange(0,len(MeanOceanTemp))/12.,FileNames,RunTypeColor)

if plot_AMOC_timeseries:    
    RunTypeColor=[]
    def getAMOC(src_name):
	vn=nc.Dataset(src_name)
	AMOC=np.amax(vn.variables['MOC'][0,1,0,:,:])
	print 'AMOC='+str(AMOC)
	return AMOC    
    
    RunTypeColor=[]
    FileNames=BuildFileList('*pop.h.0*')
    for src_name in FileNames: #All history files
	print 'Reading history file: '+src_name
	AMOC=AMOC.set_value(len(AMOC),getAMOC(src_name))
	RunTypeColor=TScolor(RunTypeColor,src_name)
    WritePlot(AMOC,'AMOC',np.arange(0,len(AMOC))/12.,FileNames,RunTypeColor)    
    
if plot_calving_flux:
    RunTypeColor=[]
    FileNames=BuildFileList('glc.log.*.gz')
    for src_name in BuildFileList('*cpl.hi.0*-12-*'): #All coupler history files
	print 'Reading history file: '+src_name
	vn=nc.Dataset(src_name)
	calving=np.sum(vn.variables['g2x_Fogg_rofi'])*196.9e6/(4.*3.141)/1000./1.e6 #convert from sum of steridians to Sv (m^6/s)
	CalvingFlux=CalvingFlux.set_value(len(CalvingFlux),calving)
	RunTypeColor=TScolor(RunTypeColor,src_name) 
    WritePlot(CalvingFlux,'Calving_flux',np.arange(0,len(CalvingFlux)),FileNames,RunTypeColor)





