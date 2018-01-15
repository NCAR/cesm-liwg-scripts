## Description

Some scripts I wrote to convert CLM column or vector output to something gridded. The two most relevant high-level scripts are:

* vector2gridded2d.py  
uses PCT_GLC_MEC to calculate grid cell mean values ( this should give 2D gridded equal to the ICE_* variables)

* vector2gridded3d.py  
yields 3d (levelled by MEC) fields. Can be used further to interpolate to any surface. 

Both scripts share some underlying code:

* DataReaderCesmVector.py  
this code implements a 'datareader' class, that (1) stores the vector data in memory, and (2) has member functions that operate on that data. 
In order to single out a certain landunit types (e.g. GLC_MEC) , the following piece of code is essential: 

      for lev in range(GLC_NEC):
         mask = (self.coltype==(400+lev+1))
         idx, = np.where(mask)
         ix = self.ixy[idx]-1
         iy = self.jxy[idx]-1
         var_out[:,iy,ix,lev] = self.data_cesm[:,idx]

The hardcoded 400 signals GLC_MEC columns and may need to be changed depending on your application.


## Procedure for remapping CLM vector output to RACMO grid

1. convert vector to 3d gridded data: vector2gridded3d.py
2. horizontal interpolation: convert_grid.bsh OR cdo remapnn command below
3. vertical interpolation: interpolate_to_racmo_topo.py
4. (optional) merge all together  
`cdo -O merge *_downscaled.nc downscaled_merged.nc`

## Horizontal remapping example (step 2)
Example for Cheyenne
```
alias cdo='/glade/u/home/lvank/bin/cdo-1.8.2/src/cdo'
FILE=gridded3d_QICE_MELT_b.e20.BHIST.f09_g17.20thC.190_ramp204_reset.002_1961-1990_ymonmean.nc
NEWFILE=$(echo $FILE | sed 's/ymonmean/ymonmean_racmo/')
cdo remapnn,grid_ZGRN11_new_CDO_1.8.txt $FILE $NEW
```
Batch script is called `convert_grid.bsh`
