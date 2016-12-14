;
;  DESCRIPTION
;     Calculates the Greenland integrated sum of rain and snow using different ice masks. 
;
;  AUTHOR
;     Leo van Kampenhout
;
;  DATE
;     09-DEC-2016
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

;fid_rain = addfile("RAIN_REPARTITIONED_b.e15.B1850.f09_g16.pi_control.all.119_timmean.nc","r")
;fid_snow = addfile("SNOW_REPARTITIONED_b.e15.B1850.f09_g16.pi_control.all.119_timmean.nc","r")
;rain = fid_rain->RAIN_REPARTITIONED(0,:,:)
;snow = fid_snow->SNOW_REPARTITIONED(0,:,:)

fid_rain = addfile("RAIN_b.e15.B1850.f09_g16.pi_control.all.119_timmean.nc","r")
fid_snow = addfile("SNOW_b.e15.B1850.f09_g16.pi_control.all.119_timmean.nc","r")
rain = fid_rain->RAIN(0,:,:)
snow = fid_snow->SNOW(0,:,:)


fid_cesm16 = addfile("../cesm_5km/CESM_ICE_1DEG.nc","r")
mask_cesm16 = fid_cesm16->ICE_MASK(:,:)/100.

fid_cesm20 = addfile("../cesm_4km/CESM_ICE_1DEG.nc","r")
mask_cesm20 = fid_cesm20->ICE_MASK(:,:)/100.

;fid_racmo23 = addfile("../racmo23/RACMO23_ICE_1DEG.nc","r") 
;fid_racmo23 = addfile("../racmo23/gismask_RACMO23_0.9x1.25_bil.nc","r")
fid_racmo23 = addfile("../racmo23/gismask_RACMO23_0.9x1.25_con.nc","r")
mask_rac23 = fid_racmo23->GrIS_mask(:,:) 

;fid_racmo24 = addfile("../racmo24/RACMO24_ICE_1DEG.nc","r") 
fid_racmo24 = addfile("../racmo24/gismask_RACMO24_0.9x1.25_con.nc","r")
mask_rac24 = fid_racmo24->Promicemask(:,:) 

fid_area = addfile("area.nc","r")
area = fid_area->area

; factor to convert mm/s to Gt/yr
; account for area multiplication km^2
fac = 86400 * 365 / 1.e12 * 1.e6

do i = 1, 4, 1
   print("--------------")
   if (i .eq. 1) then
      print("CESM 1.6 ice mask")
      my_mask = mask_cesm16
   else if (i .eq. 2) then
      print("CESM 2.0 ice mask ")
      my_mask = mask_cesm20
   else if (i .eq. 3) then
      print("RACMO 2.3 ice mask")
      my_mask = mask_rac23
   else if (i .eq. 4) then
      print("RACMO 2.4 ice mask")
      my_mask = mask_rac24
   else
      print("ERROR: index not reached")
   end if
   end if
   end if
   end if
   gis_mul_area = my_mask*area
   gis_area=sum(gis_mul_area)
   print("GIS area = "+sprintf("%8.2f",gis_area)+" km^2")
   gis_rain=sum(rain*gis_mul_area)*fac
   gis_snow=sum(snow*gis_mul_area)*fac
   print("snow = "+gis_snow + " Gt")
   print("rain = "+gis_rain + " Gt")
   print("total = "+(gis_rain+gis_snow)+" Gt")
end do

end
