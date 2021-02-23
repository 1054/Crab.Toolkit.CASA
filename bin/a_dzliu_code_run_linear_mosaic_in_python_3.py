#!/usr/bin/env python
# 
# Run this code with python3
# 
# Note that the beams are different among images!
# 

import os, sys, re, glob, shutil, copy
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.coordinates import SkyCoord, FK5
from regions import DS9Parser, read_ds9, write_ds9, RectangleSkyRegion, CircleSkyRegion, PixelRegion, LinePixelRegion, ds9_objects_to_string
from reproject import reproject_interp
from scipy.interpolate import griddata
sys.path.append('/Users/dzliu/Cloud/Github/Crab.Toolkit.CASA/lib/python')
from dzliu_linear_mosaic_in_python import dzliu_linear_mosaic


# Find all images
Region_files = glob.glob('ds9_regions_of_mosaic_pointings_in_*_boxes_unpadded.reg')
DataSet_names = [re.sub(r'ds9_regions_of_mosaic_pointings_in_(.*)_boxes_unpadded.reg', r'\1', os.path.basename(t)) for t in Region_files] # [os.path.basename(os.path.dirname(os.path.dirname((t)))) for t in DataSet_dirs]


# Check input
if len(Region_files) == 0:
    print('Error! File "ds9_regions_of_mosaic_pointings_in_*_boxes_unpadded.reg" not found! Please run \'a_dzliu_code_print_phasecenters_in_casa.py\' or \'a_dzliu_code_run_dividing_mosaic_and_imaging_in_casa.py\' first!')
    sys.exit(255)


# Start processing
print('getcwd: %s'%(os.getcwd()))


# loop datasets
for Region_file, DataSet_name in list(zip(Region_files, DataSet_names)):
    
    print('DataSet: '+DataSet_name)

    # read ds9 regions
    list_of_regions = read_ds9(Region_file)
    
    # find images
    list_of_images = glob.glob('Level_4_Data_Images_Divide_Mosaic/%s_Mosaic_*/output_*_clean.image.pbcor.fits'%(DataSet_name))
    list_of_images = sorted(list_of_images)
    
    # run 
    output_mosaic_image = 'Level_4_Data_Images_Divide_Mosaic/%s_Mosaicked.fits'%(DataSet_name)
    dzliu_linear_mosaic(list_of_images, output_mosaic_image)



# done
print('Done!')








