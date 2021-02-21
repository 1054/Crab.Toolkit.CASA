#
# Run this code inside CASA as:
#     execfile('')
# or
#     exec(open('').read())
#

import os, sys, re, glob
sys.path.append(os.path.expanduser('~/Cloud/Github/Crab.Toolkit.CASA/lib/python'))
from dzliu_clean_utils import get_synbeam_and_imcell, get_mosaic_imsize_and_phasecenter


# User-defined parameters
DataSet_dirs = ['Level_3_Split/DataSet_Merged_A/calibrated.ms', 'Level_3_Split/DataSet_Merged_B/calibrated.ms']
DataSet_names = ['DataSet_Merged_A', 'DataSet_Merged_B']
Target_names = ['Tune147', 'Tune139']


# Start processing
print('getcwd: %s'%(os.getcwd()))

for DataSet_dir, DataSet_name, Target_name in list(zip(DataSet_dirs, DataSet_names, Target_names)):
    
    vis = DataSet_dir
    field = Target_name
    output_ds9_region_file = 'ds9_regions_of_mosaic_pointings_in_'+DataSet_name+'.reg'
    
    synbeam, imcell = get_synbeam_and_imcell(vis)
    
    imcell = '0.20arcsec'
    
    imsize, phasecenter = get_mosaic_imsize_and_phasecenter(vis, imcell, galaxy_name=field, output_ds9_region_file=output_ds9_region_file)
    
print('Done!')







