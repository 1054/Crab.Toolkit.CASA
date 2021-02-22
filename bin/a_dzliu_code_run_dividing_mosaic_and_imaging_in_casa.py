# 
# Run this code inside CASA as:
#     execfile('')
# or
#     exec(open('').read())
# 

import os, sys, re, glob
sys.path.append(os.path.expanduser('~/Cloud/Github/Crab.Toolkit.CASA/lib/python'))
from dzliu_clean_utils import (get_datacolumn, get_synbeam_and_imcell, get_mosaic_imsize_and_phasecenter, get_field_IDs_in_mosaic, get_field_names, print_params, 
                               cleanup_tclean_products, apply_pbcor_to_tclean_image, export_tclean_products_as_fits_files)
from taskinit import casalog, tb


# User-defined parameters
DataSet_dirs = ['Level_3_Split/DataSet_Merged_A/calibrated.ms'] # sorted(glob.glob('Level_3_Split/DataSet_01/split_GOODS-S_sci.ms'))
DataSet_names = ['DataSet_Merged_A'] # [os.path.basename(os.path.dirname(os.path.dirname((t)))) for t in DataSet_dirs]
Target_names = ['Tune147']
ncol = 5
nrow = 10
imcell = '0.25arcsec'
dry_run = False


# def print2
def print2(message):
    print(message)
    casalog.post(message, 'INFO')


# Start processing
print2('getcwd: %s'%(os.getcwd()))

clean_params = {}
clean_params['selectdata'] = True
clean_params['spw'] = ''
clean_params['gridder'] = 'mosaic' # 'standard'
clean_params['specmode'] = 'mfs'
clean_params['outframe'] = 'LSRK'
clean_params['deconvolver'] = 'hogbom'
clean_params['usemask'] = 'pb' # construct a 1/0 mask at the 0.2 level
clean_params['pbmask'] = 0.2 # data outside this pbmask will not be fitted
clean_params['threshold'] = '0.30mJy' # 2-sigma
clean_params['pblimit'] = 0.2
clean_params['pbcor'] = False
clean_params['restoration'] = True
clean_params['restoringbeam'] = 'common'
clean_params['interactive'] = False

for DataSet_dir, DataSet_name, Target_name in list(zip(DataSet_dirs, DataSet_names, Target_names)):

    vis = os.path.abspath(DataSet_dir)
    field = Target_name
    current_dir = os.getcwd()

    print2('DataSet: '+DataSet_name)
    
    # 
    imsize_list, phasecenter_list = get_mosaic_imsize_and_phasecenter(vis, imcell, galaxy_name=field, divide_into_ncol_and_nrow=(ncol,nrow), padding_by_primary_beam=1.5)
    if dry_run:
        raise NotImplementedError('This is a dry run.')

    # loop to make divided mosaic
    for i in range(len(imsize_list)):
        field_IDs = get_field_IDs_in_mosaic(vis, cell=imcell, imsize=imsize_list[i], phasecenter=phasecenter_list[i], galaxy_name=field)
        if len(field_IDs) == 0:
            continue
        outputvis = 'Level_3_Split_Divide_Mosaic/%s_Mosaic_%d_%d.ms'%(DataSet_name, i%ncol, int(i/ncol))
        if not os.path.isdir(os.path.dirname(outputvis)):
            os.makedirs(os.path.dirname(outputvis))
        if not os.path.isdir(outputvis):
            split_params = {}
            split_params['vis'] = vis
            split_params['outputvis'] = outputvis
            split_params['field'] = ','.join([str(t) for t in field_IDs])
            split_params['datacolumn'] = get_datacolumn(vis)
            print_params(split_params, 'split')
            split(**split_params)
        if not os.path.isdir(outputvis):
            raise Exception('Error! Failed to run split!')
        #break
    
    # loop to fixvis
    for i in range(len(imsize_list)):
        vis = 'Level_3_Split_Divide_Mosaic/%s_Mosaic_%d_%d.ms'%(DataSet_name, i%ncol, int(i/ncol))
        outputvis = 'Level_3_Split_Divide_Mosaic/%s_Mosaic_%d_%d.fixvis.ms'%(DataSet_name, i%ncol, int(i/ncol))
        if not os.path.isdir(vis):
            continue
        if not os.path.isdir(outputvis):
            fixvis_params = {}
            fixvis_params['vis'] = vis
            fixvis_params['outputvis'] = outputvis
            fixvis_params['field'] = field
            fixvis_params['phasecenter'] = phasecenter_list[i]
            print_params(fixvis_params, 'fixvis')
            fixvis(**fixvis_params)

    # loop to run tclean
    for i in range(len(imsize_list)):
        vis = 'Level_3_Split_Divide_Mosaic/%s_Mosaic_%d_%d.fixvis.ms'%(DataSet_name, i%ncol, int(i/ncol))
        imagename = 'Level_4_Data_Images_Divide_Mosaic/%s_Mosaic_%d_%d/output_%s_dirty'%(DataSet_name, i%ncol, int(i/ncol), field)
        if not os.path.isdir(vis):
            continue
        if not os.path.isdir(os.path.dirname(imagename)):
            os.makedirs(os.path.dirname(imagename))
        if not os.path.isdir(imagename+'.image'):
            clean_params['vis'] = vis
            clean_params['field'] = field
            clean_params['phasecenter'] = phasecenter_list[i]
            clean_params['cell'] = imcell
            clean_params['imsize'] = imsize_list[i]
            clean_params['imagename'] = imagename
            clean_params['niter'] = 0
            print2('Making dirty continuum image ...')
            print2('tclean('+', '.join("{!s}={!r}".format(k, clean_params[k]) for k in clean_params.keys())+')')
            cleanup_tclean_products(imagename)
            tclean(**clean_params)
            if not os.path.isdir(imagename+'.image'):
                raise Exception('Error! Failed to make the dirty continuum image.')
            if not os.path.isfile(imagename+'.image.fits'):
                export_tclean_products_as_fits_files(imagename, velocity=False)
            if not os.path.isfile(imagename+'.image.pbcor.fits'):
                apply_pbcor_to_tclean_image(imagename, velocity=False)
        #break


print('Done!')








