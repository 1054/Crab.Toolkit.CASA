
# RUN THIS CODE INSIDE CASA

import os, sys, re, copy, glob, shutil, datetime
if not (os.path.expanduser('~/Cloud/Github/Crab.Toolkit.CASA/lib/python') in sys.path):
    sys.path.insert(1, os.path.expanduser('~/Cloud/Github/Crab.Toolkit.CASA/lib/python'))
import dzliu_clean_utils
#reload(dzliu_clean_utils)
from dzliu_clean_utils import (get_datacolumn, get_field_IDs_in_mosaic, get_mstransform_params_for_spectral_line,
    get_synbeam_and_imcell, apply_pbcor_to_tclean_image, export_tclean_products_as_fits_files,
    cleanup_tclean_products, _print_params, load_params_from_dot_last_file)
from collections import OrderedDict
import numpy as np

# galaxy info
galaxy_name = 'PJ011646.8'
galaxy_redshift = 2.1249
line_name = 'CO32'
line_restfreq_GHz = 345.7959899
line_width_kms = 1200.0
cont_width_kms = 1200.0 # to add to both line wings
chan_width_kms = 30.0
phasecenter = 'J2000 19.194875deg -24.617194deg'
fov_arcsec = 9.0
run_tclean_dir = 'run_tclean_perplanebeam_with_cleanmask'

if not os.path.isdir(run_tclean_dir):
    os.mkdir(run_tclean_dir)


# load tclean params
prev_tclean_params_file = 'run_tclean_perplanebeam/PJ011646.8_CO32_nat.image.tclean.last'
tclean_params = load_params_from_dot_last_file(prev_tclean_params_file)

tclean_params['usemask'] = 'user'
tclean_params['mask'] = os.path.join(run_tclean_dir, 'cleanmask.mask')  # 'mask.mask' # 'out_signalmask_freq_axis_CO32.fits'


# load rms
prev_rms_file = tclean_params['imagename']+'.residual.rms.txt' # 'run_tclean_perplanebeam/PJ011646.8_CO32_nat.residual.rms.txt'
rms = float(np.loadtxt(prev_rms_file))


# prepare cleanmask
inp_clean_mask_file = 'out_signalmask_freq_axis_CO32.fits'
prev_image_file = tclean_params['imagename']+'.fits'
print('>>>', os.path.exists(tclean_params['mask']))
if not os.path.exists(tclean_params['mask']):
    importfits(fitsimage = tclean_params['imagename']+'.image.fits', imagename = os.path.join(run_tclean_dir, 'cleanmasktemplateimage.image'), defaultaxes=True, defaultaxesvalues=['0.0h', '0.0deg', '100GHz', 'I'], overwrite=True)
    importfits(fitsimage = inp_clean_mask_file, imagename = os.path.join(run_tclean_dir, 'cleanmaskinputmask.image'), defaultaxes=True, defaultaxesvalues=['0.0h', '0.0deg', '100GHz', 'I'], overwrite=True)
    makemask(mode = 'copy', inpimage = os.path.join(run_tclean_dir, 'cleanmasktemplateimage.image'), inpmask = os.path.join(run_tclean_dir, 'cleanmaskinputmask.image'), output = tclean_params['mask'])


# set up image_name
image_name = os.path.join(run_tclean_dir, galaxy_name+'_'+line_name+'_nat')
tclean_params['imagename'] = image_name


# imaging natural multiscale
tclean_params['weighting'] = 'natural'
tclean_params['niter'] = 3000
tclean_params['threshold'] = 3.0 * rms
tclean_params['deconvolver'] = 'multiscale'
tclean_params['scales'] = [0,3,10,30]
if not os.path.exists(image_name+'_multiscale.image.fits'):
    #tclean_params['imagename'] = image_name
    _print_params(tclean_params, 'tclean')
    cleanup_tclean_products(image_name)
    tclean(**tclean_params)
    shutil.copy2('tclean.last', image_name+'.image.tclean.last')
    shutil.copy2('tclean.last', image_name+'_multiscale'+'.image.tclean.last')
    for data_type in ['.image', '.mask', '.model', '.pb', '.psf', '.residual', '.sumwt', '.weight']:
        if os.path.isdir(image_name+'_multiscale'+data_type):
            shutil.rmtree(image_name+'_multiscale'+data_type)
        shutil.copytree(image_name+data_type, image_name+'_multiscale'+data_type)
    export_tclean_products_as_fits_files(image_name+'_multiscale')
    apply_pbcor_to_tclean_image(image_name+'_multiscale')


# imaging natural singlescale
tclean_params['niter'] = 10000
tclean_params['threshold'] = 1.0 * rms
tclean_params['deconvolver'] = 'hogbom'
tclean_params['calcpsf'] = False
tclean_params['mask'] = ''
if True:
    #tclean_params['imagename'] = image_name
    _print_params(tclean_params, 'tclean')
    cleanup_tclean_products(image_name)
    for data_type in ['.mask', '.model', '.pb', '.psf', '.residual', '.sumwt', '.weight']:
        shutil.copytree(image_name+'_multiscale'+data_type, image_name+data_type)
    
    tclean(**tclean_params)
    
    shutil.copy2('tclean.last', image_name+'.image.tclean.last')
    shutil.copy2('tclean.last', image_name+'_singlescale'+'.image.tclean.last')
    for data_type in ['.image', '.mask', '.model', '.pb', '.psf', '.residual', '.sumwt', '.weight']:
        if os.path.isdir(image_name+'_singlescale'+data_type):
            shutil.rmtree(image_name+'_singlescale'+data_type)
        shutil.copytree(image_name+data_type, image_name+'_singlescale'+data_type)
    export_tclean_products_as_fits_files(image_name+'_singlescale')
    apply_pbcor_to_tclean_image(image_name+'_singlescale')

    export_tclean_products_as_fits_files(image_name)
    apply_pbcor_to_tclean_image(image_name)







