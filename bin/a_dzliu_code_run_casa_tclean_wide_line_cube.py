
# RUN THIS SCRIPT INSIDE CASA

# Input a *.ms or *.uvfits data file, 
# Run CASA tclean and output image cube and other files.

import os, sys, re, shutil
import numpy as np
from collections import OrderedDict
sys.path.append(os.path.expanduser('~/Cloud/Github/Crab.Toolkit.CASA/lib/python'))
sys.path.append(os.path.expanduser('~/Cloud/Github/Crab.Toolkit.CASA/lib/python/analysis_scripts'))
from dzliu_clean_utils import (get_datacolumn, get_synbeam_and_imcell, get_mstransform_params_for_spectral_line, 
                               apply_pbcor_to_tclean_image, export_tclean_products_as_fits_files, cleanup_tclean_products)
import analysisUtils as au
from taskinit import casalog #, tb, ms, iatool
#import casadef
#def version_tuple(version_str):
#    return tuple(map(int, (version_str.split("."))))
#def version_less_than(version_str, compared_version_str):
#    return version_tuple(version_str) < version_tuple(compared_version_str)
#def version_greater_equal(version_str, compared_version_str):
#    return version_tuple(version_str) >= version_tuple(compared_version_str)
#if version_less_than(casadef.casa_version, '6.0.0'):
#    #from __main__ import default, inp, saveinputs
#    ##import task_tclean; task_tclean.tclean # this is what tclean.py calls
#    ##import tclean
#    ##import tclean_cli
#    from tclean_cli import tclean_cli_
#    tclean = tclean_cli_()
#    from mstransform_cli import mstransform_cli_
#    mstransform = mstransform_cli_()
#    from uvcontsub_cli import uvcontsub_cli_
#    uvcontsub = uvcontsub_cli_()
#    from exportfits_cli import exportfits_cli_
#    exportfits = exportfits_cli_()
#    from importfits_cli import importfits_cli_
#    importfits = importfits_cli_()
#    from concat_cli import concat_cli_
#    concat = concat_cli_()
#    from split_cli import split_cli_
#    split = split_cli_()
#    from imstat_cli import imstat_cli_
#    imstat = imstat_cli_()
#    from simobserve_cli import simobserve_cli_
#    simobserve = simobserve_cli_()
#else:
#    # see CASA 6 updates here: https://alma-intweb.mtk.nao.ac.jp/~eaarc/UM2018/presentation/Nakazato.pdf
#    from casatasks import tclean, mstransform, uvcontsub, exportfits, importfits, concat, split, imstat
#    #from casatasks import sdbaseline
#    #from casatools import ia



#########################
# User input parameters #
#########################
#input_vis = 'input_vis.ms'
#input_field = ''
#output_prefix = 'output_wide_line'
input_vis = os.getenv('CASA_COMMANDLINE_VIS')
input_field = os.getenv('CASA_COMMANDLINE_FIELD')
output_prefix = os.getenv('CASA_COMMANDLINE_OUT')
if input_vis is None:
    raise Exception('Error! input_vis is None!')
if input_field is None:
    input_field = ''
if output_prefix is None:
    output_prefix = re.sub(r'\.(ms|uvfits)$', r'_tclean_line', input_vis)
phase_center = ''
max_imsize = 1000
do_dirty = True
do_clean = True # down to 3-sigma, no mask
do_natural = False # down to 2-sigma, with aperture mask if applicable
do_briggs = False
do_superuniform = False
overwrite = False






#########################
# Define some functions #
#########################
# global casalog_origin2
global casalog_origin2
global casalog_origin
casalog_origin2 = 'dzliu_run_casa_tclean_wide_line_cube'
casalog_origin = 'dzliu_run_casa_tclean_wide_line_cube'
def set_casalog_origin(origin):
    global casalog_origin2
    global casalog_origin
    casalog_origin2 = casalog_origin
    casalog_origin = origin
    casalog.origin(casalog_origin)
def restore_casalog_origin():
    global casalog_origin2
    global casalog_origin
    casalog_origin = casalog_origin2
    casalog.origin(casalog_origin)
def print2(message):
    print(message)
    casalog.post(message, 'INFO')



####################
# start processing #
####################

# set_casalog_origin
set_casalog_origin(casalog_origin)


# get abspath
input_vis = os.path.abspath(input_vis)


# make working_dir
if output_prefix.find(os.sep) >= 0:
    working_dir = os.path.dirname(output_prefix)
    output_prefix = os.path.basename(output_prefix)
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)
else:
    working_dir = os.getcwd()


# change to working dir
current_dir = os.getcwd()
print2("os.chdir(%r)"%(working_dir))
os.chdir(working_dir)


# if user has input a uvfits, then importuvfits
if input_vis.endswith('.uvfits'):
    working_uvfits = input_vis
    working_vis = re.sub(r'\.uvfits$', r'.ms', os.path.basename(working_uvfits))
    # check if the ms data directory and file content already exist
    if not (os.path.isdir(working_vis) and os.path.isfile(working_vis+os.sep+'table.dat')):
        if os.path.isdir(working_vis):
            shutil.rmtree(working_vis)
        importuvfits_params = OrderedDict()
        importuvfits_params['fitsfile'] = working_uvfits
        importuvfits_params['vis'] = working_vis
        importuvfits_params['antnamescheme'] = 'new' # for VLA/EVLA/CARMA only; ’new’ or ’old’; ’VA04’ or ’04’ for VLA ant 4
        print2('Running importuvfits('+', '.join(["{!s}={!r}".format(k, importuvfits_params[k]) for k in importuvfits_params.keys()])+')')
        importuvfits(**importuvfits_params)
        shutil.copy2('importuvfits.last', working_vis+'.importuvfits.last')
    if not (os.path.isdir(working_vis) and os.path.isfile(working_vis+os.sep+'table.dat')):
        raise Exception('Error! Could not produce "%s"!'%(os.path.abspath(working_vis)))


# perform continuum subtraction <TODO> need a line channel selection
#if not (os.path.isdir(working_vis+'.contsub') and os.path.isdir(working_vis+'.contsub'+os.sep+'table.dat')):
#    if os.path.isdir(working_vis+'.contsub'):
#        shutil.rmtree(working_vis+'.contsub')
#    line_left_wing_freq_GHz = line_obsfreq_GHz - line_width_kms/2.0/2.99792458e5*line_obsfreq_GHz
#    line_right_wing_freq_GHz = line_obsfreq_GHz + line_width_kms/2.0/2.99792458e5*line_obsfreq_GHz
#    uvcontsub_params = OrderedDict()
#    uvcontsub_params['vis'] = working_vis
#    uvcontsub_params['fitspw'] = '0:%.7f~%.7fGHz'%(line_left_wing_freq_GHz, line_right_wing_freq_GHz)
#    uvcontsub_params['excludechans'] = True
#    uvcontsub_params['want_cont'] = True
#    print2('Running uvcontsub('+', '.join(["{!s}={!r}".format(k, uvcontsub_params[k]) for k in uvcontsub_params.keys()])+')')
#    uvcontsub(**uvcontsub_params)
#    shutil.copy2('uvcontsub.last', uvcontsub_params['vis']+'.contsub'+'.uvcontsub.last')
#    if not (os.path.isdir(working_vis+'.contsub') and os.path.isfile(working_vis+'.contsub'+os.sep+'table.dat')):
#        raise Exception('Error! Could not produce "%s"!'%(os.path.abspath(working_vis+'.contsub')))


# pick cell size for imaging
print2('au.pickCellSize')
au_cellsize, au_imsize, _ = au.pickCellSize(working_vis,
                                            imsize=True,
                                            npix=5,
                                            intent='',
                                            pblevel=0.5,
                                            )
imsize = au_imsize
cell = str(au_cellsize)+'arcsec'
print2('imsize: %s'%(imsize))
print2('cell: %s'%(cell))

if imsize[0] > max_imsize:
    print2('Warning! imsize[0] %s is larger than max_imsize %d, setting it to the max size.'%(imsize[0], max_imsize))
    imsize[0] = max_imsize
if imsize[1] > max_imsize:
    print2('Warning! imsize[1] %s is larger than max_imsize %d, setting it to the max size.'%(imsize[1], max_imsize))
    imsize[1] = max_imsize


# prepare tclean params
tclean_params = OrderedDict()
tclean_params['vis'] = working_vis
tclean_params['imagename'] = output_prefix+'_dirty_cube'
tclean_params['field'] = input_field
tclean_params['spw'] = ''
tclean_params['phasecenter'] = phase_center
tclean_params['cell'] = cell
tclean_params['imsize'] = imsize
#tclean_params['reffreq'] = '%.7fGHz'%(line_obsfreq_GHz)
#tclean_params['restfreq'] = '%.7fGHz'%(line_obsfreq_GHz)
tclean_params['width'] = 1
tclean_params['specmode'] = 'cube'
#tclean_params['restoringbeam'] = 'common'
tclean_params['weighting'] = 'natural'
tclean_params['niter'] = 0

# run tclean to make dirty image
if not os.path.isdir(tclean_params['imagename']+'.image'):
    print2('Running tclean('+', '.join(["{!s}={!r}".format(k, tclean_params[k]) for k in tclean_params.keys()])+')')
    cleanup_tclean_products(tclean_params['imagename']+'.image')
    tclean(**tclean_params)
    shutil.copy2('tclean.last', tclean_params['imagename']+'.image'+'.tclean.last')
# apply pbcorr
if not os.path.isdir(tclean_params['imagename']+'.image.pbcor'):
    apply_pbcor_to_tclean_image(tclean_params['imagename']+'.image')
# export as fits files
if not os.path.isfile(tclean_params['imagename']+'.image.fits'):
    export_tclean_products_as_fits_files(tclean_params['imagename']+'.image')
# analyze dirty image and get rms
imstat_dict = imstat(tclean_params['imagename']+'.image', axes=[0,1])
rms = np.nanmedian(imstat_dict['rms'])
with open('rms_dirty_cube.txt', 'w') as fp:
    fp.write('%.7g\n'%(rms))
print2('RMS from dirty image: %.7g'%(rms))



# set threshold for next cleaning
tclean_params['threshold'] = '%.7gJy'%(3.0*rms)
tclean_params['imagename'] = output_prefix+'_clean_cube'
tclean_params['niter'] = 30000

# run tclean to make clean image
if not os.path.isdir(tclean_params['imagename']+'.image'):
    print2('Running tclean('+', '.join(["{!s}={!r}".format(k, tclean_params[k]) for k in tclean_params.keys()])+')')
    cleanup_tclean_products(tclean_params['imagename']+'.image')
    tclean(**tclean_params)
    shutil.copy2('tclean.last', tclean_params['imagename']+'.image'+'.tclean.last')
# apply pbcorr
if not os.path.isdir(tclean_params['imagename']+'.image.pbcor'):
    apply_pbcor_to_tclean_image(tclean_params['imagename']+'.image')
# export as fits files
if not os.path.isfile(tclean_params['imagename']+'.image.fits'):
    export_tclean_products_as_fits_files(tclean_params['imagename']+'.image')
# analyze clean residual image and get rms
imstat_dict = imstat(tclean_params['imagename']+'.residual', axes=[0,1])
rms = np.nanmedian(imstat_dict['rms'])
with open('rms_clean_residual.txt', 'w') as fp:
    fp.write('%.7g\n'%(rms))
with open('rms.txt', 'w') as fp:
    fp.write('%.7g\n'%(rms))
print2('RMS from clean residual image: %.7g'%(rms))



# set threshold for next cleaning
tclean_params['threshold'] = '%.7gJy'%(2.0*rms)
tclean_params['imagename'] = output_prefix+'_clean_cube_natural_weighting'
tclean_params['niter'] = 50000
tclean_params['weighting'] = 'natural'

# set mask for next cleaning if applicable <TODO>
#center_aperture_radius_arcsec = clean_mask_radius_arcsec
#center_aperture_radius_pix = center_aperture_radius_arcsec/au_cellsize
#if center_aperture_radius_pix < np.min(imsize)/2.0:
#    tclean_params['usemask'] = 'user'
#    tclean_params['mask'] = 'circle[[%dpix,%dpix],%dpix]'%(int(np.round(imsize[0]+1.0/2.0)), int(np.round(imsize[1]+1.0/2.0)), int(np.ceil(center_aperture_radius_pix)))

# run tclean
if do_natural:
    # run tclean to make clean image
    if not os.path.isdir(tclean_params['imagename']+'.image'):
        print2('Running tclean('+', '.join(["{!s}={!r}".format(k, tclean_params[k]) for k in tclean_params.keys()])+')')
        cleanup_tclean_products(tclean_params['imagename']+'.image')
        tclean(**tclean_params)
        shutil.copy2('tclean.last', tclean_params['imagename']+'.image'+'.tclean.last')
    # apply pbcorr
    if not os.path.isdir(tclean_params['imagename']+'.image.pbcor'):
        apply_pbcor_to_tclean_image(tclean_params['imagename']+'.image')
    # export as fits files
    if not os.path.isfile(tclean_params['imagename']+'.image.fits'):
        export_tclean_products_as_fits_files(tclean_params['imagename']+'.image')
    # analyze clean residual image and get rms
    imstat_dict = imstat(tclean_params['imagename']+'.residual', axes=[0,1])
    rms = np.nanmedian(imstat_dict['rms'])
    print2('RMS from clean residual image: %.7g'%(rms))



# set threshold for next cleaning
tclean_params['threshold'] = '%.7gJy'%(2.0*rms)
tclean_params['imagename'] = output_prefix+'_clean_cube_briggs_weighting'
tclean_params['niter'] = 50000
tclean_params['weighting'] = 'briggs'
tclean_params['robust'] = 0.5

# run tclean
if do_briggs:
    # run tclean to make clean image
    if not os.path.isdir(tclean_params['imagename']+'.image'):
        print2('Running tclean('+', '.join(["{!s}={!r}".format(k, tclean_params[k]) for k in tclean_params.keys()])+')')
        cleanup_tclean_products(tclean_params['imagename']+'.image')
        tclean(**tclean_params)
        shutil.copy2('tclean.last', tclean_params['imagename']+'.image'+'.tclean.last')
    # apply pbcorr
    if not os.path.isdir(tclean_params['imagename']+'.image.pbcor'):
        apply_pbcor_to_tclean_image(tclean_params['imagename']+'.image')
    # export as fits files
    if not os.path.isfile(tclean_params['imagename']+'.image.fits'):
        export_tclean_products_as_fits_files(tclean_params['imagename']+'.image')
    # analyze clean residual image and get rms
    imstat_dict = imstat(tclean_params['imagename']+'.residual', axes=[0,1])
    rms = np.nanmedian(imstat_dict['rms'])
    print2('RMS from clean residual image: %.7g'%(rms))



# set threshold for next cleaning
tclean_params['threshold'] = '%.7gJy'%(2.0*rms)
tclean_params['imagename'] = output_prefix+'_clean_cube_superuniform_weighting'
tclean_params['niter'] = 50000
tclean_params['weighting'] = 'superuniform'

# run tclean
if do_superuniform:
    # run tclean to make clean image
    if not os.path.isdir(tclean_params['imagename']+'.image'):
        print2('Running tclean('+', '.join(["{!s}={!r}".format(k, tclean_params[k]) for k in tclean_params.keys()])+')')
        cleanup_tclean_products(tclean_params['imagename']+'.image')
        tclean(**tclean_params)
        shutil.copy2('tclean.last', tclean_params['imagename']+'.image'+'.tclean.last')
    # apply pbcorr
    if not os.path.isdir(tclean_params['imagename']+'.image.pbcor'):
        apply_pbcor_to_tclean_image(tclean_params['imagename']+'.image')
    # export as fits files
    if not os.path.isfile(tclean_params['imagename']+'.image.fits'):
        export_tclean_products_as_fits_files(tclean_params['imagename']+'.image')
    # analyze clean residual image and get rms
    imstat_dict = imstat(tclean_params['imagename']+'.residual', axes=[0,1])
    rms = np.nanmedian(imstat_dict['rms'])
    print2('RMS from clean residual image: %.7g'%(rms))













