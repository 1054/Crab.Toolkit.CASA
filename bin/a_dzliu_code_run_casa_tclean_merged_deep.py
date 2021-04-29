
# RUN THIS SCRIPT INSIDE CASA

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
target = 'GS30274'
redshift = 2.225
line_name = 'CO32'
line_width_kms = 800.0
line_free_kms = 800.0 # line free channels at each side of the line
line_restfreq_GHz = 345.7959899
line_obsfreq_GHz = line_restfreq_GHz/(1.0+redshift)
channel_width_kms = 20.0
clean_mask_radius_arcsec = 2.0 #<TODO># galaxy must be at the image center
phase_center = ''
max_imsize = 1000
#max_line_width_kms = 1.875/(line_restfreq_GHz/(1.0+redshift))*2.99792458e5

vis_list = [\
    'Level_3_Split/DataSet_01/calibrated.ms',
    'Level_3_Split/DataSet_02/calibrated.ms',
]
target_list = [\
    target,
    target,
]





#########################
# Define some functions #
#########################
def print2(message):
    print(message)
    casalog.post(message, 'INFO')



####################
# start processing #
####################
working_dir = 'Level_4_Data_Images/'+target+'/DataSet_Merged'

# make working dir
if not os.path.isdir(working_dir):
    os.makedirs(working_dir)

# mstransform, uvcontsub and concat vis
concatenated_vis = os.path.join(working_dir, 'concatenated.ms')
if not os.path.isdir(concatenated_vis):
    concat_params = OrderedDict()
    concat_params['vis'] = []
    concat_params['concatvis'] = 'merged.ms'
    concat_params['freqtol'] = '%.4fMHz'%(chan_width_kms*0.8/2.99792458e5*line_obsfreq_GHz*1e3)
    concat_params['respectname'] = False
    concat_params['copypointing'] = False
    for ivis, vis in enumerate(vis_list):
        this_vis = os.path.join(working_dir, 'mstransformed_vis_%d_.ms'%(ivis+1))
        if not os.path.isdir(this_vis):
            mstransform_params = \
                get_mstransform_params_for_spectral_line(\
                    vis, 
                    this_vis, 
                    field=target_list[ivis], 
                    redshift=redshift, 
                    rest_freq_GHz=line_restfreq_GHz, 
                    line_width_kms=line_width_kms+line_free_kms+line_free_kms,
                    chan_width_kms=channel_width_kms,
                    force_integer_chan_width=False,
                    )
            print2('Running mstransform('+', '.join(["{!s}={!r}".format(k, mstransform_params[k]) for k in mstransform_params.keys()])+')')
            mstransform(**mstransform_params)
            shutil.copy2('mstransform.last', mstransform_params['outputvis']+'.mstransform.last')
            # 
        if not os.path.isdir(this_vis+'.contsub'):
            line_left_wing_freq_GHz = line_obsfreq_GHz - line_width_kms/2.0/2.99792458e5*line_obsfreq_GHz
            line_right_wing_freq_GHz = line_obsfreq_GHz + line_width_kms/2.0/2.99792458e5*line_obsfreq_GHz
            uvcontsub_params = OrderedDict()
            uvcontsub_params['vis'] = this_vis
            uvcontsub_params['fitspw'] = '0:%.7f~%.7fGHz'%(line_left_wing_freq_GHz, line_right_wing_freq_GHz)
            uvcontsub_params['excludechans'] = True
            uvcontsub_params['want_cont'] = True
            print2('Running uvcontsub('+', '.join(["{!s}={!r}".format(k, uvcontsub_params[k]) for k in uvcontsub_params.keys()])+')')
            uvcontsub(**uvcontsub_params)
            shutil.copy2('uvcontsub.last', uvcontsub_params['vis']+'.contsub'+'.uvcontsub.last')
            # 
        concat_params['vis'].append(this_vis+'.contsub')
    # 
    print2('Running concat('+', '.join(["{!s}={!r}".format(k, concat_params[k]) for k in concat_params.keys()])+')')
    concat(**concat_params)
    shutil.copy2('concat.last', concat_params['concatvis']+'.concat.last')
    # 
    if not os.path.isdir(concatenated_vis):
        raise Exception('Error! Could not produce "%s"!'%(os.path.abspath(concatenated_vis)))


# change to working dir
current_dir = os.getcwd()
print2("os.chdir(%r)"%(working_dir))
os.chdir(working_dir)


# use vis basename
vis = os.path.basename(concatenated_vis)


# pick cell size
print2('au.pickCellSize')
au_cellsize, au_imsize, _ = au.pickCellSize(vis,
                                            imsize=True,
                                            npix=5,
                                            intent='',
                                            pblevel=0.5,
                                            )
imsize = au_imsize
cell = str(au_cellsize)+'arcsec'
print2('imsize: %d'%(imsize))
print2('cell: %s'%(cell))

if imsize > max_imsize:
    print2('Warning! imsize %d is larger than max_imsize %d, setting it to the max size.'%(imsize, max_imsize))


# prepare tclean params
tclean_params = OrderedDict()
tclean_params['vis'] = vis
tclean_params['imagename'] = '%s_%s_line_dirty_cube'%(target, line_name)
tclean_params['field'] = ''
tclean_params['spw'] = ''
tclean_params['phasecenter'] = phase_center
tclean_params['cell'] = cell
tclean_params['imsize'] = imsize
tclean_params['reffreq'] = line_obsfreq_GHz
tclean_params['restfreq'] = line_obsfreq_GHz
tclean_params['width'] = 1
tclean_params['specmode'] = 'cube'
tclean_params['restoringbeam'] = 'common'
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
print2('RMS from dirty image: %.7g'%(rms))


# set threshold for next cleaning
tclean_params['threshold'] = '%.7gJy'%(3.0*rms)
tclean_params['imagename'] = '%s_%s_line_clean_cube'%(target, line_name)
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
with open('rms.txt', 'w') as fp:
    fp.write('%.7g\n'%(rms))
print2('RMS from clean residual image: %.7g'%(rms))


# set threshold for next cleaning
tclean_params['threshold'] = '%.7gJy'%(2.0*rms)
tclean_params['imagename'] = '%s_%s_line_clean_cube_natural_weighting'%(target, line_name)
tclean_params['niter'] = 50000
tclean_params['weighting'] = 'natural'
center_aperture_radius_arcsec = clean_mask_radius_arcsec
center_aperture_radius_pix = center_aperture_radius_arcsec/au_cellsize
if center_aperture_radius_pix < imsize/2.0:
    tclean_params['usemask'] = 'user'
    tclean_params['mask'] = 'circle[[%dpix,%dpix],%dpix]'%(int(np.round(imsize+1.0/2.0)), int(np.round(imsize+1.0/2.0)), int(np.ceil(center_aperture_radius_pix)))

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
rms = np.nanmedian(imstat_dict['rms'])



# set threshold for next cleaning
tclean_params['threshold'] = '%.7gJy'%(2.0*rms)
tclean_params['imagename'] = '%s_%s_line_clean_cube_briggs_weighting'%(target, line_name)
tclean_params['niter'] = 50000
tclean_params['weighting'] = 'briggs'
tclean_params['robust'] = 0.5

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
rms = np.nanmedian(imstat_dict['rms'])



# set threshold for next cleaning
tclean_params['threshold'] = '%.7gJy'%(2.0*rms)
tclean_params['imagename'] = '%s_%s_line_clean_cube_superuniform_weighting'%(target, line_name)
tclean_params['niter'] = 50000
tclean_params['weighting'] = 'superuniform'

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













