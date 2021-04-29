
# Run this script in CASA

import os, sys, re, shutil
import numpy as np
sys.path.append(os.path.expanduser('~/Cloud/Github/Crab.Toolkit.CASA/lib/python'))
sys.path.append(os.path.expanduser('~/Cloud/Github/Crab.Toolkit.CASA/lib/python/analysis_scripts'))
from dzliu_clean_utils import (get_datacolumn, get_spw_for_spectral_line, cleanup_tclean_products, apply_pbcor_to_tclean_image, export_tclean_products_as_fits_files)
import analysisUtils as au

target = 'SDSSJ0901+1814'
redshift = 2.2586
# linefreq = 35.37445584 # 115.2712018/(1.+redshift)
# ch0 35298.519 GHz, ChanWid 0.5 MHz, TotBW 128 MHz, CtrFreq 35362.2692 MHz. 
vis = 'Level_3_Split/DataSet_Merged_Ka_Band/merged_'+target+'.ms'
outname = 'Level_4_Data_Images/'+target+'/DataSet_Merged_Ka_Band/merged_'+target
fitspw = '0:4~24;55~62'
fitorder = 0
want_cont = True


# make output dir
if not os.path.isdir(os.path.dirname(outname)):
    os.makedirs(os.path.dirname(outname))


# uvcontsub
print('vis: %r'%(vis))
if not os.path.isdir(vis+'.contsub'):
    print('Running uvcontsub(vis=%r)'%(vis))
    uvcontsub()
    shutil.copy2('uvcontsub.last', vis+'.contsub'+'.uvcontsub.last')
if os.path.isdir(vis+'.contsub'):
    vis = vis+'.contsub'
    print('vis: %r'%(vis))



# tclean preparation
au_cellsize, au_imsize, _ = au.pickCellSize(vis,
                                            imsize=True,
                                            npix=5,
                                            intent='',
                                            pblevel=0.2,
                                            )
imsize = au_imsize
cell = str(au_cellsize)+'arcsec'
specmode = 'cube'
reffreq = '%.9fGHz'%(115.2712018/(1.+redshift))
outframe = 'LSRK'
veltype = 'radio'
restfreq = '115.2712018GHz'
restoration = True
restoringbeam = 'common'


# make dirty image
threshold = '0mJy'
niter = 0
imagename = outname + '_cube_dirty'
if not os.path.isdir(imagename+'.image'):
    print('Making dirty image ...')
    #inp(tclean)
    tclean()
    shutil.copy2('tclean.last', imagename+'.image.tclean.last')
if not os.path.isdir(imagename+'.image.pbcor'):
    apply_pbcor_to_tclean_image(imagename+'.image')
if not os.path.isfile(imagename+'.image.fits'):
    export_tclean_products_as_fits_files(imagename+'.image')


# analyze dirty image
imstat_dict = imstat(imagename+'.image', axes=[0,1])
rms = np.median(imstat_dict['rms'])


# make clean image 
threshold = '%gJy'%(3.5*rms)
niter = 30000
datacolumn = get_datacolumn(vis)
print('rms: %s'%(rms))
print('threshold: %s'%(threshold))
imagename = outname + '_cube_clean_shallow'
if not os.path.isdir(imagename+'.image'):
    print('Making clean image ...')
    cleanup_tclean_products(imagename+'.image')
    inp(tclean)
    tclean()
    shutil.copy2('tclean.last', imagename+'.image.tclean.last')
if not os.path.isdir(imagename+'.image.pbcor'):
    apply_pbcor_to_tclean_image(imagename+'.image')
if not os.path.isfile(imagename+'.image.fits'):
    export_tclean_products_as_fits_files(imagename+'.image')


# analyze clean residual image
imstat_dict = imstat(imagename=imagename+'.residual', axes=[0,1])
rms = np.median(imstat_dict['rms'])


# make clean image deep in mask
threshold = '%gJy'%(2.0*rms)
niter = 30000
usemask = 'user'
mask = 'ellipse[[9h01m22.5001, +18d14m31.677], [8.190arcsec, 11.241deg], 0deg]'
weighting = 'briggs'
robust = 0.5
restart = True
previmagename = imagename
datacolumn = get_datacolumn(vis)
print('rms: %s'%(rms))
print('threshold: %s'%(threshold))
print('mask: %s'%(mask))
previmagename = imagename
imagename = outname + '_cube_clean_deep'
if not os.path.isdir(imagename+'.image'):
    print('Making clean image ...')
    cleanup_tclean_products(imagename+'.image')
    for suffix in ['.image', '.model', '.residual', '.psf', '.pb', '.sumwt']:
        shutil.copytree(previmagename+suffix, imagename+suffix)
    inp(tclean)
    tclean()
    shutil.copy2('tclean.last', imagename+'.image.tclean.last')
if not os.path.isdir(imagename+'.image.pbcor'):
    apply_pbcor_to_tclean_image(imagename+'.image')
if not os.path.isfile(imagename+'.image.fits'):
    export_tclean_products_as_fits_files(imagename+'.image')





