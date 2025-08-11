
# RUN THIS INSIDE CASA

import os, sys, re, copy, glob, shutil, datetime
if not (os.path.expanduser('~/Cloud/Github/Crab.Toolkit.CASA/lib/python') in sys.path):
    sys.path.insert(1, os.path.expanduser('~/Cloud/Github/Crab.Toolkit.CASA/lib/python'))
import dzliu_clean_utils
#reload(dzliu_clean_utils)
from dzliu_clean_utils import (get_datacolumn, get_field_IDs_in_mosaic, get_mstransform_params_for_spectral_line, 
    get_synbeam_and_imcell, apply_pbcor_to_tclean_image, export_tclean_products_as_fits_files, 
    cleanup_tclean_products, _print_params)
from collections import OrderedDict
import numpy as np

# ALMA projects to combine
# calibrated ms data must exist under '../../../uvdata/'+project+'/s*/g*/m*/calibrated/calibrated.ms'
list_of_projects = [\
    '2018.1.00543.S', 
    '2015.1.00228.S',
    '2015.1.01379.S',
]

# galaxy info
galaxy_name = 'K20-ID7'
galaxy_redshift = 2.224
line_name = 'CO32'
line_restfreq_GHz = 345.7959899
line_width_kms = 800.0
cont_width_kms = 1200.0 # to add to both line wings
chan_width_kms = 30.0
phasecenter = 'J2000 53.131083deg -27.773108deg'
fov_arcsec = 25.0
contsub_lines = [
    dict(name=line_name, restfreq=line_restfreq_GHz, linewidth=line_width_kms), 
    #dict(name='CI21', restfreq=809.34197, linewidth=800.),
    #dict(name='CO76', restfreq=806.6518060, linewidth=1200.),
]

run_tclean_dir = 'run_tclean_perplanebeam'

if not os.path.isdir(run_tclean_dir):
    os.makedirs(run_tclean_dir)

# loop projects
vis_list_to_combine = []
for project in list_of_projects:

    vis_to_combine = 'input_data_to_combine_'+project

    if len(glob.glob(vis_to_combine+'*.ms')) == 0:

        # find ms
        list_of_files = []
        find_file_patterns = []
        if len(list_of_files) == 0:
            find_file = '/nfs/irdata02/datamine_radio/alma_archive/uvdata/'+project+'/s*/g*/m*/calibrated/calibrated.ms'
            list_of_files = glob.glob(find_file)
            find_file_patterns.append(find_file)
        if len(list_of_files) == 0:
            find_file = '../../'+project+'_PI_*/Level_1_Raw/'+project+'/s*/g*/m*/calibrated/calibrated_final.ms'
            list_of_files = glob.glob(find_file)
            find_file_patterns.append(find_file)
        if len(list_of_files) == 0:
            find_file = '../../'+project+'_PI_*/Level_1_Raw/'+project+'/s*/g*/m*/calibrated/calibrated.ms'
            list_of_files = glob.glob(find_file)
            find_file_patterns.append(find_file)
        if len(list_of_files) == 0:
            find_file = '../../'+project+'_PI_*/Level_1_Raw/'+project+'/s*/g*/m*/calibrated/uid___A002_*.ms.split.cal'
            list_of_files = glob.glob(find_file)
            find_file_patterns.append(find_file)
        if len(list_of_files) == 0:
            find_file = '../../'+project+'_PI_*/Level_1_Raw/'+project+'/s*/g*/m*/calibrated/uid___A002_*.ms.split'
            list_of_files = glob.glob(find_file)
            find_file_patterns.append(find_file)
        if len(list_of_files) == 0:
            find_file = '../../'+project+'_PI_*/Level_1_Raw/'+project+'/s*/g*/m*/calibrated/uid___A002_*.ms'
            list_of_files = glob.glob(find_file)
            find_file_patterns.append(find_file)
        if len(list_of_files) == 0:
            find_file = '../../../uvdata/'+project+'/s*/g*/m*/calibrated/calibrated.ms'
            list_of_files = glob.glob(find_file)
            find_file_patterns.append(find_file)
        if len(list_of_files) == 0:
            raise Exception('Error! Could not find calibrated ms in following paths:\n'+'\n'.join(find_file_patterns))
        print('list_of_files', list_of_files)
        
        vis_idx = 0
        for vis in list_of_files:
            
            # find field
            field_ids = get_field_IDs_in_mosaic(\
                vis = vis, 
                cell = '0.1arcsec', imsize = [50, 50], # representing a 5 arcsec size searching box
                phasecenter = phasecenter, 
                padding_by_primary_beam = 0.0, 
                verbose = True, 
            )

            print('field_ids', field_ids)

            if len(field_ids) == 0:
                print('Warning! Skipping vis {!r} because no matched field IDs.'.format(vis))
                continue


            # set outputvis name
            if len(list_of_files) != 1:
                vis_to_transform = vis_to_combine+'.%d.ms'%(vis_idx+1)
            else:
                vis_to_transform = vis_to_combine+'.ms'

            # mstransform
            if True:
                mstransform_params = get_mstransform_params_for_spectral_line(\
                    vis = vis, 
                    outputvis = vis_to_transform, 
                    field = ','.join([str(t) for t in field_ids]), 
                    redshift = galaxy_redshift, 
                    rest_freq_GHz = line_restfreq_GHz, 
                    line_width_kms = line_width_kms + 2.0 * cont_width_kms, 
                    chan_width_kms = chan_width_kms, 
                    force_integer_chan_width = False, 
                    exclude_continuum_spw = True, 
                    check_same_chan_width = True, 
                    verbose = True, 
                )
            if len(mstransform_params) == 0:
                continue
            _print_params(mstransform_params, 'mstransform')
            mstransform(**mstransform_params)
            shutil.copy2('mstransform.last', vis_to_transform+'.mstransform.last')


            vis_list_to_combine.append(vis_to_transform)

            vis_idx += 1

    else:

        vis_list_to_combine.extend(glob.glob(vis_to_combine+'*.ms'))

print('vis_list_to_combine', vis_list_to_combine)
# raise NotImplementedError()


# concat
vis_combined = 'combined.ms'
chan_width_MHz = (chan_width_kms/2.99792458e5)* (line_restfreq_GHz / (1.0+galaxy_redshift))
if not os.path.exists(vis_combined):
    concat_params = {'vis': vis_list_to_combine, 'concatvis': vis_combined, 'freqtol':'%.3fMHz'%(chan_width_MHz/10.0), 'copypointing':False}
    _print_params(concat_params, 'concat')
    concat(**concat_params)
    shutil.copy2('concat.last', vis_combined+'.concat.last')


# statwt
do_statwt = False
if do_statwt and not os.path.exists(vis_combined+'.statwt.done'):
    freq1 = (1.0 - (line_width_kms/2.0+cont_width_kms)/2.99792458e5) * (line_restfreq_GHz / (1.0+galaxy_redshift))
    freq2 = (1.0 - (line_width_kms/2.0)/2.99792458e5) * (line_restfreq_GHz / (1.0+galaxy_redshift))
    freq3 = (1.0 + (line_width_kms/2.0)/2.99792458e5) * (line_restfreq_GHz / (1.0+galaxy_redshift))
    freq4 = (1.0 + (line_width_kms/2.0+cont_width_kms)/2.99792458e5) * (line_restfreq_GHz / (1.0+galaxy_redshift))
    fitspw = '*:%.6f~%.6fGHz;%.6f~%.6fGHz'%(freq1, freq2, freq3, freq4)
    statwt_params = {'vis': vis_combined, 'spw': '', 'fitspw': fitspw, 'datacolumn': get_datacolumn(vis_combined)}
    _print_params(statwt_params, 'statwt')
    statwt(**statwt_params)
    shutil.copy2('statwt.last', vis_combined+'.statwt.last')
    with open(vis_combined+'.statwt.done', 'w') as fp:
        fp.write('done at %s\n'%(datetime.datetime.now()))


# uvcontsub
vis_combined_contsub = 'combined.ms.contsub'
if not os.path.exists(vis_combined_contsub):
    fullfreq1 = (1.0 - (line_width_kms/2.0+cont_width_kms)/2.99792458e5) * (line_restfreq_GHz / (1.0+galaxy_redshift))
    fullfreq2 = (1.0 + (line_width_kms/2.0+cont_width_kms)/2.99792458e5) * (line_restfreq_GHz / (1.0+galaxy_redshift))
    freqarr = np.linspace(fullfreq1, fullfreq2, num=1000)
    maskarr = np.full(len(freqarr), fill_value=False)
    freqdelta = (freqarr[1]-freqarr[0])
    contsub_freq_pairs = []
    for contsub_line in contsub_lines:
        freq1 = (1.0 - (contsub_line['linewidth']/2.0)/2.99792458e5) * (contsub_line['restfreq'] / (1.0+galaxy_redshift))
        freq2 = (1.0 + (contsub_line['linewidth']/2.0)/2.99792458e5) * (contsub_line['restfreq'] / (1.0+galaxy_redshift))
        maskarr[np.logical_and(freqarr>=(freq1-0.5*freqdelta), freqarr<=(freq2+0.5*freqdelta))] = True # mask line channels
    freq1 = None
    freq2 = None
    fitspw = ''
    for k in range(len(maskarr)):
        if maskarr[k]: # this channel is line
            if freq1 is None:
                freq1 = freqarr[k]
        else: # this channel is not line
            if freq1 is not None:
                freq2 = freqarr[k-1]
                if fitspw != '':
                    fitspw += ';'
                fitspw += '%.6f~%.6fGHz'%(freq1, freq2)
                freq1 = None
    fitspw = '*:'+fitspw
    #fitspw = '*:%.6f~%.6fGHz'%(freq1, freq2)
    uvcontsub_params = {'vis': vis_combined, 'spw': '', 'fitspw': fitspw, 'fitorder': 0, 'excludechans': True, 'combine': 'spw', 'want_cont': True}
    _print_params(uvcontsub_params, 'uvcontsub')
    uvcontsub(**uvcontsub_params)
    shutil.copy2('uvcontsub.last', vis_combined_contsub+'.uvcontsub.last')


# imaging prep
synbeam, imcell = get_synbeam_and_imcell(vis_combined_contsub)
imsize = int(np.ceil(fov_arcsec/float(imcell.replace('arcsec',''))/2)*2)
tclean_params = OrderedDict()
tclean_params['vis'] = vis_combined_contsub
tclean_params['imagename'] = ''
tclean_params['datacolumn'] = get_datacolumn(vis_combined_contsub)
tclean_params['field'] = ''
tclean_params['spw'] = ''
tclean_params['imsize'] = imsize
tclean_params['cell'] = imcell
tclean_params['phasecenter'] = phasecenter
tclean_params['specmode'] = 'cube'
tclean_params['reffreq'] = '%.6fGHz'%(line_restfreq_GHz / (1.0+galaxy_redshift)) # precision 1 kHz
tclean_params['restfreq'] = '%.6fGHz'%(line_restfreq_GHz / (1.0+galaxy_redshift)) # precision 1 kHz
tclean_params['width'] = 1
tclean_params['outframe'] = 'LSRK'
tclean_params['veltype'] = 'radio'
tclean_params['gridder'] = 'mosaic'
tclean_params['restoringbeam'] = '' # 'common'
tclean_params['weighting'] = ''
tclean_params['robust'] = 0.5
tclean_params['niter'] = 0
tclean_params['threshold'] = 0.0
tclean_params['usemask'] = 'pb'


# imaging dirty
image_name = os.path.join(run_tclean_dir, galaxy_name+'_'+line_name+'_dirty')
tclean_params['weighting'] = 'natural'
tclean_params['niter'] = 0
tclean_params['threshold'] = 0.0
if not os.path.exists(image_name+'.image.fits'):
    if os.path.exists(image_name+'.image.rms.txt'):
        os.remove(image_name+'.image.rms.txt')
    tclean_params['imagename'] = image_name
    cleanup_tclean_products(image_name)
    _print_params(tclean_params, 'tclean')
    tclean(**tclean_params)
    shutil.copy2('tclean.last', image_name+'.image'+'.tclean.last')
    export_tclean_products_as_fits_files(image_name)
    #apply_pbcor_to_tclean_image(image_name)

# estimate rms
if not os.path.exists(image_name+'.image.rms.txt'):
    imstat_results = imstat(image_name+'.image')
    rms = imstat_results['rms']
    with open(image_name+'.image.rms.txt', 'w') as fp:
        fp.write('%.10e\n'%(rms))
with open(image_name+'.image.rms.txt', 'r') as fp:
    rms = float(fp.readline())

# imaging natural
image_name = os.path.join(run_tclean_dir, galaxy_name+'_'+line_name+'_nat')
tclean_params['weighting'] = 'natural'
tclean_params['niter'] = 3000
tclean_params['threshold'] = 2.0 * rms
if not os.path.exists(image_name+'.image.fits'):
    if os.path.exists(image_name+'.residual.rms.txt'):
        os.remove(image_name+'.residual.rms.txt')
    tclean_params['imagename'] = image_name
    cleanup_tclean_products(image_name)
    _print_params(tclean_params, 'tclean')
    tclean(**tclean_params)
    shutil.copy2('tclean.last', image_name+'.image'+'.tclean.last')
    export_tclean_products_as_fits_files(image_name)
    apply_pbcor_to_tclean_image(image_name)

# estimate rms
if not os.path.exists(image_name+'.residual.rms.txt'):
    imstat_results = imstat(image_name+'.residual')
    rms = imstat_results['rms']
    with open(image_name+'.residual.rms.txt', 'w') as fp:
        fp.write('%.10e\n'%(rms))
else:
    with open(image_name+'.residual.rms.txt', 'r') as fp:
        rms = float(fp.readline())

# imaging briggs
image_name = os.path.join(run_tclean_dir, galaxy_name+'_'+line_name+'_briggs')
tclean_params['weighting'] = 'briggs'
tclean_params['robust'] = 0.5
tclean_params['niter'] = 3000
tclean_params['threshold'] = 2.0 * rms
if not os.path.exists(image_name+'.image.fits'):
    if os.path.exists(image_name+'.residual.rms.txt'):
        os.remove(image_name+'.residual.rms.txt')
    tclean_params['imagename'] = image_name
    cleanup_tclean_products(image_name)
    _print_params(tclean_params, 'tclean')
    tclean(**tclean_params)
    shutil.copy2('tclean.last', image_name+'.image'+'.tclean.last')
    export_tclean_products_as_fits_files(image_name)
    apply_pbcor_to_tclean_image(image_name)

# estimate rms
if not os.path.exists(image_name+'.residual.rms.txt'):
    imstat_results = imstat(image_name+'.residual')
    rms = imstat_results['rms']
    with open(image_name+'.residual.rms.txt', 'w') as fp:
        fp.write('%.10e\n'%(rms))
else:
    with open(image_name+'.residual.rms.txt', 'r') as fp:
        rms = float(fp.readline())

# imaging uniform
image_name = os.path.join(run_tclean_dir, galaxy_name+'_'+line_name+'_uni')
tclean_params['weighting'] = 'uniform'
tclean_params['niter'] = 3000
tclean_params['threshold'] = 3.0 * rms
if not os.path.exists(image_name+'.image.fits'):
    if os.path.exists(image_name+'.residual.rms.txt'):
        os.remove(image_name+'.residual.rms.txt')
    tclean_params['imagename'] = image_name
    cleanup_tclean_products(image_name)
    _print_params(tclean_params, 'tclean')
    tclean(**tclean_params)
    shutil.copy2('tclean.last', image_name+'.image'+'.tclean.last')
    export_tclean_products_as_fits_files(image_name)
    apply_pbcor_to_tclean_image(image_name)

# estimate rms
if not os.path.exists(image_name+'.residual.rms.txt'):
    imstat_results = imstat(image_name+'.residual')
    rms = imstat_results['rms']
    with open(image_name+'.residual.rms.txt', 'w') as fp:
        fp.write('%.10e\n'%(rms))
else:
    with open(image_name+'.residual.rms.txt', 'r') as fp:
        rms = float(fp.readline())

# imaging uniform
image_name = os.path.join(run_tclean_dir, galaxy_name+'_'+line_name+'_superuni')
tclean_params['weighting'] = 'superuniform'
tclean_params['niter'] = 3000
tclean_params['threshold'] = 3.0 * rms
if not os.path.exists(image_name+'.image.fits'):
    if os.path.exists(image_name+'.residual.rms.txt'):
        os.remove(image_name+'.residual.rms.txt')
    tclean_params['imagename'] = image_name
    cleanup_tclean_products(image_name)
    _print_params(tclean_params, 'tclean')
    tclean(**tclean_params)
    shutil.copy2('tclean.last', image_name+'.image'+'.tclean.last')
    export_tclean_products_as_fits_files(image_name)
    apply_pbcor_to_tclean_image(image_name)

# estimate rms
if not os.path.exists(image_name+'.residual.rms.txt'):
    imstat_results = imstat(image_name+'.residual')
    rms = imstat_results['rms']
    with open(image_name+'.residual.rms.txt', 'w') as fp:
        fp.write('%.10e\n'%(rms))
else:
    with open(image_name+'.residual.rms.txt', 'r') as fp:
        rms = float(fp.readline())

# all done
print('All done!')

