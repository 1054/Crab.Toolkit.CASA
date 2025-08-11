
# RUN THIS INSIDE CASA
# Last update: 2023-08-03
# By: Daizhong Liu

import os, sys, re, copy, glob, shutil, datetime
if not (os.path.expanduser('~/Cloud/Github/Crab.Toolkit.CASA/lib/python/analysis_scripts') in sys.path):
    sys.path.insert(1, os.path.expanduser('~/Cloud/Github/Crab.Toolkit.CASA/lib/python/analysis_scripts'))
import analysisUtils as au
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
    '2021.1.00024.S',
    '2017.1.00138.S', #mosaic
    '2021.1.00547.S', #mosaic
]

uvdata_dirs = [
    '/nfs/irdata02/datamine_radio/alma_archive/uvdata',
    #'/nfs/irdata07/dzliu/datamining/uvdata_by_a3cosmos',
]

mem_ous_ids = [ # these are only Band 3
    'uid://A001/X133d/X2026', 'uid://A001/X133d/X202a',
    'uid://A001/X2f6/X45e', 
	'uid://A001/X2fb/X818',
	'uid://A001/X1590/X35b7', 'uid://A001/X1590/X35b3',
    'uid://A001/X1288/X647', 'uid://A001/X1288/X643', 'uid://A001/X1288/X64b', 
    'uid://A001/X158f/X745', 'uid://A001/X158f/X749', 'uid://A001/X158f/X74d', 
]

def check_mem_ous_ids(file_list):
    out_file_list = []
    for file_name in file_list:
        matched = False
        for mem_ous_id in mem_ous_ids:
            mem_ous_id_str = re.sub(r'[^a-zA-Z0-9_]', r'_', mem_ous_id)
            if file_name.find(mem_ous_id_str) >= 0:
                matched = True
                break
        if matched:
            out_file_list.append(file_name)
    return out_file_list

# galaxy info
galaxy_name = 'K20-ID7'
galaxy_redshift = 2.224
product_name = 'cont100GHz'
phasecenter = 'J2000 53.131083deg -27.773108deg'
fov_arcsec = 60.0
contsub_lines = [
    dict(name='CO32', restfreq=345.7959899, linewidth=1000.0),
    #dict(name='CI21', restfreq=809.34197, linewidth=800.),
    #dict(name='CO76', restfreq=806.6518060, linewidth=1200.),
]

run_tclean_dir = 'run_tclean_for_continuum'

if not os.path.isdir(run_tclean_dir):
    os.makedirs(run_tclean_dir)


# loop projects
vis_list_to_combine = []
for project in list_of_projects:

    # output base name
    vis_to_combine = 'input_data_to_combine_for_continuum_'+project
    print('project', project)

    # find ms in an alma project
    list_of_files = []
    find_file_patterns = []
    for uvdata_dir in uvdata_dirs:
        if len(list_of_files) == 0:
            find_file = uvdata_dir+'/'+project+'/s*/g*/m*/calibrated/calibrated.ms'
            list_of_files = glob.glob(find_file)
            list_of_files = check_mem_ous_ids(list_of_files)
            find_file_patterns.append(find_file)
        if len(list_of_files) == 0:
            find_file = uvdata_dir+'/'+project+'/s*/g*/m*/calibrated/calibrated_*.ms'
            list_of_files = glob.glob(find_file)
            list_of_files = check_mem_ous_ids(list_of_files)
            find_file_patterns.append(find_file)
    if len(list_of_files) == 0:
        raise Exception('Error! Could not find "%s"'%(find_file))
    print('list_of_files', list_of_files)

    # loop each ms and process the uvdata
    if len(glob.glob(vis_to_combine+'*.ms')) != len(list_of_files):

        for ivis, vis in enumerate(list_of_files):
            
            # set outputvis name
            if len(list_of_files) != 1:
                vis_to_flag = vis_to_combine+'.%d.ms.to.flag'%(ivis+1)
                vis_to_concat = vis_to_combine+'.%d.ms'%(ivis+1)
            else:
                vis_to_flag = vis_to_combine+'.ms.to.flag'
                vis_to_concat = vis_to_combine+'.ms'

            # check outputvis
            #if os.path.isdir(vis_to_concat):
            #    shutil.rmtree(vis_to_concat)
            if not os.path.isdir(vis_to_concat):
    
                # find field
                print('get_field_IDs_in_mosaic', vis)
                field_ids = get_field_IDs_in_mosaic(\
                    vis = vis,
                    cell = '0.1arcsec', imsize = [50, 50], # representing a 5 arcsec size searching box
                    phasecenter = phasecenter,
                    padding_by_primary_beam = 0.0,
                    verbose = True,
                )

                print('field_ids', field_ids)

                # find science spw
                science_spws = au.getScienceSpws(vis)

                print('science_spws', science_spws)

                # split
                #if os.path.isdir(vis_to_flag):
                #    shutil.rmtree(vis_to_flag)
                if not os.path.isdir(vis_to_flag):
                    split_params = {}
                    split_params['vis'] = vis
                    split_params['outputvis'] = vis_to_flag
                    split_params['field'] = ','.join([str(t) for t in field_ids])
                    split_params['spw'] = science_spws
                    split_params['timebin'] = '10s'
                    split_params['keepflags'] = False
                    split_params['datacolumn'] = get_datacolumn(vis)
                    _print_params(split_params, 'split')
                    split(**split_params)
                    shutil.copy2('split.last', vis_to_flag+'.split.last')

                # flag line channels
                if not os.path.isdir(vis_to_concat):
                    spw_selection_str = ''
                    if len(contsub_lines) > 0:
                        tb.open(vis_to_flag+'/SPECTRAL_WINDOW')
                        spw_names = tb.getcol('NAME')
                        spw_chan_freqs = [tb.getcell('CHAN_FREQ', ispw) for ispw in range(len(spw_names))]
                        tb.close()
                        # loop each spw
                        for ispw in range(len(spw_names)):
                            # loop lines to flag
                            freqarr = np.array(spw_chan_freqs)
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
                            print('ispw %d, flagging lines with fitspw = "%s"'%(ispw, fitspw))
                            # append to spw_selection_str
                            if ispw > 0 and spw_selection_str != '':
                                spw_selection_str += ','
                            spw_selection_str += '%d:%s'(ispw, fitspw)
                    print('spw_selection_str', spw_selection_str)
                    #raise NotImplementedError()
                    # 
                    # do flagging
                    if spw_selection_str != '':
                        flagdata_params = {}
                        flagdata_params['vis'] = vis_to_flag
                        flagdata_params['mode'] = 'manual'
                        flagdata_params['field'] = ''
                        flagdata_params['spw'] = spw_selection_str
                        _print_params(flagdata_params, 'flagdata')
                        flagdata(**flagdata_params)
                        shutil.copy2('flagdata.last', vis_to_flag+'.flagdata.last')
                    # 
                    # rename
                    shutil.copytree(vis_to_flag, vis_to_concat)

            vis_list_to_combine.append(vis_to_concat)

    else:

        vis_list_to_combine.extend(glob.glob(vis_to_combine+'*.ms'))

print('vis_list_to_combine', vis_list_to_combine)


# concat
vis_combined = 'combined_for_continuum.ms'
if not os.path.exists(vis_combined):
    concat_params = {'vis': vis_list_to_combine, 'concatvis': vis_combined, 'freqtol':'0Hz'}
    _print_params(concat_params, 'concat')
    concat(**concat_params)
    shutil.copy2('concat.last', vis_combined+'.concat.last')


# statwt
if not os.path.exists(vis_combined+'.statwt.done'):
    statwt_params = {'vis': vis_combined, 'spw': '', 'datacolumn': get_datacolumn(vis_combined)}
    _print_params(statwt_params, 'statwt')
    statwt(**statwt_params)
    shutil.copy2('statwt.last', vis_combined+'.statwt.last')
    with open(vis_combined+'.statwt.done', 'w') as fp:
        fp.write('done at %s\n'%(datetime.datetime.now()))


vis_combined_for_continuum = 'combined_for_continuum.ms'


# imaging prep
synbeam, imcell = get_synbeam_and_imcell(vis_combined_for_continuum)
imsize = int(np.ceil(fov_arcsec/float(imcell.replace('arcsec',''))))
tclean_params = OrderedDict()
tclean_params['vis'] = vis_combined_for_continuum
tclean_params['imagename'] = ''
tclean_params['datacolumn'] = get_datacolumn(vis_combined_for_continuum)
tclean_params['field'] = ''
tclean_params['spw'] = ''
tclean_params['imsize'] = imsize
tclean_params['cell'] = imcell
tclean_params['phasecenter'] = phasecenter
tclean_params['specmode'] = 'mfs'
#tclean_params['reffreq'] = '%.6fGHz'%(line_restfreq_GHz / (1.0+galaxy_redshift)) # precision 1 kHz
#tclean_params['restfreq'] = '%.6fGHz'%(line_restfreq_GHz / (1.0+galaxy_redshift)) # precision 1 kHz
#tclean_params['width'] = 1
tclean_params['outframe'] = 'LSRK'
tclean_params['veltype'] = 'radio'
tclean_params['gridder'] = 'mosaic'
tclean_params['restoringbeam'] = '' # per chan
tclean_params['weighting'] = ''
tclean_params['robust'] = 0.5
tclean_params['niter'] = 0
tclean_params['threshold'] = 0.0
tclean_params['usemask'] = 'pb'

if not os.path.isdir(run_tclean_dir):
    os.makedirs(run_tclean_dir)

# imaging dirty
image_name = os.path.join(run_tclean_dir, galaxy_name+'_'+product_name+'_dirty')
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
image_name = os.path.join(run_tclean_dir, galaxy_name+'_'+product_name+'_nat')
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
with open(image_name+'.residual.rms.txt', 'r') as fp:
    rms = float(fp.readline())

# imaging briggs
image_name = os.path.join(run_tclean_dir, galaxy_name+'_'+product_name+'_briggs')
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
with open(image_name+'.residual.rms.txt', 'r') as fp:
    rms = float(fp.readline())

# imaging uniform
image_name = os.path.join(run_tclean_dir, galaxy_name+'_'+product_name+'_uni')
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
with open(image_name+'.residual.rms.txt', 'r') as fp:
    rms = float(fp.readline())

# imaging uniform
image_name = os.path.join(run_tclean_dir, galaxy_name+'_'+product_name+'_superuni')
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
with open(image_name+'.residual.rms.txt', 'r') as fp:
    rms = float(fp.readline())

# all done
print('All done!')

