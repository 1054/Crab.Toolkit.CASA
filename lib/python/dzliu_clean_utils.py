# 
# These functions must be run inside CASA.
# 
# This script contains following functions (not a complete list):
#     get_datacolumn
#     get_antenn_diameter
#     get_ref_frequency
#     get_optimized_imsize
#     get_field_phasecenters
#     get_mosaic_imsize_and_phasecenter
#     cleanup_tclean_products
#     apply_pbcor_to_tclean_image
#     export_tclean_products_as_fits_files
#     imsmooth_tclean_image
# 
# Last update: 
#     2020-12-23 copied functions from "dzliu_clean.py"
#     2021-01-04 updated tclean functions
# 

import os, sys, re, copy, shutil
import numpy as np

from taskinit import casalog, tb #, ms, iatool
#from taskinit import casac
#tb = casac.table
#from __casac__.table import table as tb
#from recipes import makepb, pixelmask2cleanmask
import casadef
def version_tuple(version_str):
    return tuple(map(int, (version_str.split("."))))
def version_less_than(version_str, compared_version_str):
    return version_tuple(version_str) < version_tuple(compared_version_str)
def version_greater_equal(version_str, compared_version_str):
    return version_tuple(version_str) >= version_tuple(compared_version_str)
if version_less_than(casadef.casa_version, '6.0.0'):
    #from __main__ import default, inp, saveinputs
    ##import task_tclean; task_tclean.tclean # this is what tclean.py calls
    ##import tclean
    ##import tclean_cli
    from tclean_cli import tclean_cli_
    tclean = tclean_cli_()
    from mstransform_cli import mstransform_cli_
    mstransform = mstransform_cli_()
    from exportfits_cli import exportfits_cli_
    exportfits = exportfits_cli_()
    from concat_cli import concat_cli_
    concat = concat_cli_()
    from split_cli import split_cli_
    split = split_cli_()
    from imstat_cli import imstat_cli_
    imstat = imstat_cli_()
    from impbcor_cli import impbcor_cli_
    impbcor = impbcor_cli_()
    from imsmooth_cli import imsmooth_cli_
    imsmooth = imsmooth_cli_()
else:
    # see CASA 6 updates here: https://alma-intweb.mtk.nao.ac.jp/~eaarc/UM2018/presentation/Nakazato.pdf
    from casatasks import tclean, mstransform, exportfits, concat, split, imstat, impbcor, imsmooth
    #from casatasks import sdbaseline
    #from casatools import ia



# 
# def print2
# 
def print2(message):
    print(message)
    casalog.post(message, 'INFO')



# 
# get datacolumn
# 
def get_datacolumn(vis):
    #
    # Requires CASA module/function tb.
    #
    casalog.origin('get_datacolumn')
    #
    tb.open(vis)
    if 'CORRECTED_DATA' in tb.colnames():
        datacolumn = 'CORRECTED'
    else:
        datacolumn = 'DATA'
    tb.close()
    #
    return datacolumn



#
# def get antenna diameter
#
def get_antenn_diameter(vis):
    #
    # Requires CASA module/function tb.
    #
    casalog.origin('get_antenn_diameter')
    #
    tb.open(vis+os.sep+'ANTENNA')
    ant_names = tb.getcol('NAME')
    ant_diams = tb.getcol('DISH_DIAMETER') # meter
    tb.close()
    #
    minantdiam = np.min(ant_diams) # meter
    print2('ant_diams = %s'%(ant_diams))
    print2('minantdiam = %s [m]'%(minantdiam))
    #
    return minantdiam



#
# def get_ref_frequency
#
def get_ref_frequency(vis, spw_list=None):
    return get_ref_frequency_Hz(vis, spw_list=spw_list)
#
def get_ref_frequency_Hz(vis, spw_list=None):
    #
    # Requires CASA module/function tb.
    #
    casalog.origin('get_central_channel_frequency')
    #
    tb.open(vis+os.sep+'SPECTRAL_WINDOW')
    #spw_chan_freq = tb.getcol('CHAN_FREQ') # a list of list, Hz
    spw_ref_frequency = tb.getcol('REF_FREQUENCY') # Hz
    tb.close()
    #
    if spw_list is not None:
        ref_frequency = spw_ref_frequency[spw_list]
    else:
        ref_frequency = spw_ref_frequency[0] # if no spw is specified, use the first spw REF_FREQUNCY
    # 
    return ref_frequency



#
# def get_chan_width_MHz
#
def get_chan_width(vis, spw_list=None):
    return get_chan_width_MHz(vis, spw_list=spw_list)
# 
def get_chan_width_MHz(vis, spw_list=None):
    #
    # Requires CASA module/function tb.
    #
    casalog.origin('get_central_channel_frequency')
    #
    tb.open(vis+os.sep+'SPECTRAL_WINDOW')
    spw_id = np.arange(0,tb.nrows()) # 
    spw_chan_width = np.array([tb.getcell('CHAN_WIDTH', i)[0] for i in spw_id]) # a list of list, Hz
    tb.close()
    #
    if spw_list is not None:
        chan_width_array_MHz = np.array([spw_chan_width[t] for t in spw_list]) / 1e6
        return chan_width_array_MHz
    else:
        chan_width_MHz = spw_chan_width[0] / 1e6 # if no spw is specified, use the first spw CHAN_WIDTH
        return chan_width_MHz
    # 
    #return chan_width_MHz



# 
# def get_optimized_imsize which is a factor of 2,3,5,7
# 
def get_optimized_imsize(imsize, return_decomposed_factors = False):
    # 
    # No CASA module/function used here.
    # 
    #casalog.origin('get_optimized_imsize')
    # 
    # try to make imsize be even and only factorizable by 2,3,5,7
    imsize = int(imsize)
    decomposed_factors = []
    # 
    # if imsize is 1, then return it
    if imsize == 1:
        if return_decomposed_factors == True:
            return 1, [1]
        else:
            return 1
    # 
    # make it even
    if imsize % 2 != 0:
        imsize += 1
    # 
    # factorize by 2,3,5,7
    for k in [2, 3, 5, 7]:
        while imsize % k == 0:
            imsize = int(imsize / k)
            decomposed_factors.append(k)
    # 
    # make the non-factorizable number factorizable by 2, 3, 5, or 7
    while imsize != 1 and int( np.prod( [ (imsize % k) for k in [2, 3, 5, 7] ] ) ) != 0:
        # as long as it is factorizable by any of the [2, 3, 5, 7], the mod ("%") will be zero, so the product will also be zero
        #print('imsize', imsize, '(imsize % k)', [ (imsize % k) for k in [2, 3, 5, 7] ], 
        #                               np.prod( [ (imsize % k) for k in [2, 3, 5, 7] ] ) )
        imsize += 1
        #print('imsize', imsize, '(imsize % k)', [ (imsize % k) for k in [2, 3, 5, 7] ], 
        #                               np.prod( [ (imsize % k) for k in [2, 3, 5, 7] ] ) )
        
    # 
    imsize2, decomposed_factors2 = get_optimized_imsize(imsize, return_decomposed_factors = True)
    # 
    imsize = imsize2
    # 
    decomposed_factors.extend(decomposed_factors2)
    # 
    if return_decomposed_factors == True:
        return np.prod(decomposed_factors), decomposed_factors
    else:
        return np.prod(decomposed_factors)



#
# def get field phasecenters
#
def get_field_phasecenters(vis, galaxy_name = '', column_name = 'DELAY_DIR'):
    """
    Get Measurement Set phase centers ('DELAY_DIR'). 
    
    Return 3 lists: matched_field_name, matched_field_indices, and matched_field_phasecenters. 
    
    The 3rd return is a list of two lists: a RA_deg and a Dec_deg list.
    
    If galaxy_name is '', then all field phase centers will be returned.
    """
    #
    # Requires CASA module/function tb.
    #
    casalog.origin('get_field_phasecenters')
    #
    tb.open(vis+os.sep+'FIELD')
    field_names = tb.getcol('NAME')
    field_phasecenters = [tb.getcell(column_name, i) for i in range(tb.nrows())] # rad,rad
    tb.close()
    #
    if galaxy_name != '':
        galaxy_name_cleaned = re.sub(r'[^a-zA-Z0-9]', r'', galaxy_name).lower() #<TODO># What if someone use "_" as a field name?
    else:
        galaxy_name_cleaned = '' # if the user has input an empty string, then we will get all fields in this vis data.
    #
    matched_field_name = ''
    matched_field_indices = []
    matched_field_phasecenters = []
    for i, field_name in enumerate(field_names):
        # find galaxy_name in field_names:
        field_name_cleaned = re.sub(r'[^a-zA-Z0-9]', r'', field_name).lower()
        field_RA_rad, field_Dec_rad = field_phasecenters[i]
        field_RA_rad = field_RA_rad[0]
        field_Dec_rad = field_Dec_rad[0]
        if field_RA_rad < 0:
            field_RA_rad += 2.0 * np.pi
        field_RA_deg = field_RA_rad / np.pi * 180.0
        field_Dec_deg = field_Dec_rad / np.pi * 180.0
        if galaxy_name_cleaned == '' or field_name_cleaned.startswith(galaxy_name_cleaned):
            matched_field_name = field_name
            matched_field_indices.append(i)
            matched_field_phasecenters.append([field_RA_deg, field_Dec_deg])
    #
    if '' == matched_field_name:
        raise ValueError('Error! Target source %s was not found in the "FIELD" table of the input vis "%s"!'%(galaxy_name, vis))
    #
    matched_field_indices = np.array(matched_field_indices)
    matched_field_phasecenters = np.array(matched_field_phasecenters).T # two columns, nrows
    return matched_field_name, matched_field_indices, matched_field_phasecenters




#
# def get mosaic width and height in degree
#
def get_mosaic_imsize_and_phasecenter(vis, cell, galaxy_name='', ref_freq_Hz=None, padding_by_primary_beam=0.5, verbose=True):
    """
    cell is the same as pixel_size, and overrides pixel_size. 
    pixel_size can be a string or a float number. If it is a float number, a unit of arcsec is assumed.
    """
    # 
    # get antdiam
    minantdiam = get_antenn_diameter(vis)
    # 
    # get field and phase centers
    matched_field_name, matched_field_indices, matched_field_phasecenters = get_field_phasecenters(vis, galaxy_name=galaxy_name)
    matched_field_min_RA_deg = np.min(matched_field_phasecenters[0, :])
    matched_field_max_RA_deg = np.max(matched_field_phasecenters[0, :])
    matched_field_min_Dec_deg = np.min(matched_field_phasecenters[1, :])
    matched_field_max_Dec_deg = np.max(matched_field_phasecenters[1, :])
    matched_field_center_RA_deg = np.mean(matched_field_phasecenters[0, :])
    matched_field_center_Dec_deg = np.mean(matched_field_phasecenters[1, :])
    phasecenter = 'J2000 %.8fdeg %.8fdeg'%(matched_field_center_RA_deg, matched_field_center_Dec_deg)
    #print('matched_field_phasecenters.shape', matched_field_phasecenters.shape)
    #print('matched_field_phasecenters:', matched_field_phasecenters)
    #raise NotImplementedError()
    # 
    # calc primary beam
    if ref_freq_Hz is None:
        ref_freq_Hz = get_ref_frequency(vis)
    pribeam = 1.13  * (2.99792458e8 / ref_freq_Hz / minantdiam / np.pi * 180.0 ) # in units of degrees, see -- https://help.almascience.org/index.php?/Knowledgebase/Article/View/90
    if verbose:
        print2('minantdiam = %s [meter]'%(minantdiam))
        print2('pribeam = %s [arcsec]'%(pribeam * 3600.0))
        print2('matched_field_phasecenters = %s'%(matched_field_phasecenters))
        print2('matched_field_min_RA_deg = %s'%(matched_field_min_RA_deg))
        print2('matched_field_max_RA_deg = %s'%(matched_field_max_RA_deg))
        print2('matched_field_min_Dec_deg = %s'%(matched_field_min_Dec_deg))
        print2('matched_field_max_Dec_deg = %s'%(matched_field_max_Dec_deg))
    # 
    # calc mosaic width and height, half primary beam padding at both sides are considered.
    imsize_RA_deg = (matched_field_max_RA_deg - matched_field_min_RA_deg) * np.cos(np.deg2rad((matched_field_max_Dec_deg+matched_field_min_Dec_deg)/2.0))
    imsize_Dec_deg = (matched_field_max_Dec_deg - matched_field_min_Dec_deg)
    imsize_RA_deg = imsize_RA_deg + 2.0 * padding_by_primary_beam * pribeam # padding this size at each side
    imsize_Dec_deg = imsize_Dec_deg + 2.0 * padding_by_primary_beam * pribeam # padding this size at each side
    if verbose:
        print2('imsize_RA = %s [arcsec]'%(imsize_RA_deg * 3600.0))
        print2('imsize_Dec = %s [arcsec]'%(imsize_Dec_deg * 3600.0))
    # 
    # get imcell_arcsec (pixel_size) from the given cell
    pixel_size = re.sub(r'[^0-9.a-zA-Z+-]', r'', str(cell))
    if re.match(r'^([0-9.eE+-]+)(asec|arcsec)$', pixel_size): 
        imcell_arcsec = float(re.sub(r'^([0-9.eE+-]+)(asec|arcsec)$', r'\1', pixel_size))
    elif re.match(r'^([0-9.eE+-]+)(amin|arcmin)$', pixel_size): 
        imcell_arcsec = float(re.sub(r'^([0-9.eE+-]+)(amin|arcmin)$', r'\1', pixel_size)) * 60.
    elif re.match(r'^([0-9.eE+-]+)(deg|degree)$', pixel_size): 
        imcell_arcsec = float(re.sub(r'^([0-9.eE+-]+)(deg|degree)$', r'\1', pixel_size)) * 3600.
    elif re.match(r'^([0-9.eE+-]+)$', pixel_size): 
        imcell_arcsec = float(re.sub(r'^([0-9.eE+-]+)$', r'\1', pixel_size)) # in default we assume arcsec unit
    else:
        print2('Error! The input pixel_size could not be understood. It should be a string with a unit, e.g. \'1.0arcsec\'.')
        raise Exception('Error! The input pixel_size could not be understood.')
    if verbose:
        print2('imcell = %s [arcsec]'%(imcell_arcsec))
    # 
    # 
    imsize_RA = imsize_RA_deg / (imcell_arcsec / 3600.0) # pixels
    imsize_Dec = imsize_Dec_deg / (imcell_arcsec / 3600.0) # pixels
    imsize = [get_optimized_imsize(imsize_RA), get_optimized_imsize(imsize_Dec)]
    if verbose:
        print2('imsize_RA = %s [pixel]'%(imsize_RA))
        print2('imsize_Dec = %s [pixel]'%(imsize_Dec))
    return imsize, phasecenter



#
# def get spw for spectral line
#
def get_spw_for_spectral_line(vis, redshift=None, rest_freq_GHz=None, line_width_kms=None, return_dict=False, verbose=True):
    #
    casalog.origin('get_spw_for_spectral_line')
    # 
    if redshift is None or rest_freq_GHz is None:
        print2('Error! Please input redshift and rest_freq_GHz!')
        raise Exception('Error! Please input redshift and rest_freq_GHz!')
    # 
    rest_freq_Hz = rest_freq_GHz * 1e9
    line_freq_GHz = rest_freq_GHz / (1.0+redshift)
    if line_width_kms is None:
        print2('Setting line_width_kms to default value 1000.0 km/s.')
        line_width_kms = 1000.0 # km/s
    line_width_MHz = line_width_kms/2.99792458e5*line_freq_GHz*1000. # km/s
    line_freq_range_Hz = [line_freq_GHz*1e9 - line_width_MHz/2.0*1e6, line_freq_GHz*1e9 + line_width_MHz/2.0*1e6]
    if verbose:
        print2('line_width_MHz = %s'%(line_width_MHz))
        print2('line_freq_range_Hz = %s'%(line_freq_range_Hz))
    # 
    #au.getScienceSpws(vis)
    #
    tb.open(vis+os.sep+'SPECTRAL_WINDOW')
    spw_id = np.arange(0,tb.nrows()) # 
    spw_name = tb.getcol('NAME') # 
    spw_nchan = tb.getcol('NUM_CHAN') # 
    spw_chan_freq = np.array([tb.getcell('CHAN_FREQ', i)[0] for i in spw_id]) # a list of list, Hz
    spw_chan_width = np.array([tb.getcell('CHAN_WIDTH', i)[0] for i in spw_id]) # a list of list, Hz, assuming all channels have the same width in a spw
    spw_ref_freq = tb.getcol('REF_FREQUENCY') # Hz
    tb.close()
    # 
    spw_selection_dict = {}
    spw_selection_str = ''
    for i in spw_id:
        nchan = spw_nchan[i]
        ch0 = spw_chan_freq[i] # 
        chstep = spw_chan_width[i] # 
        chlast = ch0 + (nchan-1.) * chstep
        if verbose:
            print2('spw %s, ch0 %s, chlast %s, chstep %s, nchan %s'%(i, ch0, chlast, chstep, nchan))
        if (line_freq_range_Hz[1] > min(ch0, chlast)) and (line_freq_range_Hz[0] < max(ch0, chlast)) and nchan > 1:
            if chstep > 0:
                # line in spw
                chleft = np.ceil((line_freq_range_Hz[0] - ch0) / chstep)
                chright = np.floor((line_freq_range_Hz[1] - ch0) / chstep)
            else:
                chleft = np.ceil((line_freq_range_Hz[1] - ch0) / chstep)
                chright = np.floor((line_freq_range_Hz[0] - ch0) / chstep)
            # 
            if chleft < 0:
                chleft = 0
            if chright > nchan-1:
                chright = nchan-1
            # 
            spw_selection_dict[i] = '%d~%d'%(chleft, chright)
            if spw_selection_str != '':
                spw_selection_str += ','
            spw_selection_str += '%d:%d~%d'%(i, chleft, chright)
    # 
    if verbose:
        print2('spw_selection_str = %r'%(spw_selection_str))
    # 
    if return_dict:
        return spw_selection_str, spw_selection_dict
    else:
        return spw_selection_str



#
# def get_mstransform_params_for_spectral_line
#
def get_mstransform_params_for_spectral_line(
        vis, 
        outputvis, 
        field=None, 
        redshift=None, 
        rest_freq_GHz=None, 
        line_width_kms=None, 
        chan_width_kms=None, 
        force_integer_chan_width=True, 
        verbose=True, 
    ):
    #
    casalog.origin('get_spw_for_spectral_line')
    # 
    if field is None or redshift is None or rest_freq_GHz is None or line_width_kms is None or chan_width_kms is None:
        print2('Error! Please input field, redshift, rest_freq_GHz, line_width_kms and chan_width_kms!')
        raise Exception('Error! Please input field, redshift, rest_freq_GHz, line_width_kms and chan_width_kms!')
    # 
    rest_freq_Hz = rest_freq_GHz * 1e9
    line_freq_GHz = rest_freq_GHz / (1.0+redshift)
    line_freq_Hz = line_freq_GHz * 1e9
    line_freq_MHz = line_freq_GHz * 1e3
    line_width_MHz = line_width_kms/2.99792458e5*line_freq_GHz*1000. # km/s -> MHz
    line_freq_range_Hz = [line_freq_GHz*1e9 - line_width_MHz/2.0*1e6, line_freq_GHz*1e9 + line_width_MHz/2.0*1e6]
    # 
    #au.getScienceSpws(vis)
    #
    tb.open(vis+os.sep+'SPECTRAL_WINDOW')
    spw_id = np.arange(0,tb.nrows()) # 
    spw_name = tb.getcol('NAME') # 
    spw_nchan = tb.getcol('NUM_CHAN') # 
    spw_chan_freq = np.array([tb.getcell('CHAN_FREQ', i)[0] for i in spw_id]) # a list of list, Hz
    spw_chan_width = np.array([tb.getcell('CHAN_WIDTH', i)[0] for i in spw_id]) # a list of list, Hz, assuming all channels have the same width in a spw
    spw_ref_freq = tb.getcol('REF_FREQUENCY') # Hz
    tb.close()
    # 
    mstransform_params = {}
    # 
    spw_selection_str, spw_selection_dict = get_spw_for_spectral_line(vis, redshift=redshift, rest_freq_GHz=rest_freq_GHz, line_width_kms=line_width_kms, return_dict=True)
    # 
    output_chan_width_MHz = (chan_width_kms/2.99792458e5*line_freq_MHz)
    if force_integer_chan_width:
        chan_width_MHz = get_chan_width_MHz(vis, spw_list=list(spw_selection_dict.keys()))
        if verbose:
            print2('chan_width_MHz: %s, spw_selection_dict: %s'%(chan_width_MHz, spw_selection_dict))
        chan_width_MHz = np.abs(chan_width_MHz)
        min_chan_width_MHz = chan_width_MHz[0]
        for i,ispw in enumerate(list(spw_selection_dict.keys())):
            if chan_width_MHz[i] > output_chan_width_MHz:
                if verbose:
                    print2('Discarding spw %d because its channel width %s MHz is broader than the output channel width %s MHz'%(ispw, chan_width_MHz[i], output_chan_width_MHz))
                del spw_selection_dict[ispw]
            else:
                if min_chan_width_MHz > chan_width_MHz[i]:
                    min_chan_width_MHz = chan_width_MHz[i]
        output_chan_width_MHz = np.round(output_chan_width_MHz/min_chan_width_MHz)*min_chan_width_MHz
    # 
    output_nchan = int(np.ceil(line_width_MHz / output_chan_width_MHz))
    # 
    mstransform_params['vis'] = vis
    mstransform_params['outputvis'] = outputvis
    mstransform_params['field'] = field
    mstransform_params['spw'] = ','.join([str(t) for t in spw_selection_dict.keys()])
    mstransform_params['width'] = '%.6fMHz'%(output_chan_width_MHz)
    mstransform_params['regridms'] = True
    mstransform_params['mode'] = 'frequency'
    mstransform_params['restfreq'] = '%.6fMHz'%(line_freq_MHz)
    mstransform_params['start'] = '%.6fMHz'%(line_freq_MHz - output_chan_width_MHz*(output_nchan-1.)/2.)
    mstransform_params['nchan'] = output_nchan
    mstransform_params['outframe'] = 'LSRK'
    mstransform_params['datacolumn'] = get_datacolumn(vis)
    mstransform_params['combinespws'] = True
    mstransform_params['keepflags'] = False
    #mstransform_params['keepmms'] = False
    # 
    if verbose:
        print2('mstransform_params: mstransform('+', '.join("{!s}={!r}".format(k, mstransform_params[k]) for k in mstransform_params.keys())+')')
    # 
    return mstransform_params



# 
# 
# 
def cleanup_tclean_products(imagename, suffix_list=None, cleanup_mask=True, cleanup_fits=True, exit_on_error=True):
    if imagename.endswith('.image'):
        imagename = re.sub(r'\.image$', r'', imagename)
    if suffix_list is None:
        suffix_list = ['.image', '.image.pbcor', '.model', '.pb', '.psf', '.residual', '.weight', '.sumwt', 
                       '.mask', 
                       '.image.tt0', 'image.pbcor.tt0', '.model.tt0', '.pb.tt0', '.psf.tt0', '.residual.tt0', '.weight.tt0', '.sumwt.tt0', 
                       '.image.tt1', 'image.pbcor.tt1', '.model.tt1', '.pb.tt1', '.psf.tt1', '.residual.tt1', '.weight.tt1', '.sumwt.tt1', 
                       '.image.tt2', 'image.pbcor.tt2', '.model.tt2', '.pb.tt2', '.psf.tt2', '.residual.tt2', '.weight.tt2', '.sumwt.tt2', 
                       '.alpha', '.alpha.error'] #<TODO># depends on CASA version and tclean cube type
    if cleanup_mask:
        suffix_list.append('.mask')
    for suffix in suffix_list:
        if os.path.isdir(imagename+suffix):
            shutil.rmtree(imagename+suffix)
            if not os.path.isdir(imagename+suffix):
                print2('Deleted "%s"'%(imagename+suffix))
            else:
                if exit_on_error:
                    raise Exception('Error! Failed to cleanup tclean product data directory %s'%(imagename+suffix))
        if cleanup_fits:
            if os.path.isfile(imagename+suffix+'.fits'):
                os.remove(imagename+suffix+'.fits')
                if not os.path.isfile(imagename+suffix+'.fits'):
                    print2('Deleted "%s"'%(imagename+suffix+'.fits'))
                else:
                    if exit_on_error:
                        raise Exception('Error! Failed to cleanup tclean product data file %s'%(imagename+suffix+'.fits'))


def apply_pbcor_to_tclean_image(imagename, cutoff=0.1, overwrite=True, exit_on_error=True):
    if imagename.endswith('.image'):
        imagename = re.sub(r'\.image$', r'', imagename)
    elif imagename.endswith('.image.tt0'):
        imagename = re.sub(r'\.image\.tt0$', r'', imagename)
    # 
    for suffix in ['', '.tt0']:
        infile2 = imagename+'.image'+suffix
        pbimage2 = imagename+'.pb'+suffix
        outfile2 = imagename+'.image.pbcor'+suffix
        #if not os.path.isdir(infile2) or not os.path.isdir(pbimage2):
        #    raise Exception('Error! Data not found: "%s" or "%s". Please check your input imagename and pbimage'%(infile2, pbimage2))
        if os.path.isdir(infile2) and os.path.isdir(pbimage2):
            if os.path.isdir(outfile2):
                if not overwrite:
                    raise Exception('Found existing data "%s"! Please clean it up first!'%(outfile2))
                else:
                    print2('Found existing data "%s", overwriting it.'%(outfile2))
                    shutil.rmtree(outfile2)
            print2('Running CASA task: impbcor(imagename=%r, pbimage=%r, outfile=%r, mode=%r, cutoff=%s)'%(infile2, pbimage2, outfile2, 'divide', cutoff))
            impbcor(imagename=infile2, pbimage=pbimage2, outfile=outfile2, mode='divide', cutoff=cutoff)
            if os.path.isdir(outfile2):
                print2('Output to "%s"'%(outfile2))
            else:
                if exit_on_error:
                    raise Exception('Error! Failed to run CASA impbcor and output "%s"'%(outfile2))
            # 
            export_tclean_products_as_fits_files(imagename, suffix_list=['.image.pbcor'+suffix])


def export_tclean_products_as_fits_files(imagename, dropstokes=True, suffix_list=None, overwrite=True, exit_on_error=True):
    if imagename.endswith('.image'):
        imagename = re.sub(r'\.image$', r'', imagename)
    elif imagename.endswith('.image.tt0'):
        imagename = re.sub(r'\.image\.tt0$', r'', imagename)
    if suffix_list is None:
        suffix_list = ['.image', '.image.pbcor', '.model', '.pb', '.psf', '.residual', '.weight', 
                       '.mask', 
                       '.image.tt0', 'image.pbcor.tt0', '.model.tt0', '.pb.tt0', '.psf.tt0', '.residual.tt0', '.weight.tt0', 
                       '.image.tt1', 'image.pbcor.tt1', '.model.tt1', '.pb.tt1', '.psf.tt1', '.residual.tt1', '.weight.tt1', 
                       '.image.tt2', 'image.pbcor.tt2', '.model.tt2', '.pb.tt2', '.psf.tt2', '.residual.tt2', '.weight.tt2', 
                       '.alpha', '.alpha.error'] #<TODO># depends on CASA version and tclean cube type
    for suffix in suffix_list:
        infile = imagename+suffix
        outfile = imagename+suffix+'.fits'
        if os.path.isdir(infile):
            if os.path.isfile(outfile):
                if not overwrite:
                    if exit_on_error:
                        raise Exception('Found existing tclean product fits file "%s"! Please clean it up first!'%(outfile))
                    else:
                        print2('Found existing tclean product fits file "%s", will not overwrite it.'%(outfile))
                        continue
                else:
                    print2('Found existing tclean product fits file "%s", overwriting it.'%(outfile))
                    os.remove(outfile)
            print2('Running CASA task: exportfits(%r, %r, dropstokes=%s)'%(infile, outfile, dropstokes))
            exportfits(infile, outfile, dropstokes=dropstokes)
            if os.path.isfile(outfile):
                print2('Output to "%s"'%(outfile))
            else:
                raise Exception('Error! Failed to run CASA exportfits and output %s'%(imagename+suffix))
    # 


def imsmooth_tclean_image(infile, major, minor=None, pa=None, kernel='gaussian', targetres=True, outfile=None, 
                          overwrite=True, exit_on_error=True, export_fits=True):
    # check user input
    if not isinstance(infile, str):
        raise Exception('Error! Please input the tclean image data file path. It must be an existing directory.')
    if not isinstance(major, str):
        raise Exception('Error! Please input the beam major axis FWHM size as a string like \'10.00arcsec\'. It will be used as the output name if outfile is not specified.')
    if minor is None:
        minor = major
    if pa is None:
        pa = '0deg'
    # check input data
    if not os.path.isdir(infile):
        raise Exception('Error! Input image data not found "%s"!'%(infile))
    # set default outfile name
    if outfile is None:
        outfile = infile+'.conv.to.beam.'+major+'.image' #<TODO># set default outfile name here
    # check output data
    run_imsmooth = True
    if os.path.isdir(outfile):
        if not overwrite:
            if exit_on_error:
                raise Exception('Found existing tclean product fits file "%s"! Please clean it up first!'%(outfile))
            else:
                print2('Found existing tclean product fits file "%s", will not overwrite it.'%(outfile))
                run_imsmooth = False
        else:
            print2('Found existing tclean product fits file "%s", overwriting it.'%(outfile))
            shutil.rmtree(outfile)
    # run imsmooth
    if run_imsmooth:
        print2('Running CASA task: imsmooth(imagename = %r, kernel = %r, major = %r, minor = %r, pa = %r, targetres = %s, outfile = %r)'%(infile, kernel, major, minor, pa, targetres, outfile))
        imsmooth(imagename = infile, kernel = kernel, major = major, minor = minor, pa = pa, targetres = targetres, outfile = outfile)
    # check result
    if os.path.isdir(outfile):
        print2('Output to "%s"'%(outfile))
    else:
        if exit_on_error:
            raise Exception('Error! Failed to run CASA imsmooth and output %s'%(outfile))
    # check output data fits file
    run_export_fits = True
    if os.path.isfile(outfile+'.fits'):
        if not overwrite:
            if exit_on_error:
                raise Exception('Found existing tclean product fits file "%s"! Please clean it up first!'%(outfile))
            else:
                print2('Found existing tclean product fits file "%s", will not overwrite it.'%(outfile))
                run_export_fits = False
        else:
            print2('Found existing tclean product fits file "%s", overwriting it.'%(outfile))
            os.remove(outfile+'.fits')
    # export fits file
    if export_fits and run_export_fits:
        exportfits(outfile, outfile+'.fits', dropstokes=True)
        if os.path.isfile(outfile+'.fits'):
            print2('Output to "%s"'%(outfile+'.fits'))
        else:
            if exit_on_error:
                raise Exception('Error! Failed to run CASA exportfits and output %s'%(outfile+'.fits'))
    # 




