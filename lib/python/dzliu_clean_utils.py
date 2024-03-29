# -*- coding: utf-8 -*-
# 
"""Utilities for running CASA tclean.

Notes
-----
Functions in this code must be run in CASA.  

Functions
---------
This code contains following functions (not a complete list):

- get_datacolumn
- get_antenn_diameter
- get_ref_frequency
- get_spw_frequency_range
- get_chan_width
- get_optimized_imsize
- get_field_phasecenters
- get_mosaic_imsize_and_phasecenter
- get_synbeam_and_imcell
- get_spw_for_spectral_line
- get_mstransform_params_for_spectral_line
- set_common_field_name_and_reference_dir
- cleanup_tclean_products
- apply_pbcor_to_tclean_image
- export_tclean_products_as_fits_files
- imsmooth_tclean_image
- load_params_from_dot_last_file

Last updates
------------
- 2020-12-23 copied functions from "dzliu_clean.py"
- 2021-01-04 updated tclean functions
- 2021-06-10 fixed CASA 6 import issue
- 2023-11-09 added load_params_from_dot_last_file

Example
-------
Example commands to run this code::

    import os, sys, glob
    sys.path.append(os.path.expanduser('~/Cloud/Github/Crab.Toolkit.CASA/lib/python'))
    from dzliu_combine_uvfits import dzliu_combine_uvfits
    dzliu_combine_uvfits(glob.glob('*.uvfits'), 'combined.uvfits')

"""
# 
from __future__ import print_function
import os, sys, re, json, copy, shutil
import numpy as np
from collections import OrderedDict
from collections import namedtuple
VersionInfo = namedtuple('VersionInfo', ['major', 'minor', 'micro', 'patch']) # 'releaselevel', serial'
str2versioninfo = lambda x: VersionInfo(*(np.pad(np.array(str(x).replace('-','.').split('.')).astype(int), (0,4), mode='constant', constant_values=0)[0:4]))
try:
    from astropy.io import fits
except:
    import pyfits as fits
try:
    #from taskinit import casalog, tb #, ms, iatool
    #from taskinit import casac
    #tb = casac.table
    #from __casac__.table import table as tb
    #from recipes import makepb, pixelmask2cleanmask
    try:
        import casadef
        casa_version_str = casadef.casa_version.replace('-','.')
        casa_version_info = str2versioninfo(casa_version_str)
    except:
        from casatools import utils as casatools_utils
        casa_version_str = casatools_utils.utils().version_string().replace('-','.')
        casa_version_info = str2versioninfo(casa_version_str)
    def _version_tuple(version_str):
        return tuple(map(int, (version_str.split("."))))
    def _version_less_than(version_str, compared_version_str):
        return _version_tuple(version_str) < _version_tuple(compared_version_str)
    def _version_greater_equal(version_str, compared_version_str):
        return _version_tuple(version_str) >= _version_tuple(compared_version_str)
    if _version_less_than(casa_version_str, '6.0.0'):
        from taskinit import casalog, tb #, ms, iatool
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
        # see CASA 6 usage here: https://casa.nrao.edu/casadocs/casa-5.6.0/introduction/casa6-installation-and-usage
        # see CASA 6 APIs: https://casadocs.readthedocs.io/en/stable/api/tt/casatools.table.html
        #from casatools import imager as imtool
        from casatools import table as casatools_table
        tb = casatools_table()
        from casatools import image as casatools_image
        iatool = casatools_image()
        from casatools import imager as casatools_imager
        imtool = casatools_imager()
        #from casatools import logsink as casatools_logsink
        #casalog = casatools_logsink() # this will create a new casalog. use the casalog in casatasks.
        from casatasks import casalog
        #from casatasks import listobs
        from casatasks import tclean, mstransform, exportfits, concat, split, imstat, impbcor, imsmooth
        #from casatasks import (concat, exportfits, feather, flagdata, imhead, immath, impbcor, importfits, 
        #    imrebin, imregrid, imsmooth, imstat, imsubimage, imtrans, imval, makemask, mstransform, 
        #    split, statwt, tclean, uvcontsub, visstat)
        casalog.post('CASA version >= 6', 'INFO')
        #from casatasks import sdbaseline
        #from casatools import ia
except:
    print('Error! Could not import taskinit and other CASA modules!')
    pass



# 
# def _print2
# 
def _print2(message):
    print(message)
    casalog.post(message, 'INFO')



# 
# def _print_params
# 
def _print_params(dict_params, prefix_str):
    print_str = prefix_str+'('+', '.join("{!s}={!r}".format(k, dict_params[k]) for k in dict_params.keys())+')'
    _print2(print_str)



# 
# json load byteified -- https://stackoverflow.com/questions/956867/how-to-get-string-objects-instead-of-unicode-from-json
# 
def _iteritems(dict_object):
    for key in dict_object:
        yield key, dict_object[key]

def _byteify(data, ignore_dicts = False):
    # if this is a unicode string, return its string representation
    if sys.version_info.major == 2 and isinstance(data, unicode):
        return data.encode('utf-8')
    # if this is a list of values, return list of byteified values
    if isinstance(data, list):
        return [ _byteify(item, ignore_dicts=True) for item in data ]
    # if this is a dictionary, return dictionary of byteified keys and values
    # but only if we haven't already byteified it
    if isinstance(data, dict) and not ignore_dicts:
        return {
            _byteify(key, ignore_dicts=True): _byteify(value, ignore_dicts=True)
            for key, value in _iteritems(data)
        }
    # if it's anything else, return it in its original form
    return data

def json_load_byteified(file_handle):
    return _byteify(
        json.load(file_handle, object_hook=_byteify),
        ignore_dicts=True
    )

def json_loads_byteified(json_text):
    return _byteify(
        json.loads(json_text, object_hook=_byteify),
        ignore_dicts=True
    )



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
    _print2('ant_diams = %s'%(ant_diams))
    _print2('minantdiam = %s [m]'%(minantdiam))
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
    casalog.origin('get_ref_frequency_Hz')
    #
    tb.open(vis+os.sep+'SPECTRAL_WINDOW')
    #spw_chan_freq = tb.getcol('CHAN_FREQ') # a list of list, Hz
    spw_ref_frequency = tb.getcol('REF_FREQUENCY') # Hz
    tb.close()
    #
    if spw_list is not None:
        if np.isscalar(spw_list):
            spw_list = [spw_list]
        ref_frequency = np.mean([spw_ref_frequency[int(t)] for t in spw_list])
    else:
        ref_frequency = spw_ref_frequency[0] # if no spw is specified, use the first spw REF_FREQUNCY
    # 
    return ref_frequency



#
# def get_spw_frequency_range
#
def get_spw_frequency_range(vis, spw_list=None):
    return get_spw_frequency_range_GHz(vis, spw_list=spw_list)
# 
def get_spw_frequency_range_GHz(vis, spw_list=None):
    #
    # Requires CASA module/function tb.
    #
    casalog.origin('get_spw_frequency_range_GHz')
    #
    tb.open(vis+os.sep+'SPECTRAL_WINDOW')
    spw_id = np.arange(0,tb.nrows()) # 
    spw_chan_freq = np.array([tb.getcell('CHAN_FREQ', i) for i in spw_id]) # a list of list, Hz
    spw_chan_freq_range = np.array([[t[0], t[-1]] for t in spw_chan_freq])
    tb.close()
    #
    if spw_list is not None:
        if np.isscalar(spw_list):
            spw_list = [spw_list]
        spw_freq_range_array_GHz = np.array([spw_chan_freq_range[int(t)] for t in spw_list]) / 1e9
        return spw_freq_range_array_GHz
    else:
        spw_freq_range_array_GHz = spw_chan_freq_range / 1e9 # return a list
        return spw_freq_range_array_GHz
    # 
    #return spw_freq_range_array_GHz



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
    casalog.origin('get_chan_width_MHz')
    #
    tb.open(vis+os.sep+'SPECTRAL_WINDOW')
    spw_id = np.arange(0,tb.nrows()) # 
    spw_chan_width = np.array([tb.getcell('CHAN_WIDTH', i)[0] for i in spw_id]) # a list of list, Hz
    tb.close()
    #
    if spw_list is not None:
        if np.isscalar(spw_list):
            spw_list = [spw_list]
        chan_width_array_MHz = np.array([spw_chan_width[int(t)] for t in spw_list]) / 1e6
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
# def get field names
#
def get_field_names(vis):
    """
    Read the FIELD table of the input measurement set and output all field names.
    """
    tb.open(vis+os.sep+'FIELD')
    field_names = tb.getcol('NAME')
    tb.close()
    #
    return field_names



#
# def get field phasecenters
#
def get_field_phasecenters(vis, galaxy_name = '', column_name = 'DELAY_DIR', return_all_names = False):
    """
    Get Measurement Set phase centers ('DELAY_DIR'). 
    
    Return 1 string and 2 lists: matched_field_name, matched_field_indices, and matched_field_phasecenters. 
    
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
    matched_field_names = []
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
            matched_field_names.append(field_name)
            matched_field_indices.append(i)
            matched_field_phasecenters.append([field_RA_deg, field_Dec_deg])
    #
    if '' == matched_field_name:
        raise ValueError('Error! Target source %s was not found in the "FIELD" table of the input vis "%s"!'%(galaxy_name, vis))
    #
    matched_field_indices = np.array(matched_field_indices)
    matched_field_phasecenters = np.array(matched_field_phasecenters).T # two columns, nrows
    if return_all_names:
        return matched_field_names, matched_field_indices, matched_field_phasecenters
    return matched_field_name, matched_field_indices, matched_field_phasecenters




#
# def get mosaic width and height in degree
#
def get_mosaic_imsize_and_phasecenter(vis, cell, galaxy_name='', ref_freq_Hz=None, padding_by_primary_beam=0.5, no_optimized_imsize=False, 
    divide_into_ncol_and_nrow = None, output_ds9_region_file = '', verbose=True):
    """
    cell is the same as pixel_size, and overrides pixel_size. 
    pixel_size can be a string or a float number. If it is a float number, a unit of arcsec is assumed.
    """
    # 
    casalog.origin('get_mosaic_imsize_and_phasecenter')
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
    matched_field_center_RA_deg = (matched_field_max_RA_deg+matched_field_min_RA_deg)/2.0 # BUGGY: np.mean(matched_field_phasecenters[0, :])
    matched_field_center_Dec_deg = (matched_field_max_Dec_deg+matched_field_min_Dec_deg)/2.0 # BUGGY: np.mean(matched_field_phasecenters[1, :])
    phasecenter = 'J2000 %.10fdeg %.10fdeg'%(matched_field_center_RA_deg, matched_field_center_Dec_deg)
    #print('matched_field_phasecenters.shape', matched_field_phasecenters.shape)
    #print('matched_field_phasecenters:', matched_field_phasecenters)
    #raise NotImplementedError()
    # 
    # calc primary beam
    if ref_freq_Hz is None:
        ref_freq_Hz = get_ref_frequency(vis)
    pribeam = 1.13  * (2.99792458e8 / ref_freq_Hz / minantdiam / np.pi * 180.0 ) # in units of degrees, see -- https://help.almascience.org/index.php?/Knowledgebase/Article/View/90
    if verbose:
        _print2('minantdiam = %s [meter]'%(minantdiam))
        _print2('pribeam = %s [arcsec]'%(pribeam * 3600.0))
        _print2('matched_field_phasecenters = %s'%(re.sub(r'[ \t\n\r]+', r' ', str(matched_field_phasecenters))))
        _print2('matched_field_min_RA_deg = %s'%(matched_field_min_RA_deg))
        _print2('matched_field_max_RA_deg = %s'%(matched_field_max_RA_deg))
        _print2('matched_field_min_Dec_deg = %s'%(matched_field_min_Dec_deg))
        _print2('matched_field_max_Dec_deg = %s'%(matched_field_max_Dec_deg))
    # 
    # calc mosaic width and height, half primary beam padding at both sides are considered.
    imsize_RA_deg = (matched_field_max_RA_deg - matched_field_min_RA_deg) * np.cos(np.deg2rad((matched_field_max_Dec_deg+matched_field_min_Dec_deg)/2.0))
    imsize_Dec_deg = (matched_field_max_Dec_deg - matched_field_min_Dec_deg)
    imsize_RA_deg_padded = imsize_RA_deg + 2.0 * padding_by_primary_beam * pribeam # padding this size at each side
    imsize_Dec_deg_padded = imsize_Dec_deg + 2.0 * padding_by_primary_beam * pribeam # padding this size at each side
    if verbose:
        _print2('imsize_RA = %s [arcsec]'%(imsize_RA_deg * 3600.0))
        _print2('imsize_Dec = %s [arcsec]'%(imsize_Dec_deg * 3600.0))
        _print2('imsize_RA = %s [arcsec] (padded)'%(imsize_RA_deg_padded * 3600.0))
        _print2('imsize_Dec = %s [arcsec] (padded)'%(imsize_Dec_deg_padded * 3600.0))
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
        _print2('Error! The input pixel_size could not be understood. It should be a string with a unit, e.g. \'1.0arcsec\'.')
        raise Exception('Error! The input pixel_size could not be understood.')
    if verbose:
        _print2('imcell = %s [arcsec]'%(imcell_arcsec))
    # 
    # 
    imsize_RA = imsize_RA_deg_padded / (imcell_arcsec / 3600.0) # pixels
    imsize_Dec = imsize_Dec_deg_padded / (imcell_arcsec / 3600.0) # pixels
    if no_optimized_imsize:
        imsize = [int(np.ceil(imsize_RA)), int(np.ceil(imsize_Dec))]
    else:
        imsize = [get_optimized_imsize(imsize_RA), get_optimized_imsize(imsize_Dec)]
    if verbose:
        _print2('imsize_RA = %s [pixel]'%(imsize_RA))
        _print2('imsize_Dec = %s [pixel]'%(imsize_Dec))
    # 
    # also write to ds9 region file
    if output_ds9_region_file != '':
        if os.path.isfile(output_ds9_region_file):
            shutil.move(output_ds9_region_file, output_ds9_region_file+'.backup')
        if output_ds9_region_file.find(os.sep) >= 0:
            if not os.path.isdir(os.path.dirname(output_ds9_region_file)):
                os.makedirs(os.path.dirname(output_ds9_region_file))
        with open(output_ds9_region_file, 'w') as fp:
            fp.write('fk5\n')
            fp.write('box(%s,%s,%s",%s",0.0) # text={%s}\n'%(\
                matched_field_center_RA_deg, matched_field_center_Dec_deg, imsize[0]*imcell_arcsec, imsize[1]*imcell_arcsec, vis))
            for i in range(len(matched_field_indices)):
                fp.write('circle(%s,%s,%s") # text={%s}\n'%(\
                    matched_field_phasecenters[0, i], matched_field_phasecenters[1, i], pribeam*3600./2., matched_field_indices[i]))
        if verbose:
            _print2('Output to "%s"'%(output_ds9_region_file))
    # 
    # if the user wants to divide the field into rows and cols, then do that
    if divide_into_ncol_and_nrow is not None:
        if np.isscalar(divide_into_ncol_and_nrow):
            divide_into_ncol_and_nrow = [divide_into_ncol_and_nrow]
        if len(divide_into_ncol_and_nrow) == 1:
            divide_into_ncol_and_nrow = [divide_into_ncol_and_nrow[0], divide_into_ncol_and_nrow[0]]
        ncol = divide_into_ncol_and_nrow[0]
        nrow = divide_into_ncol_and_nrow[1]
        _print2('Dividing into ncolxnrow %dx%d'%(ncol, nrow))
        divided_imsize_list = []
        divided_imsize_unpadded_list = []
        divided_center_RA_list = []
        divided_center_Dec_list = []
        divided_phasecenter_list = []
        for irow in range(nrow):
            for icol in range(ncol):
                divided_center_RA_deg = (matched_field_max_RA_deg - (matched_field_max_RA_deg - matched_field_min_RA_deg) / (2*ncol) * (2*icol + 1))
                divided_center_Dec_deg = (matched_field_min_Dec_deg + (matched_field_max_Dec_deg - matched_field_min_Dec_deg) / (2*nrow) * (2*irow + 1))
                divided_imsize_RA_deg = imsize_RA_deg / ncol
                divided_imsize_Dec_deg = imsize_Dec_deg / nrow
                divided_imsize_RA_deg_padded = divided_imsize_RA_deg + (imsize_RA_deg_padded - imsize_RA_deg)
                divided_imsize_Dec_deg_padded = divided_imsize_Dec_deg + (imsize_Dec_deg_padded - imsize_Dec_deg)
                divided_imsize_RA = divided_imsize_RA_deg_padded / (imcell_arcsec / 3600.0) # pixels
                divided_imsize_Dec = divided_imsize_Dec_deg_padded / (imcell_arcsec / 3600.0) # pixels
                divided_imsize_RA_unpadded = divided_imsize_RA_deg / (imcell_arcsec / 3600.0) # pixels, unpadded core size
                divided_imsize_Dec_unpadded = divided_imsize_Dec_deg / (imcell_arcsec / 3600.0) # pixels, unpadded core size
                divided_imsize = [int(np.ceil(divided_imsize_RA)), int(np.ceil(divided_imsize_Dec))]
                divided_imsize_unpadded = [int(np.ceil(divided_imsize_RA_unpadded)), int(np.ceil(divided_imsize_Dec_unpadded))]
                divided_phasecenter = 'J2000 %.10fdeg %.10fdeg'%(divided_center_RA_deg, divided_center_Dec_deg)
                divided_imsize_list.append(divided_imsize)
                divided_imsize_unpadded_list.append(divided_imsize_unpadded)
                divided_center_RA_list.append(divided_center_RA_deg)
                divided_center_Dec_list.append(divided_center_Dec_deg)
                divided_phasecenter_list.append(divided_phasecenter)
                if verbose:
                    _print2('icol %d irow %d phasecenter %s imsize %sx%s'%(icol, irow, divided_phasecenter, divided_imsize[0], divided_imsize[1]))
        if output_ds9_region_file != '':
            with open(output_ds9_region_file, 'a') as fp:
                for i in range(len(divided_imsize_list)):
                    icol = i%ncol
                    irow = int(i/ncol)
                    fp.write('box(%s,%s,%s",%s",0.0) # text={divided mosaic %d %d}\n'%(\
                        divided_center_RA_list[i], divided_center_Dec_list[i], 
                        divided_imsize_list[i][0]*imcell_arcsec, divided_imsize_list[i][1]*imcell_arcsec, icol, irow))
                    fp.write('box(%s,%s,%s",%s",0.0) # text={divided mosaic %d %d} dash=1\n'%(\
                        divided_center_RA_list[i], divided_center_Dec_list[i], 
                        divided_imsize_unpadded_list[i][0]*imcell_arcsec, divided_imsize_unpadded_list[i][1]*imcell_arcsec, icol, irow))
            if verbose:
                _print2('Output divided mosaic regions to "%s"'%(output_ds9_region_file))
        return divided_imsize_list, divided_phasecenter_list
    # 
    return imsize, phasecenter



# 
# get field IDs in mosaic
# 
def get_field_IDs_in_mosaic(vis, cell=None, imsize=None, phasecenter=None, ref_freq_Hz=None, galaxy_name='', padding_by_primary_beam=0.75, verbose=True):
    """Get field IDs in mosaic.
    
    The `padding_by_primary_beam` arg should correspond to the pblimit. pblimit 0.2 means padding_by_primary_beam 1.5/2.=0.75. 
    Here we choose a padding_by_primary_beam of 0.8, which means we require the field center must be within the input field of view and
    has a distance to the edges of at least 0.8 primary beam. This makes sure all field pointings are within the input field of view, 
    even for a small pblimit (0.2). 
    """
    # 
    casalog.origin('get_field_IDs_in_mosaic')
    # 
    if verbose:
        _print2('Getting field IDs inside phasecenter %s imsize %s cell %s in vis %r'%(phasecenter, imsize, cell, vis))
    # 
    if cell is None or imsize is None or phasecenter is None:
        _print2('Error! cell is None or imsize is None or phasecenter is None!')
        raise Exception('Error! cell is None or imsize is None or phasecenter is None!')
    # 
    input_RA, input_Dec = phasecenter.replace('J2000','').split()
    if input_RA.endswith('deg'):
        input_RA = float(input_RA.replace('deg',''))
    else:
        try:
            input_RA = float(input_RA)
        except:
            _print2('Error! Please input a float number in units of deg for the RA in the input phasecenter!')
            raise Exception('Error! Please input a float number in units of deg for the RA in the input phasecenter!')
    if input_Dec.endswith('deg'):
        input_Dec = float(input_Dec.replace('deg',''))
    else:
        try:
            input_Dec = float(input_Dec)
        except:
            _print2('Error! Please input a float number in units of deg for the Dec in the input phasecenter!')
            raise Exception('Error! Please input a float number in units of deg for the Dec in the input phasecenter!')
    # 
    # check imsize and make it 2-element list
    if np.isscalar(imsize):
        imsize = [imsize]
    if len(imsize) == 1:
        imsize = [imsize[0], imsize[0]]
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
        _print2('Error! The input pixel_size could not be understood. It should be a string with a unit, e.g. \'1.0arcsec\'.')
        raise Exception('Error! The input pixel_size could not be understood.')
    # 
    # get antdiam
    minantdiam = get_antenn_diameter(vis)
    # 
    # calc primary beam
    if ref_freq_Hz is None:
        ref_freq_Hz = get_ref_frequency(vis)
    pribeam = 1.13  * (2.99792458e8 / ref_freq_Hz / minantdiam / np.pi * 180.0 ) # in units of degrees, see -- https://help.almascience.org/index.php?/Knowledgebase/Article/View/90
    # 
    center_RA = input_RA
    center_Dec = input_Dec
    imsize_RA_deg = imsize[0]*imcell_arcsec/3600.0
    imsize_Dec_deg = imsize[1]*imcell_arcsec/3600.0
    # 
    output_field_IDs = []
    matched_field_name, matched_field_indices, matched_field_phasecenters = get_field_phasecenters(vis, galaxy_name=galaxy_name)
    if verbose:
        if galaxy_name != '':
            _print2('Found %d fields with galaxy_name %r in the input ms %r'%(len(matched_field_indices), galaxy_name, vis))
        else:
            _print2('Found %d fields in the input ms %r'%(len(matched_field_indices), vis))
    for i in range(len(matched_field_indices)):
        field_RA = matched_field_phasecenters[0, i]
        field_Dec = matched_field_phasecenters[1, i]
        dRA = np.abs(center_RA-field_RA)*np.cos(np.deg2rad(field_Dec))
        dDec = np.abs(center_Dec-field_Dec)
        # 
        if verbose:
            _print2('Checking field ID %d RA Dec %s %s with half primary beam size %s versus the input mosaic RA Dec %s %s (offset %s %s) FoV %s %s phasecenter %s imsize %s cell %s'%(\
                matched_field_indices[i], field_RA, field_Dec, pribeam/2., center_RA, center_Dec, dRA, dDec, 
                imsize_RA_deg, imsize_Dec_deg, phasecenter, imsize, cell))
        # 
        # we must make sure the whole field 1.5*pribeam (pblimit=0.2) is inside the imsize, otherwise tclean will be wrapped at the edge when there are data outside the FoV.
        if dRA < imsize_RA_deg/2.-pribeam*padding_by_primary_beam and \
           dDec < imsize_Dec_deg/2.-pribeam*padding_by_primary_beam:
            # 
            if verbose:
                _print2('Found field ID %d RA Dec %s %s with half primary beam size %s inside the input mosaic RA Dec %s %s (offset %s %s) FoV %s %s phasecenter %s imsize %s cell %s'%(\
                    matched_field_indices[i], field_RA, field_Dec, pribeam/2., center_RA, center_Dec, dRA, dDec, 
                    imsize_RA_deg, imsize_Dec_deg, phasecenter, imsize, cell))
            # 
            output_field_IDs.append(matched_field_indices[i])
    # 
    return output_field_IDs



# 
# get uvdist
# 
def get_synbeam_and_imcell(vis, ref_freq_Hz = None, oversampling = 5.0):
    """Get beam size and imcell.
    
    Returns
    -------
    synbeam : str
        Beam size in string format in units of arcsec.
    imcell : str
        Cell size in string format in units of arcsec. 
    
    """
    # 
    tb.open(vis)
    uvw = tb.getcol('UVW') # shape (3, nrows), each row is a (u,v,w) array
    ##_print2('tb.query(\'FIELD_ID in [%s] AND DATA_DESC_ID in [%s] AND STATE_ID in [%s]\', \'UVW\')'%(','.join(matched_field_indices.astype(str)), ','.join(valid_data_desc_indicies.astype(str)), ','.join(valid_state_indices.astype(str))))
    ##result = tb.query('FIELD_ID in [%s] AND DATA_DESC_ID in [%s] AND STATE_ID in [%s]'%(','.join(matched_field_indices.astype(str)), ','.join(valid_data_desc_indicies.astype(str)), ','.join(valid_state_indices.astype(str))), 'UVW')
    #_print2('tb.query(\'FIELD_ID in [%s] AND STATE_ID in [%s]\', \'UVW\')'%(','.join(matched_field_indices.astype(str)), ','.join(valid_state_indices.astype(str))))
    #result = tb.query('FIELD_ID in [%s] AND STATE_ID in [%s]'%(','.join(matched_field_indices.astype(str)), ','.join(valid_state_indices.astype(str))), 'UVW')
    #uvw = result.getcol('UVW')
    tb.close()
    # 
    if ref_freq_Hz is None:
        ref_freq_Hz = get_ref_frequency(vis)
    # 
    uvdist = np.sqrt(np.sum(np.square(uvw[0:2, :]), axis=0))
    maxuvdist = np.max(uvdist)
    _print2('maxuvdist = %s [m]'%(maxuvdist))
    L80uvdist = np.percentile(uvdist, 80) # np.max(uvdist) # now I am using 90-th percentile of baselies, same as used by 'analysisUtils.py' pickCellSize() getBaselineStats(..., percentile=...)
    _print2('L80uvdist = %s [m] (80-th percentile)'%(L80uvdist))
    # 
    synbeam = 2.99792458e8 / ref_freq_Hz / maxuvdist / np.pi * 180.0 * 3600.0 # arcsec
    synbeam = 0.574 * 2.99792458e8 / ref_freq_Hz / L80uvdist / np.pi * 180.0 * 3600.0 # arcsec # .574lambda/L80, see 'analysisUtils.py' estimateSynthesizedBeamFromASDM()
    synbeam_nprec = 3 # keep 3 valid digits
    synbeam_ndigits = (synbeam_nprec-1) - int(np.floor(np.log10(synbeam))) # keep to these digits (precision) after point, e.g., 0.1523 -> nprec 2 -> round(0.1523*100)/100 = 0.15
    synbeam = (np.round(synbeam * 10**(synbeam_ndigits))) / 10**(synbeam_ndigits)
    #oversampling = 5.0
    imcell_arcsec = synbeam / oversampling
    imcell = '%sarcsec'%(imcell_arcsec)
    _print2('synbeam = %s [arcsec]'%(synbeam))
    _print2('imcell_arcsec = %s'%(imcell_arcsec))
    return synbeam, imcell


#
# def get spw for spectral line
#
def get_spw_for_spectral_line(vis, redshift=None, rest_freq_GHz=None, line_width_kms=None, 
        exclude_continuum_spw=True, 
        reverse_selection=False, 
        return_dict=False, 
        verbose=True, 
    ):
    #
    casalog.origin('get_spw_for_spectral_line')
    # 
    if redshift is None or rest_freq_GHz is None:
        _print2('Error! Please input redshift and rest_freq_GHz!')
        raise Exception('Error! Please input redshift and rest_freq_GHz!')
    # 
    rest_freq_Hz = rest_freq_GHz * 1e9
    line_freq_GHz = rest_freq_GHz / (1.0+redshift)
    if line_width_kms is None:
        _print2('Setting line_width_kms to default value 1000.0 km/s.')
        line_width_kms = 1000.0 # km/s
    line_width_MHz = line_width_kms/2.99792458e5*line_freq_GHz*1000. # km/s
    line_freq_range_Hz = [line_freq_GHz*1e9 - line_width_MHz/2.0*1e6, line_freq_GHz*1e9 + line_width_MHz/2.0*1e6]
    if verbose:
        _print2('line_width_MHz = %s'%(line_width_MHz))
        _print2('line_freq_range_Hz = %s'%(line_freq_range_Hz))
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
        is_continuum_spw = ((nchan <= 4) or chstep>=31.25e6)
        if verbose:
            _print2('spw %s, ch0 %s, chlast %s, chstep %s, nchan %s, is_continuum_spw %s'%(i, ch0, chlast, chstep, nchan, is_continuum_spw))
        if exclude_continuum_spw and is_continuum_spw and (not reverse_selection):
            continue
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
            if (not reverse_selection):
                spw_selection_dict[str(i)] = '%d~%d'%(chleft, chright)
                if spw_selection_str != '':
                    spw_selection_str += ','
                spw_selection_str += '%d:%d~%d'%(i, chleft, chright)
            else:
                if chleft-1 < 0:
                    rev_chleft = []
                elif chleft-1 == 0:
                    rev_chleft = [str(0)]
                else:
                    rev_chleft = [str(0), '%d'%(chleft-1)]
                if chright+1 > nchan:
                    rev_chright = []
                elif chright+1 == nchan-1:
                    rev_chright = ['%d'%(nchan-1)]
                else:
                    rev_chright = ['%d'%(chright+1), '%d'%(nchan-1)]
                spw_selection_dict[str(i)] = ';'.join([t for t in ['~'.join(rev_chleft), '~'.join(rev_chright)] if t != ''])
                if spw_selection_str != '':
                    spw_selection_str += ','
                spw_selection_str += '%d:%s'%(i, spw_selection_dict[str(i)])
        else:
            if reverse_selection:
                spw_selection_dict[str(i)] = '0~%d'%(nchan-1)
                if spw_selection_str != '':
                    spw_selection_str += ','
                spw_selection_str += '%d'%(i)
    # 
    if verbose:
        _print2('spw_selection_str = %r'%(spw_selection_str))
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
        exclude_continuum_spw=True, 
        check_same_chan_width=True, 
        verbose=True, 
    ):
    """Get a mstransform parameter dict for the given spectral line.
    
    The returned dict can be directly used to call the CASA `mstransform` task, for example:: 
        
        mstransform_params = get_mstransform_params_for_spectral_line(vis, outputvis, field='my_target', redshift=3.0, rest_freq_GHz=230.538, line_width_kms=1000.0, chan_width_kms=20.0)
        mstransform(**mstransform_params)
    
    Returns
    -------
    mstransform_params : dict
        A dictionary that can be used as the input to the CASA `mstransform` task.
    
    """
    #
    casalog.origin('get_spw_for_spectral_line')
    # 
    if field is None or redshift is None or rest_freq_GHz is None or line_width_kms is None or chan_width_kms is None:
        _print2('Error! Please input field, redshift, rest_freq_GHz, line_width_kms and chan_width_kms!')
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
    spw_selection_str, spw_selection_dict = get_spw_for_spectral_line(\
        vis, redshift=redshift, rest_freq_GHz=rest_freq_GHz, line_width_kms=line_width_kms, 
        exclude_continuum_spw=exclude_continuum_spw, 
        return_dict=True, 
        )
    # 
    if spw_selection_str == '':
        _print2('Error! No spw found containing the input line!')
        return mstransform_params
    # 
    output_chan_width_MHz = (chan_width_kms/2.99792458e5*line_freq_MHz)
    # 
    # check channel width to be the same
    list_chan_width_MHz = get_chan_width_MHz(vis, spw_list=list(spw_selection_dict.keys()))
    if verbose:
        _print2('list_chan_width_MHz: %s, spw_selection_dict: %s'%(list_chan_width_MHz, spw_selection_dict))
    list_chan_width_MHz = np.abs(np.array(list_chan_width_MHz)).tolist()
    list_chan_width_MHz = dict(zip(list(spw_selection_dict.keys()), list_chan_width_MHz)) # convert to dict
    min_chan_width_MHz = None
    for i,ispw in enumerate(list(spw_selection_dict.keys())):
        if list_chan_width_MHz[ispw] > output_chan_width_MHz:
            if verbose:
                _print2('Discarding spw %s because its channel width %s MHz is broader than the output channel width %s MHz'%(ispw, list_chan_width_MHz[ispw], output_chan_width_MHz))
            del spw_selection_dict[ispw]
            del list_chan_width_MHz[ispw]
        else:
            if min_chan_width_MHz is None:
                min_chan_width_MHz = list_chan_width_MHz[ispw]
            elif min_chan_width_MHz > list_chan_width_MHz[ispw]:
                min_chan_width_MHz = list_chan_width_MHz[ispw]
    if min_chan_width_MHz is None:
        _print2('Error! No spw found with channel width matching to or better than the output channel width %s MHz!'%(output_chan_width_MHz))
        return mstransform_params
    # 
    if check_same_chan_width:
        for i,ispw in enumerate(list(spw_selection_dict.keys())):
            if not np.isclose(list_chan_width_MHz[ispw], min_chan_width_MHz, atol=0.0, rtol=0.01):
                if verbose:
                    _print2('Discarding spw %s because its channel width %s MHz is different and broader than the best channel width %s MHz'%(ispw, list_chan_width_MHz[ispw], min_chan_width_MHz))
                del spw_selection_dict[ispw]
                del list_chan_width_MHz[ispw]
    # 
    if force_integer_chan_width:
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
        _print2('mstransform_params: mstransform('+', '.join("{!s}={!r}".format(k, mstransform_params[k]) for k in mstransform_params.keys())+')')
    # 
    return mstransform_params



# 
# 
# 
def parse_RA_Dec_str(RA_Dec_str):
    """Convert RA Dec string to degrees. 
    """
    #<TODO># Not Tested!
    pattern_1 = re.compile(r'^([0-9+-]+)h([0-9]+)m([0-9.]+)s[ \t]+([0-9+-]+)d([0-9]+)m([0-9.]+)s')
    pattern_2 = re.compile(r'^([0-9+-]+):([0-9]+):([0-9.]+)[ \t]+([0-9+-]+):([0-9]+):([0-9.]+)')
    pattern_3 = re.compile(r'^([0-9+-]+):([0-9]+):([0-9.]+)[ \t]+([0-9+-]+)\.([0-9]+)\.([0-9.]+)')
    pattern_4 = re.compile(r'^([0-9.+-]+)deg[ \t]+([0-9.+-]+)deg')
    pattern_5 = re.compile(r'^([0-9.+-]+)[ \t]+([0-9.+-]+)')
    matched = False
    RAh, RAm, RAs = np.nan, np.nan, np.nan
    DEd, DEm, DEs = np.nan, np.nan, np.nan
    RAdeg = np.nan
    DEdeg = np.nan
    if not matched:
        match_1 = pattern_1.match(RA_Dec_str, re.IGNORECASE)
        if match_1:
            matched = True
            RAh, RAm, RAs, DEd, DEm, DEs = [float(t) for t in match_1.groups()]
    if not matched:
        match_2 = pattern_2.match(RA_Dec_str, re.IGNORECASE)
        if match_2:
            matched = True
            RAh, RAm, RAs, DEd, DEm, DEs = [float(t) for t in match_2.groups()]
    if not matched:
        match_3 = pattern_3.match(RA_Dec_str, re.IGNORECASE)
        if match_3:
            matched = True
            RAh, RAm, RAs, DEd, DEm, DEs = [float(t) for t in match_3.groups()]
    if not matched:
        match_4 = pattern_4.match(RA_Dec_str, re.IGNORECASE)
        if match_4:
            matched = True
            RAdeg, DEdeg = [float(t) for t in match_4.groups()]
    if not matched:
        match_5 = pattern_5.match(RA_Dec_str, re.IGNORECASE)
        if match_5:
            matched = True
            RAdeg, DEdeg = [float(t) for t in match_5.groups()]
    if matched:
        if np.isnan(RAdeg):
            RAdeg = RAh*15. + RAm/60. + RAs/3600.
            DEdeg = DEd + DEm/60. + DEs/3600.
    #<TODO># Not Tested!
    return RAdeg, DEdeg



# 
# 
# 
def set_common_field_name_and_reference_dir(vis, source_name, source_RA, source_Dec, outputvis, overwrite=False):
    """Set a common field name and reference direction for all fields in the input measurement set. 
    """
    if os.path.isdir(outputvis):
        if not os.path.isfile(outputvis+os.sep+'table.dat'):
            print('Warning! Found existing output data %r but it is incomplete. Deleting it.'%(outputvis))
            shutil.rmtree(outputvis)
    # 
    if os.path.isdir(outputvis):
        if os.path.isdir(outputvis+'.touch'):
            print('Warning! Found existing output data %r but it seems still being touched by some process. Skipping it.'%(outputvis))
            return
    # 
    if os.path.isdir(outputvis):
        if overwrite:
            shutil.rmtree(outputvis)
        else:
            print('Found existing output data %r and overwrite is False. Will not overwrite it.'%(outputvis))
            return
    # 
    os.makedirs(outputvis+'.touch')
    # 
    shutil.copytree(vis, outputvis)
    # 
    tb.open(outputvis+os.sep+'FIELD', nomodify=False)
    fplog = open(outputvis+'.modifying.field.table.txt', 'w')
    fplog.write('row, old name, new name, old reference dir a0, old reference dir d0, new reference dir a0, new reference dir d0\n')
    for i in range(tb.nrows()):
        old_field_name = tb.getcell('NAME', i)
        old_reference_dir = tb.getcell('REFERENCE_DIR', i)
        new_field_name = source_name
        new_reference_dir = [[np.deg2rad(source_RA)], [np.deg2rad(source_Dec)]]
        fplog.write('%d, %r, %r, %s, %s, %s, %s\n'%(i, old_field_name, new_field_name, \
            old_reference_dir[0][0], old_reference_dir[0][1], new_reference_dir[0][0], new_reference_dir[0][1]))
        tb.putcell('NAME', i, new_field_name)
        tb.putcell('REFERENCE_DIR', i, new_reference_dir)
    fplog.close()
    tb.flush()
    tb.close()
    # 
    #tb.open(outputvis+os.sep+'SOURCE', nomodify=False)
    #fplog = open(outputvis+'.modifying.source.table.txt', 'w')
    #fplog.write('row, old name, new name, old direction a0, old direction d0, new direction a0, new direction d0\n')
    #for i in [0]:
    #    old_name = tb.getcell('NAME', i)
    #    old_direction = tb.getcell('DIRECTION', i)
    #    new_name = new_name
    #    new_direction = [np.deg2rad(source_RA), np.deg2rad(source_Dec)]
    #    fplog.write('%d, %r, %r, %s, %s, %s, %s\n'%(i, old_name, new_name, \
    #        old_direction[0], old_direction[1], new_direction[0], new_direction[1]))
    #    tb.putcell('NAME', i, new_name)
    #    tb.putcell('DIRECTION', i, new_direction)
    #    tb.putcell('CALIBRATION_GROUP', i, 0)
    #    tb.putcell('CODE', i, '')
    #fplog.close()
    #tb.flush()
    #tb.close()
    # 
    os.rmdir(outputvis+'.touch')
    # 
    print('Output to "%s"'%(outputvis))
    # 
    return



# 
# 
# 
def cleanup_tclean_products(imagename, suffix_list=None, cleanup_mask=True, cleanup_fits=True, exit_on_error=True):
    """Delete all tclean products given a base image name.
    """
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
                _print2('Deleted "%s"'%(imagename+suffix))
            else:
                if exit_on_error:
                    raise Exception('Error! Failed to cleanup tclean product data directory %s'%(imagename+suffix))
        if cleanup_fits:
            if os.path.isfile(imagename+suffix+'.fits'):
                os.remove(imagename+suffix+'.fits')
                if not os.path.isfile(imagename+suffix+'.fits'):
                    _print2('Deleted "%s"'%(imagename+suffix+'.fits'))
                else:
                    if exit_on_error:
                        raise Exception('Error! Failed to cleanup tclean product data file %s'%(imagename+suffix+'.fits'))


def apply_pbcor_to_tclean_image(imagename, cutoff=0.1, dropstokes=True, velocity=True, overwrite=True, exit_on_error=True):
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
                    _print2('Found existing data "%s", overwriting it.'%(outfile2))
                    shutil.rmtree(outfile2)
            _print2('Running CASA task: impbcor(imagename=%r, pbimage=%r, outfile=%r, mode=%r, cutoff=%s)'%(infile2, pbimage2, outfile2, 'divide', cutoff))
            impbcor(imagename=infile2, pbimage=pbimage2, outfile=outfile2, mode='divide', cutoff=cutoff)
            if os.path.isdir(outfile2):
                _print2('Output to "%s"'%(outfile2))
            else:
                if exit_on_error:
                    raise Exception('Error! Failed to run CASA impbcor and output "%s"'%(outfile2))
            # 
            export_tclean_products_as_fits_files(imagename, suffix_list=['.image.pbcor'+suffix], dropstokes=dropstokes, velocity=velocity)


def export_tclean_products_as_fits_files(imagename, dropstokes=True, velocity=False, suffix_list=None, overwrite=True, exit_on_error=True):
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
                        _print2('Found existing tclean product fits file "%s", will not overwrite it.'%(outfile))
                        continue
                else:
                    _print2('Found existing tclean product fits file "%s", overwriting it.'%(outfile))
                    os.remove(outfile)
            _print2('Running CASA task: exportfits(%r, %r, dropstokes=%s, velocity=%s)'%(infile, outfile, dropstokes, velocity))
            exportfits(infile, outfile, dropstokes=dropstokes, velocity=velocity)
            if os.path.isfile(outfile):
                _print2('Output to "%s"'%(outfile))
            else:
                raise Exception('Error! Failed to run CASA exportfits and output %s'%(imagename+suffix))
    # 


def imsmooth_tclean_image(infile, major, minor=None, pa=None, kernel='gaussian', targetres=True, outfile=None, 
                          overwrite=True, exit_on_error=True, export_fits=True, dropstokes=True, velocity=False):
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
                _print2('Found existing tclean product fits file "%s", will not overwrite it.'%(outfile))
                run_imsmooth = False
        else:
            _print2('Found existing tclean product fits file "%s", overwriting it.'%(outfile))
            shutil.rmtree(outfile)
    # run imsmooth
    if run_imsmooth:
        _print2('Running CASA task: imsmooth(imagename = %r, kernel = %r, major = %r, minor = %r, pa = %r, targetres = %s, outfile = %r)'%(infile, kernel, major, minor, pa, targetres, outfile))
        imsmooth(imagename = infile, kernel = kernel, major = major, minor = minor, pa = pa, targetres = targetres, outfile = outfile)
    # check result
    if os.path.isdir(outfile):
        _print2('Output to "%s"'%(outfile))
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
                _print2('Found existing tclean product fits file "%s", will not overwrite it.'%(outfile))
                run_export_fits = False
        else:
            _print2('Found existing tclean product fits file "%s", overwriting it.'%(outfile))
            os.remove(outfile+'.fits')
    # export fits file
    if export_fits and run_export_fits:
        exportfits(outfile, outfile+'.fits', dropstokes=dropstokes, velocity=velocity)
        if os.path.isfile(outfile+'.fits'):
            _print2('Output to "%s"'%(outfile+'.fits'))
        else:
            if exit_on_error:
                raise Exception('Error! Failed to run CASA exportfits and output %s'%(outfile+'.fits'))
    # 


def load_params_from_dot_last_file(dot_last_file):
    dict_params = {}
    with open(dot_last_file, 'r') as fp:
        for line in fp:
            if line.startswith('#'): 
                continue
            line_str = line.strip()
            if line_str == '':
                continue
            line_match = re.match(r'^(.+) *# *(.+) *$', line_str)
            if line_match:
                line_str = line_match.group(1)
            line_match = re.match(r'^([^ ]+) *= *(.+)$', line_str)
            if line_match:
                dict_params[line_match.group(1)] = eval(line_match.group(2))
    if 'taskname' in dict_params:
        del dict_params['taskname']
    return dict_params




