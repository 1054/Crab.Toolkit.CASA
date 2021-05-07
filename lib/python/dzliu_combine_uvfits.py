#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
"""Utilities for combining uvfits.

Notes
-----
Functions in this code must be run in CASA.  

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
try:
    from astropy.io import fits
except:
    import pyfits as fits
try:
    from taskinit import casalog, tb #, ms, iatool
    #from taskinit import casac
    #tb = casac.table
    #from __casac__.table import table as tb
    #from recipes import makepb, pixelmask2cleanmask
    import casadef
    def _version_tuple(version_str):
        return tuple(map(int, (version_str.split("."))))
    def _version_less_than(version_str, compared_version_str):
        return _version_tuple(version_str) < _version_tuple(compared_version_str)
    def _version_greater_equal(version_str, compared_version_str):
        return _version_tuple(version_str) >= _version_tuple(compared_version_str)
    if _version_less_than(casadef.casa_version, '6.0.0'):
        from mstransform_cli import mstransform_cli_
        mstransform = mstransform_cli_()
        from concat_cli import concat_cli_
        concat = concat_cli_()
        from importuvfits_cli import importuvfits_cli_
        importuvfits = importuvfits_cli_()
        from exportuvfits_cli import exportuvfits_cli_
        exportuvfits = exportuvfits_cli_()
    else:
        # see CASA 6 updates here: https://alma-intweb.mtk.nao.ac.jp/~eaarc/UM2018/presentation/Nakazato.pdf
        from casatasks import mstransform, concat, importuvfits, exportuvfits
        #from casatasks import sdbaseline
        #from casatools import ia
except:
    print('Error! Could not import taskinit and other CASA modules!')
    pass



# 
# def _print2
# 
def _print2(message):
    """Print message on screen as well as in casalog.
    """
    print(message)
    casalog.post(message, 'INFO')



# 
# def _print_params
# 
def _print_params(prefix_str, dict_params):
    print_str = prefix_str+'('+', '.join("{!s}={!r}".format(k, dict_params[k]) for k in dict_params.keys())+')'
    _print2(print_str)



# 
# def _save_params
# 
class _numpy_array_encoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

def _save_params(dict_params, json_file):
    if json_file.find(os.sep) >= 0:
        output_dir = os.path.dirname(json_file)
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
    with open(json_file, 'w') as fp:
        json.dump(dict_params, fp, indent = 4, cls = _numpy_array_encoder)



# 
# def _get_spw_info_dict
# 
def _get_spw_info_dict(vis, target_frequency=None, output_json_file=None, only_line_spw=True):
    tb.open(vis+os.sep+'SPECTRAL_WINDOW')
    spw_info_dict = {}
    spw_count = tb.nrows()
    #print(tb.colnames()) # ['MEAS_FREQ_REF', 'CHAN_FREQ', 'REF_FREQUENCY', 'CHAN_WIDTH', 'EFFECTIVE_BW', 'RESOLUTION', 'FLAG_ROW', 'FREQ_GROUP', 'FREQ_GROUP_NAME', 'IF_CONV_CHAIN', 'NAME', 'NET_SIDEBAND', 'NUM_CHAN', 'TOTAL_BANDWIDTH']
    #spw_info_dict['count'] = spw_count
    for i in range(spw_count):
        spw_name = tb.getcell('NAME', i)
        num_chan = tb.getcell('NUM_CHAN', i)
        chan_freq = tb.getcell('CHAN_FREQ', i)
        chan_width = tb.getcell('CHAN_WIDTH', i)
        ref_freq = tb.getcell('REF_FREQUENCY', i)
        if not np.isscalar(spw_name): 
            spw_name = spw_name[0]
        if not np.isscalar(num_chan): 
            num_chan = num_chan[0]
        if not np.isscalar(chan_width): 
            chan_width = np.mean(chan_width)
        if not np.isscalar(ref_freq): 
            ref_freq = np.mean(ref_freq)
        min_freq = np.min(chan_freq)
        max_freq = np.max(chan_freq)
        # if only_line_spw is True, make sure this spw is a spectral line spw with more than one channel and in ALMA FDM mode
        if only_line_spw:
            if num_chan == 1 or chan_width >= 31.25e6 or spw_name.find('WVR')>=0 or spw_name.find('AVG')>=0:
                continue
        # if target_frequency is provided, make sure this spw covers this targetting frequency
        if target_frequency is not None:
            if (min_freq - target_frequency) * (max_freq - target_frequency) > 0:
                continue
        spw_info_dict[i] = {
            'NAME': spw_name, 
            'NUM_CHAN': num_chan, 
            'CHAN_FREQ': chan_freq, 
            'CHAN_WIDTH': chan_width, 
            'REF_FREQUENCY': ref_freq,
            'MIN_FREQ': min_freq, 
            'MAX_FREQ': max_freq, 
        }
    tb.close()
    if len(spw_info_dict) > 0 and output_json_file is not None:
        _save_params(spw_info_dict, output_json_file)
    return spw_info_dict



# 
# def _get_datacolumn
# 
def _get_datacolumn(vis):
    casalog.origin('get_datacolumn')
    tb.open(vis)
    if 'CORRECTED_DATA' in tb.colnames():
        datacolumn = 'CORRECTED'
    else:
        datacolumn = 'DATA'
    tb.close()
    return datacolumn



# 
# def _get_field_info_dict
# 
def _get_field_info_dict(vis, target_ra, target_dec, separation_limit=2.0, output_json_file=None):
    tb.open(vis+os.sep+'FIELD')
    field_info_dict = {}
    field_count = tb.nrows()
    #print(tb.colnames()) # ['DELAY_DIR', 'PHASE_DIR', 'REFERENCE_DIR', 'CODE', 'FLAG_ROW', 'NAME', 'NUM_POLY', 'SOURCE_ID', 'TIME', 'EPHEMERIS_ID', 'PhaseDir_Ref', 'DelayDir_Ref', 'RefDir_Ref']
    #field_info_dict['count'] = field_count
    for i in range(field_count):
        field_name = tb.getcell('NAME', i)
        delay_dir = tb.getcell('DELAY_DIR', i).ravel()
        phase_dir = tb.getcell('PHASE_DIR', i).ravel()
        reference_dir = tb.getcell('REFERENCE_DIR', i).ravel()
        # if target_ra and target_dec are provided, make sure this field covers this targetting coordiante
        if target_ra is not None and target_dec is not None:
            target_ra_radian = np.deg2rad(target_ra)
            target_dec_radian = np.deg2rad(target_dec)
            separation_limit_radian = np.deg2rad(separation_limit/3600.0)
            if ((delay_dir[0]-target_ra_radian)**2 + (delay_dir[1]-target_dec_radian)**2) > separation_limit_radian**2:
                continue
        field_info_dict[i] = {
            'NAME': field_name, 
            'DELAY_DIR': delay_dir, 
            'PHASE_DIR': phase_dir, 
            'REFERENCE_DIR': reference_dir, 
        }
    tb.close()
    if len(field_info_dict) > 0 and output_json_file is not None:
        _save_params(field_info_dict, output_json_file)
    return field_info_dict






# 
# def dzliu_combine_uvfits
# 
def dzliu_combine_uvfits(
        list_of_uvfits, 
        output_uvfits, 
        target_frequency, 
        target_ra, 
        target_dec, 
        separation_limit = 4.0, 
        overwrite = False, 
    ):
    """Combine a list of input uvfits files. 
    
    Args
    ----
    list_of_uvfits : list
        A list of uvfits as the input.
    output_uvfits : str
        The output uvfits file path.
    target_frequency : float
        The targetting frequency which will be covered by the spectral windows of the data to combine.
    target_ra : float
        The targetting R.A. coordinate in degrees.
    target_dec : float
        The targetting Decl. coordinate in degrees.
    separation_limit : float
        The separation limit from each dataset phase center to the target R.A. and Decl. coordinate.
    overwrite : bool
        Whether to overwrite the output file or not.
    
    Returns
    -------
    None
    
    Notes
    -----
    Note that currently we do not apply any phase shift and primary beam correction. 
    We simply assume that the primary beam is much larger than the `separation_limit`.
    
    """
    # 
    casalog.origin('dzliu_combine_uvfits')
    
    # set working dir
    working_dir = re.sub(r'\.uvfits$', r'', output_uvfits) + '_combining_uvfits_working_dir'
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)
    
    # importuvfits
    list_of_input_ms = []
    for i in range(len(list_of_uvfits)):
        fitsfile = list_of_uvfits[i]
        vis = working_dir+os.sep+'input_vis_%d.ms'%(i+1)
        if not os.path.isdir(vis) and fitsfile.endswith('.ms'):
            # if input is ms, directly copy ms
            _print2('Copying MS from "%s" to "%s"'%(fitsfile, vis))
            shutil.copytree(fitsfile, vis)
            if not os.path.isdir(vis):
                raise Exception('Error! Failed to copy the input MS to "%s"!'%(vis))
            else:
                _print2('Output to "%s"!'%(vis))
        elif not os.path.isdir(vis):
            importuvfits_params = OrderedDict()
            importuvfits_params['fitsfile'] = fitsfile
            importuvfits_params['vis'] = vis
            importuvfits_params['antnamescheme'] = 'new'
            _print_params('importuvfits', importuvfits_params)
            _save_params(importuvfits_params, vis+'.importuvfits.params.json')
            importuvfits(**importuvfits_params)
            if not os.path.isdir(vis):
                raise Exception('Error! Failed to run importuvfits and produce "%s"!'%(vis))
            else:
                _print2('Output to "%s"!'%(vis))
        else:
            _print2('Found existing "%s". Will not overwrite it!'%(vis))
        list_of_input_ms.append(vis)
    
    # find out common frequency range and coarest channel width
    list_of_spw_info_dict = []
    list_of_field_info_dict = []
    list_of_input_ms_dict = []
    for i in range(len(list_of_input_ms)):
        vis = list_of_input_ms[i]
        spw_info_dict = _get_spw_info_dict(vis, target_frequency, output_json_file=vis+'.spw.info.json', only_line_spw=True)
        field_info_dict = _get_field_info_dict(vis, target_ra, target_dec, separation_limit, output_json_file=vis+'.field.info.json')
        list_of_spw_info_dict.append(spw_info_dict)
        list_of_field_info_dict.append(field_info_dict)
        if len(spw_info_dict) > 0 and len(field_info_dict) > 0:
            input_ms_dict = {}
            input_ms_dict['spw_info_dict'] = spw_info_dict
            input_ms_dict['field_info_dict'] = field_info_dict
            input_ms_dict['target_frequency'] = target_frequency
            input_ms_dict['min_freq'] = np.min([spw_info_dict[k]['MIN_FREQ'] for k in spw_info_dict.keys()])
            input_ms_dict['max_freq'] = np.max([spw_info_dict[k]['MAX_FREQ'] for k in spw_info_dict.keys()])
            input_ms_dict['chan_width'] = np.max([spw_info_dict[k]['CHAN_WIDTH'] for k in spw_info_dict.keys()]) #<TODO># taking the coarest channel width
            input_ms_dict['spw'] = ','.join([str(k).strip() for k in spw_info_dict.keys()])
            input_ms_dict['field'] = ','.join([str(k).strip() for k in field_info_dict.keys()])
            input_ms_dict['vis'] = vis
            list_of_input_ms_dict.append(input_ms_dict)
    
    if len(list_of_input_ms_dict) == 0:
        _print2('Error! Could not find any input uvfits that contain the target frequency %s and RA Dec %s %s.'%(target_frequency, target_ra, target_dec))
        for i in range(len(list_of_input_ms)):
            _print_params('list_of_spw_info_dict[%d]: '%(i), list_of_spw_info_dict[i])
            _print_params('list_of_field_info_dict[%d]: '%(i), list_of_field_info_dict[i])
        raise Exception('Error! Could not find any input uvfits that contain the target frequency %s and RA Dec %s %s.'%(target_frequency, target_ra, target_dec))
    
    common_min_freq = np.max([input_ms_dict['min_freq'] for input_ms_dict in list_of_input_ms_dict])
    common_max_freq = np.min([input_ms_dict['max_freq'] for input_ms_dict in list_of_input_ms_dict])
    if common_min_freq >= common_max_freq:
        _print2('Error! Could not find a common intersected frequency range among the input uvfits that contain the target frequency %s.'%(target_frequency))
        for i in range(len(list_of_input_ms_dict)):
            _print_params('list_of_input_ms_dict[%d]: '%(i), list_of_input_ms_dict[i])
        raise Exception('Error! Could not find a common intersected frequency range among the input uvfits that contain the target frequency %s.'%(target_frequency))
    
    common_chan_width = np.max([np.abs(input_ms_dict['chan_width']) for input_ms_dict in list_of_input_ms_dict])
    
    # mstransform
    list_of_mstransformed_ms = []
    for i in range(len(list_of_input_ms_dict)):
        input_ms_dict = list_of_input_ms_dict[i]
        spw_info_dict = input_ms_dict['spw_info_dict']
        field_info_dict = input_ms_dict['field_info_dict']
        vis = input_ms_dict['vis']
        outputvis = working_dir+os.sep+'mstransformed_vis_%d.ms'%(i+1)
        if not os.path.isdir(outputvis):
            # check if the input vis has more than one spws and these spws do not have the same channel number
            #_print2('len(spw_info_dict): '+str(len(spw_info_dict)))
            #_print2(str([spw_info_dict[k]['NUM_CHAN'] for k in spw_info_dict.keys()]))
            #_print2(str(np.diff(np.array([spw_info_dict[k]['NUM_CHAN'] for k in spw_info_dict.keys()]))))
            if len(spw_info_dict) > 1 and np.max(np.diff(np.array([spw_info_dict[k]['NUM_CHAN'] for k in spw_info_dict.keys()]))) > 1:
                _print2('Warning! More than one spws are in the ms and they do not have the same channel number. ')
                argmin_chan_width = np.argmin(np.array([np.abs(spw_info_dict[k]['CHAN_WIDTH']) for k in spw_info_dict.keys()])).ravel()
                ispw_min_chan_width = np.array(list(spw_info_dict.keys()))[argmin_chan_width]
                if len(ispw_min_chan_width) > 1:
                    argmax_num_chan = np.argmax(np.array([np.abs(spw_info_dict[k]['NUM_CHAN']) for k in ispw_min_chan_width])).ravel()[0]
                    ispw_min_chan_width = ispw_min_chan_width[argmax_num_chan]
                else:
                    ispw_min_chan_width = ispw_min_chan_width[0]
                _print2('We will take the minimum chan_width spw %s, and exclude others.'%(ispw_min_chan_width))
                ispw_to_exclude = []
                for ispw in spw_info_dict.keys():
                    _print2('spw_info_dict[%d]: min_freq: %s, max_freq: %s, num_chan: %s, name: %r'%(ispw, 
                            spw_info_dict[ispw]['MIN_FREQ'], 
                            spw_info_dict[ispw]['MAX_FREQ'], 
                            spw_info_dict[ispw]['NUM_CHAN'], 
                            spw_info_dict[ispw]['NAME']))
                    if spw_info_dict[ispw]['CHAN_WIDTH'] > spw_info_dict[ispw_min_chan_width]['CHAN_WIDTH'] or \
                       spw_info_dict[ispw]['NUM_CHAN'] != spw_info_dict[ispw_min_chan_width]['NUM_CHAN']:
                        ispw_to_exclude.append(ispw)
                _print2('ispw_to_exclude: %s'%(str(ispw_to_exclude)))
                for ispw in ispw_to_exclude:
                    del spw_info_dict[ispw]
                input_ms_dict['spw'] = ','.join([str(k).strip() for k in spw_info_dict.keys()])
                input_ms_dict['spw_info_dict'] = spw_info_dict
                pass
            # 
            mstransform_params = OrderedDict()
            mstransform_params['vis'] = vis
            mstransform_params['outputvis'] = outputvis
            mstransform_params['datacolumn'] = _get_datacolumn(vis)
            mstransform_params['spw'] = input_ms_dict['spw']
            mstransform_params['field'] = input_ms_dict['field']
            mstransform_params['combinespws'] = True # this requires that each spw has the same channel number
            mstransform_params['nspw'] = 1 # output single spw ms
            mstransform_params['regridms'] = True
            mstransform_params['mode'] = 'frequency'
            mstransform_params['start'] = '%.0fHz'%(common_min_freq)
            mstransform_params['width'] = '%.0fHz'%(common_chan_width)
            mstransform_params['nchan'] = int(np.round((common_max_freq-common_min_freq)/common_chan_width))
            mstransform_params['mode'] = 'frequency'
            _print_params('mstransform', mstransform_params)
            _save_params(mstransform_params, outputvis+'.mstransform.params.json')
            mstransform(**mstransform_params)
            if not os.path.isdir(outputvis):
                raise Exception('Error! Failed to run mstransform and produce "%s"!'%(outputvis))
            else:
                _print2('Output to "%s"!'%(outputvis))
        else:
            _print2('Found existing "%s". Will not overwrite it!'%(outputvis))
        # 
        list_of_mstransformed_ms.append(outputvis)
    
    # concat
    concatvis = working_dir+os.sep+'concat.ms'
    if not os.path.isdir(concatvis):
        concat_params = OrderedDict()
        concat_params['vis'] = list_of_mstransformed_ms
        concat_params['concatvis'] = concatvis
        concat_params['dirtol'] = '%.3farcsec'%(max(separation_limit, 0.001))
        concat_params['freqtol'] = '%.0fHz'%(common_chan_width/10.0)
        _print_params('concat', concat_params)
        _save_params(concat_params, concatvis+'.concat.params.json')
        concat(**concat_params)
        if not os.path.isdir(concatvis):
            raise Exception('Error! Failed to run concat and produce "%s"!'%(concatvis))
        else:
            _print2('Output to "%s"!'%(concatvis))
    else:
        _print2('Found existing "%s". Will not overwrite it!'%(concatvis))

    # exportuvfits
    vis = concatvis
    fitsfile = output_uvfits
    if not os.path.isfile(fitsfile):
        exportuvfits_params = OrderedDict()
        exportuvfits_params['vis'] = vis
        exportuvfits_params['fitsfile'] = fitsfile
        _print_params('exportuvfits', exportuvfits_params)
        _save_params(exportuvfits_params, fitsfile+'.exportuvfits.params.json')
        exportuvfits(**exportuvfits_params)
        if not os.path.isfile(fitsfile):
            raise Exception('Error! Failed to run exportuvfits and produce "%s"!'%(fitsfile))
        else:
            _print2('Output to "%s"!'%(fitsfile))
    else:
        _print2('Found existing "%s". Will not overwrite it!'%(fitsfile))
    
    # Done
    _print2('Done!')









############
#   main   #
############

dzliu_main_func_name = 'dzliu_combine_uvfits' # make sure this is the right main function in this script file

if __name__ == '__main__':
    if 'casa' in globals():
        # We are already in CASA and are called via execfile
        dzliu_main_func = globals()[dzliu_main_func_name]
        dzliu_main_func(globals())
    else:
        print('Please run this in CASA via:')
        print('(Python2)')
        print('    execfile(\'%s\')'%(os.path.basename(__file__)))
        print('(Python3)')
        print('    from %s import %s'%(re.sub(r'\.py$', r'', os.path.basename(__file__)), dzliu_main_func_name) )
        print('    %s(globals())'%(dzliu_main_func_name) )
        raise Exception('Please see message above.')




