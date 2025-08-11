# -*- coding: utf-8 -*-
# 
"""Utilities for analyzing u-v visibility data within CASA.

Notes
-----
Functions in this code must be run in CASA.  

Functions
---------
This code contains following functions (not a complete list):

- TBD
- extract_field_selections
- extract_line_spw_selections
- flag_line_data
- estimate_tclean_params
- selfcal_continuum_data
- timebin_uvdata

Last updates
------------
- 2021-08-16 start
- 2025-06-26 updated

Example
-------
Example commands to run this code (in CASA environment)::

    import os, sys, glob
    sys.path.append(os.path.expanduser('~/Cloud/Github/Crab.Toolkit.CASA/lib/python'))
    from dzliu_uvdata_utils import plot_uvdist_vs_amp
    plot_uvdist_vs_amp(vis, plotfile)

"""

import os, sys, re, json, copy, glob, time, datetime, shutil
import numpy as np
import analysisUtils as aU
import importlib
importlib.reload(aU)
import inspect
import pprint



# 
# Funciton to check if we are inside CASA or not
# 
def check_casa(f_globals = None):
    if f_globals is None:
        f_globals = globals()
    if 'casac' in f_globals: # CASA 4
        return True
    elif '__casac__' in f_globals: # CASA 5
        return True
    elif 'casatasks' in f_globals: # CASA 6
        return True
    return False



# 
# Find caller scope
# 
caller_globals = {}
for iscope, scope in enumerate(inspect.stack()):
    if hasattr(scope, 'filename'):
        # casa 6.
        scope_filename = scope.filename
        scope_frame = scope.frame
    else:
        # casa 4.6.0
        scope_filename = scope[3]
        scope_frame = scope[0]
    #print('inspect.stack()', iscope, scope_filename)
    if check_casa(scope_frame.f_globals):
        caller_globals = scope_frame.f_globals
        break

if len(caller_globals) == 0:
    raise Exception('Error! Could not find CASA in inspect.stack()!')

# for iscope, scope in enumerate(inspect.stack()):
#     if hasattr(scope, 'filename'):
#         # casa 6.
#         scope_filename = scope.filename
#         scope_frame = scope.frame
#     else:
#         # casa 4.6.0
#         scope_filename = scope[3]
#         scope_frame = scope[0]
#     print('inspect.stack()', iscope, scope_filename)
#     if scope_filename.find('importlib')>=0 or scope_filename.find('<')>=0 or scope_filename in ['runcode', 'runsource', 'push', 'interact']:
#         continue
#     if iscope == 0:
#         continue
#     caller_globals = scope_frame.f_globals
#     print('caller_globals', caller_globals.keys())
#     break

#print('caller_globals', caller_globals.keys())
if check_casa(caller_globals):
    for key in caller_globals:
        if key not in globals(): 
            if key.find('casa')>=0 or key in [
                'tb', 'ia', 'me', 'qa', 
                'flagdata', 'split', 'mstransform', 'tclean', 'imstat', 'gaincal', 'plotms', 'applycal'
            ]:
                #print('caller_globals[{!r}]'.format(key))
                globals()[key] = caller_globals[key]



# 
# class for json dump function (cls=NpEncoder)
# -- for json dump and load 
# -- https://stackoverflow.com/questions/27050108/convert-numpy-type-to-python
# 
class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)



# 
# Function to convert ra dec degrees to hmsdms
# 
def deg2hms(ra_deg, dec_deg):
    ra_hours = ra_deg/360*24
    ra_minutes = (ra_hours-int(ra_hours))*60
    ra_seconds = (ra_minutes-int(ra_minutes))*60
    dec_sign = np.sign(dec_deg)
    dec_minutes = (dec_deg/dec_sign-int(dec_deg/dec_sign))*60
    dec_seconds = (dec_minutes-int(dec_minutes))*60
    return 'J2000 {:02d}h{:02d}m{:.4f}s {:+02d}d{:02d}m{:.4f}s'.format(
            int(ra_hours), int(ra_minutes), ra_seconds, 
            int(dec_sign)*int(dec_deg/dec_sign), int(dec_minutes), dec_seconds)

def hms2deg(ra_dec_str):
    if ra_dec_str.startswith('J2000'):
        ra_dec_str = ra_dec_str.replace('J2000', '')
    ra_str, dec_str = ra_dec_str.strip().split()
    ra_match = re.match(r'^([0-9]+)[h:]([0-9]+)[m:]([0-9.]+)[s]*$', ra_str)
    dec_match = re.match(r'^([+-]*)([0-9]+)[d:]([0-9]+)[m:]([0-9.]+)[s]*$', dec_str)
    if not ra_match or not dec_match:
        raise ValueError('Error! Could not convert ra dec hms to deg {}'.format(ra_dec_str))
    ra_deg = (float(ra_match.group(1)) + float(ra_match.group(2))/60.0 + float(ra_match.group(3))/3600.0) / 24 * 360
    dec_deg = float(dec_match.group(1)+'1') * (
                float(dec_match.group(2)) + float(dec_match.group(3))/60.0 + float(dec_match.group(4))/3600.0)
    return ra_deg, dec_deg



# 
# Find fields that contains a given coordinate within the primary beam
# 
def extract_field_selections(
        vis, 
        target_ra_dec, 
        field_name = None, 
        pb_factor = 1.6, # how many factor of FWHM to consider as the matching circle diameter
        pb_fwhm = None, # arcsec
        ant_diam = 12.0, # meters
        return_field_mosaic_dict = False, 
        return_string = True, 
        verbose = False, 
    ):
    # using global tb
    field_mosaic = {}
    field_mosaic['id'] = []
    field_mosaic['name'] = []
    field_mosaic['ra'] = []
    field_mosaic['dec'] = []
    field_mosaic['reference_ra'] = []
    field_mosaic['reference_dec'] = []
    field_mosaic['phase_ra'] = []
    field_mosaic['phase_dec'] = []
    field_mosaic['delay_ra'] = []
    field_mosaic['delay_dec'] = []
    field_mosaic['time'] = []
    tb.open(os.path.join(vis, 'FIELD'), nomodify=True)
    for i in range(tb.nrows()):
        name = tb.getcell('NAME', i)
        if field_name is not None:
            if name != field_name:
                continue
        reference_dir = tb.getcell('REFERENCE_DIR', i)
        phase_dir = tb.getcell('PHASE_DIR', i)
        delay_dir = tb.getcell('DELAY_DIR', i)
        reference_dec = np.rad2deg(reference_dir[1][0])
        reference_ra = np.rad2deg(reference_dir[0][0])
        if reference_ra < 0.0:
            reference_ra += 360.0
        phase_dec = np.rad2deg(phase_dir[1][0])
        phase_ra = np.rad2deg(phase_dir[0][0])
        if phase_ra < 0.0:
            phase_ra += 360.0
        delay_dec = np.rad2deg(delay_dir[1][0])
        delay_ra = np.rad2deg(delay_dir[0][0])
        if delay_ra < 0.0:
            delay_ra += 360.0
        ra, dec = delay_ra, delay_dec
        field_mosaic['id'].append(i)
        field_mosaic['name'].append(name)
        field_mosaic['ra'].append(ra)
        field_mosaic['dec'].append(dec)
        field_mosaic['reference_ra'].append(reference_ra)
        field_mosaic['reference_dec'].append(reference_dec)
        field_mosaic['phase_ra'].append(phase_ra)
        field_mosaic['phase_dec'].append(phase_dec)
        field_mosaic['delay_ra'].append(delay_ra)
        field_mosaic['delay_dec'].append(delay_dec)
        field_mosaic['time'].append(tb.getcell('TIME', i))
    tb.close()
    if verbose:
        #print('field_mosaic: \n{}'.format(pprint.pformat(field_mosaic, indent=4)))
        print('field_mosaic: {}'.format(field_mosaic))
    # 
    # compute pb
    if pb_fwhm is None:
        spw_dict = aU.getScienceSpws(vis, intent='OBSERVE_TARGET#ON_SOURCE', returnFreqRanges=True)
        #print('spw_dict.values()', spw_dict.values())
        min_freq = np.min(np.array(list(spw_dict.values())).ravel())
        pb_fwhm = 1.14*1.22*(3e8/min_freq)*3600*180/(ant_diam*3.1415926)
    match_radius = pb_factor * pb_fwhm / 2.0 # arcsec
    if verbose:
        print('match_radius: {} arcsec'.format(match_radius))
    # 
    # parse target ra dec
    if str(target_ra_dec).find('m')>0 or str(target_ra_dec).find(':')>0:
        target_ra, target_dec = hms2deg(target_ra_dec)
    else:
        target_ra, target_dec = target_ra_dec
    # 
    # find valid fields that contains the given target
    valid_fields = []
    for i in range(len(field_mosaic['id'])):
        field_ra, field_dec = field_mosaic['ra'][i], field_mosaic['dec'][i]
        offset = np.sqrt(((field_ra-target_ra)*np.cos(np.deg2rad(target_dec)))**2 + (field_dec-target_dec)**2) * 3600.0
        if verbose:
            print('Checking field {} ra dec {:.7f} {:.7f} vs target {:.7f} {:.7f}, offset {:.1f}, match radius {:.1f}'.format(
                   i, field_ra, field_dec, target_ra, target_dec, offset, match_radius))
        if offset < match_radius:
            valid_fields.append(str(field_mosaic['id'][i]))
    if verbose:
        print('Valid fields: {}'.format(valid_fields))
    if return_string:
        valid_fields = ','.join(valid_fields)
    if return_field_mosaic_dict:
        return valid_fields, field_mosaic
    return valid_fields



# 
# Function to print parameter dict
# 
def params2str(params):
    parstr = ''
    for key in params:
        if parstr != '':
            parstr += ', '
        if isinstance(params[key], (str, np.str_)):
            parstr += "{}='{}'".format(key, params[key])
        else:
            parstr += '{}={}'.format(key, params[key])
    return parstr



# 
# Function to collapse index array to channel selecting string
# Example: [0,1,2,3,4,10,11,12,13,20] -> '0~4;10~13;20'
# 
def indexarray2selectingstr(index_array):
    block_indices = np.argwhere(np.diff(index_array)>1).ravel()
    if len(block_indices) == 0:
        selecting_str = '{}~{}'.format(index_array[0], index_array[-1])
    else:
        block_ends = block_indices
        block_starts = block_ends + 1
        if block_starts[0] != 0:
            block_starts = np.concatenate([[0], block_starts])
        if block_ends[-1] != len(index_array)-1:
            block_ends = np.concatenate([block_ends, [len(index_array)-1]])
        # 
        selecting_str = ''
        for k in range(len(block_starts)):
            if k > 0:
                selecting_str += ';'
            selecting_str += '{}~{}'.format(index_array[block_starts[k]], index_array[block_ends[k]])
    return selecting_str
    

# 
# Known (sub)mm strong lines
# 
KNOWN_LINE_DICT = {
    # Band 3
    'sio54': 217.10498e9,
    'co21': 230.53800e9, 
    '13co21': 220.39868e9, 
    'c18o21': 219.56035e9, 
    'h30alpha': 231.900928e9, 
    # Band 8
    'co43': 461.04077e9, 
    'ci10': 492.16068e9, 
    # 
    # 'co65': 691.47308,
    # 'co54': 576.26793,
    # 'co43': 461.04077,
    # 'co32': 345.79599,
    # 'co21': 230.53800,
    # 'co10': 115.27120,
    # '13co65': 661.06728,
    # '13co54': 550.92629,
    # '13co43': 440.76517,
    # '13co32': 330.58797,
    # '13co21': 220.39868,
    # '13co10': 110.20135,
    # 'c18o65': 658.55328,
    # 'c18o54': 548.83101,
    # 'c18o43': 439.08877,
    # 'c18o32': 329.33055,
    # 'c18o21': 219.56035,
    # 'c18o10': 109.78217,
    # 'hcn10': 88.63185,  # J=1-0, F=2-1
    # 'hcn21': 177.26111,  # J=2-1, F=2-1
    # 'hcn32': 265.88618,
    # 'hcn43': 354.50548,
    # 'hcn54': 443.11616,
    # 'hcn65': 531.71639,
    # 'hcn76': 620.30410,
    # 'h13cn10': 86.33992140,
    # 'h13cn21': 172.67785120,
    # 'h13cn32': 259.01179760,
    # 'h13cn43': 345.33976930,
    # 'h13cn54': 431.65977480,
    # 'h13cn65': 517.96982100,
    # 'h13cn76': 604.26791400,
    # 'cs10': 48.99095,
    # 'cs21': 97.98095,
    # 'cs32': 146.96903,
    # 'cs43': 195.95421,
    # 'cs54': 244.93556,
    # 'cs65': 293.91209,
    # 'cs76': 342.88285,
    # 'cs87': 391.84689,
    # 'cs98': 440.80323,
    # 'cs109': 489.75092,
    # 'cs1110': 538.68900,
    # 'cs1211': 587.61649,
    # 'cs1312': 636.53246,
    # 'cs1413': 685.43592,
    # '13cs10': 46.24756320,
    # '13cs21': 92.49430800,
    # '13cs32': 138.73933500,
    # '13cs43': 184.98177200,
    # '13cs54': 231.22068520,
    # '13cs65': 277.45540500,
    # '13cs76': 323.68497300,
    # '13cs87': 369.90855050,
    # '13cs98': 416.12527510,
    # '13cs109': 462.33429010,
    # '13cs1110': 508.53473910,
    # '13cs1211': 554.72576570,
    # '13cs1312': 600.90648000,
    # '13cs1413': 647.07615000,
    # 'hcop10': 89.18852,
    # 'hcop21': 178.37506,
    # 'hcop32': 267.55763,
    # 'hcop43': 356.73422,
    # 'hcop54': 445.90287,
    # 'hcop65': 535.06158,
    # 'hcop76': 624.20836,
    # 'h13cop10': 86.75428840,
    # 'h13cop21': 173.50670030,
    # 'h13cop32': 260.25533900,
    # 'h13cop43': 346.99834400,
    # 'h13cop54': 433.73383270,
    # 'h13cop65': 520.45988430,
    # 'h13cop76': 607.17464560,
    # 'hnc10': 90.66357,
    # 'hnc21': 181.32476,
    # 'hnc32': 271.98114,
    # 'hnc43': 362.63030,
    # 'hnc54': 453.26992,
    # 'hnc65': 543.89755,
    # 'hnc76': 634.51083,
    # 'hnco43': 87.925237,  # 4( 0, 4)- 3( 0, 3)
    # 'hnco54': 109.905749,  # 5( 0, 5)- 4( 0, 4)
    # 'hn13c10': 87.09085000,
    # 'hn13c21': 174.17940800,
    # 'hn13c32': 261.26331010,
    # 'hn13c43': 348.34026950,
    # 'hn13c54': 435.40796260,
    # 'hn13c65': 522.46407300,
    # 'hn13c76': 609.50628400,
    # 'ci10': 492.16065,  # 3P1-3P0
    # 'ci21': 809.34197,  # 3P2-3P1
    # 'sio10': 43.42376,
    # 'sio21': 86.84696,
    # 'sio32': 130.26861,
    # 'sio43': 173.68831,
    # 'sio54': 217.10498,
    # 'sio65': 260.51802,
    # 'sio76': 303.92696,
    # 'sio87': 347.33063,
    # 'sio98': 390.72845,
    # 'sio109': 434.11955,
    # 'sio1110': 477.50310,
    # 'sio1211': 520.87820,
    # 'sio1312': 564.24396,
    # 'sio1413': 607.59942,
    # 'sio1514': 650.94359,
    # 'sio1615': 694.27543,
    # 'hi21cm': 1.420405751,
    # 'h19alpha': 888.047022,
    # 'h21alpha': 662.404162,
    # 'h24alpha': 447.540278,
    # 'h25alpha': 396.900834,
    # 'h26alpha': 353.622747,
    # 'h27alpha': 316.415425,
    # 'h28alpha': 284.250571,
    # 'h29alpha': 256.302035,
    # 'h30alpha': 231.900928,
    # 'h31alpha': 210.501771,
    # 'h32alpha': 191.656728,
    # 'h33alpha': 174.995805,
    # 'h34alpha': 160.211511,
    # 'h35alpha': 147.046878,
    # 'h36alpha': 135.286032,
    # 'h38alpha': 115.274399,
    # 'h39alpha': 106.737357,
    # 'h40alpha': 99.022952,
    # 'h41alpha': 92.034434,
    # 'h42alpha': 85.688390,
    # 'h43alpha': 79.912651,
    # 'h44alpha': 74.644562,
    # 'h45alpha': 69.829551,
    # 'h53alpha': 42.951968,
    # 'h54alpha': 40.630498,
    # 'h55alpha': 38.473358,
    # 'h56alpha': 36.466260,
    # 'h57alpha': 34.596383,
    # 'h58alpha': 32.852196,
}


# 
# Find lines in science spw, return a dictionary of spw selections of the lines and continuum channels
# 
def extract_line_spw_selections(
        vis,
        vlsrk = None, # km/s
        z = None, 
        fwzi = 500.0, # km/s
        line_dict = None, 
        verbose = False, 
    ):
    # 
    # find line spw
    c = 2.99792458e5 # km/s
    if z is None:
        z = 0.0
    if vlsrk is None:
        vlsrk = c*z
    if line_dict is None:
        line_dict = KNOWN_LINE_DICT
    else:
        for key in KNOWN_LINE_DICT:
            if key not in line_dict:
                line_dict[key] = KNOWN_LINE_DICT[key]
    # 
    spw_dict = aU.getScienceSpws(vis, intent='OBSERVE_TARGET#ON_SOURCE', returnFreqRanges=True)
    spw_chan_widths = aU.getScienceSpwChanwidths(vis)
    # 
    science_spws = [str(t) for t in list(spw_dict.keys())]
    continuum_spw_selections = []
    continuum_selections = []
    line_spw_selections = []
    line_selections = []
    individual_line_selections = {}
    for ispw, spw in enumerate(list(spw_dict.keys())):
        freq_start, freq_end = spw_dict[spw]
        chan_width = spw_chan_widths[ispw]
        freq_array = np.arange(freq_start, freq_end+chan_width, chan_width)
        line_mask = np.full(len(freq_array), fill_value=0, dtype=int) # 1 if it is a line channel
        if verbose:
            print('spw {}, freq_array: {} {}'.format(spw, np.min(freq_array), np.max(freq_array)))
        for line_name, line_freq in line_dict.items():
            freq_lower = line_freq / (1 + vlsrk/c) * (1 - fwzi / c)
            freq_upper = line_freq / (1 + vlsrk/c) * (1 + fwzi / c)
            index_lower = np.argwhere(np.abs(freq_array-freq_lower)<1.5*np.abs(chan_width)).ravel()
            index_upper = np.argwhere(np.abs(freq_array-freq_upper)<1.5*np.abs(chan_width)).ravel()
            if len(index_lower) > 0 and len(index_upper) > 0:
                index_lower = index_lower[0]
                index_upper = index_upper[-1]
            elif len(index_lower) > 0 and len(index_upper) == 0:
                index_lower = index_lower[0]
                index_upper = len(freq_array)-1
            elif len(index_lower) == 0 and len(index_upper) > 0:
                index_lower = 0
                index_upper = index_upper[-1]
            else:
                index_lower = 999
                index_upper = 0
            if verbose:
                print('spw {}, line {}, freq {} {}, index {} {}'.format(spw, line_name, freq_lower, freq_upper, index_lower, index_upper))
            if index_lower <= index_upper: # found line in this spw
                line_mask[index_lower:index_upper+1] = 1
                # 
                if line_name not in individual_line_selections:
                    individual_line_selections[line_name] = []
                individual_line_selections[line_name].append('{}:{}~{}'.format(spw, index_lower, index_upper))
        # 
        continuum_indicies = np.argwhere(line_mask==0).ravel()
        #print('spw {}, continuum_indicies: {}'.format(spw, continuum_indicies))
        if len(continuum_indicies) == len(line_mask): # all are continuum channels
            continuum_selections.append('{}'.format(spw))
            continuum_spw_selections.append(str(spw))
        elif len(continuum_indicies) > 0: # select some channels
            continuum_selection = indexarray2selectingstr(continuum_indicies)
            continuum_selections.append('{}:{}'.format(spw, continuum_selection))
            continuum_spw_selections.append(str(spw))
        # 
        line_indicies = np.argwhere(line_mask>0).ravel()
        if verbose:
            print('spw {}, line_indicies: {}'.format(spw, line_indicies))
        if len(line_indicies) == len(line_mask): # all are line channels
            line_selections.append('{}'.format(spw))
            line_spw_selections.append(str(spw))
        elif len(line_indicies) > 0: # select some channels
            line_selection = indexarray2selectingstr(line_indicies)
            line_selections.append('{}:{}'.format(spw, line_selection))
            line_spw_selections.append(str(spw))
    # 
    spw_selections = dict(
        continuum=continuum_selections, 
        continuum_spws=continuum_spw_selections,
        line=line_selections,
        line_spws=line_spw_selections,
        science=science_spws,
        individual_lines=individual_line_selections,
    )
    if verbose:
        print('spw_selections: \n{}'.format(pprint.pformat(spw_selections, indent=4)))
    # 
    return spw_selections



# 
# Flag line data in a uvdata
# 
def flag_line_data(
        vis, 
        outputvis, 
        field, 
        spw, 
        collapse = False, # collapse to make a continuum data set (keeping individual spws)
        contspw = '', # continuum spw, if collapse is True
        width = 3840, 
        verbose = False, 
        overwrite = False, 
    ):
    if verbose:
        print('flag line data with field {!r}, spw {!r}'.format(field, spw))
    if os.path.exists(outputvis): 
        if overwrite:
            print('Caution! Overwriting existing {!r}'.format(outputvis))
            if os.path.exists(outputvis+'.backup'):
                shutil.rmtree(outputvis+'.backup')
            shutil.move(outputvis, outputvis+'.backup')
        else:
            print('Found existing {!r} and overwrite is False. Skipping it.'.format(outputvis))
            return
    if os.path.exists(outputvis+'.tmp'):
        shutil.rmtree(outputvis+'.tmp')
    shutil.copytree(vis, outputvis+'.tmp')
    datacolumn = 'data'
    if 'CORRECTED_DATA' in aU.dataColumns(outputvis+'.tmp'):
        datacolumn = 'corrected'
    # 
    if spw == '':
        # nothing to flag
        print('spw is empty, nothing to flag.')
    else:
        flagdata(vis=outputvis+'.tmp', mode='manual', field=field, spw=spw, datacolumn=datacolumn, flagbackup=False)
    # 
    if collapse:
        if verbose:
            print('flag_line_data')
            print('mstransform with width {} to make a continuum data set'.format(width))
        #mstransform(vis=outputvis+'.tmp', outputvis=outputvis+'.tmpB', width=width, datacolumn=datacolumn, reindex=False)
        #--> mstransform width is not chanbin!
        mstransform(vis=outputvis+'.tmp', outputvis=outputvis+'.tmpB', 
                    chanaverage=True, chanbin=width, 
                    datacolumn=datacolumn, keepflags=False, reindex=False)
        #split(vis=outputvis+'.tmp', outputvis=outputvis+'.tmpB', spw=contspw, width=width, datacolumn=datacolumn)
        #--> split cannot handle channel edge flags well! if some edge channels have flag=1, then the output collapsed uvw flag is also 1!
        #--> bug!! -- not bug, my fault. Now checking if spw=='', do not flagdata.
        shutil.rmtree(outputvis+'.tmp')
        shutil.move(outputvis+'.tmpB', outputvis)
    else:
        shutil.move(outputvis+'.tmp', outputvis)



# 
# Estimate tclean parameters
# 
def estimate_tclean_params(
        vis,
        fov = None,
        field = '', 
        intent = 'OBSERVE_TARGET#ON_SOURCE',
        spw = '',
        maxBaselinePercentile = 95, 
        pblevel = 0.2, 
        npix = 5, 
        cell = None, 
        maximsize = None, 
    ):
    # 
    tclean_params = {}
    if cell is None:
        cell = aU.pickCellSize(vis, 
                               intent=intent, 
                               maxBaselinePercentile=maxBaselinePercentile, 
                               npix=npix, 
                               cellstring=True,
                               ) 
                           # spw=spw, cannot specify spw=spw because aU.pickCellSize issue ('meanfreq = mymsmd.meanfreq(int(spw))')
    else:
        cell = str(cell)
    cellsize = float(cell.replace('arcsec',''))
    # determine fov and imsize
    if fov is None:
        if field and (field != ''):
            if re.match(r'^[0-9,]+$', field):
                field_ints = list(map(int, field.split(',')))
                result = aU.plotmosaic(vis, doplot=False, pblevel=pblevel, intent=intent, spw=spw, field=field_ints, verbose=False) # sourceid is for ephemeris object
            else:
                result = aU.plotmosaic(vis, doplot=False, pblevel=pblevel, intent=intent, spw=spw, sourceid=field, verbose=False) # sourceid is for ephemeris object
            tclean_params['field'] = field
        else: # assuming that this vis data only contains science target fields
            result = aU.plotmosaic(vis, doplot=False, pblevel=pblevel, intent=intent, spw=spw, verbose=False)
        centralField, raMax, raMin, decMax, decMin = result
        # ra dec Min Max are in arcsec
        imsize = [int(np.ceil(abs(raMax-raMin)/cellsize)), # no need `*np.cos(np.deg2rad((decMax+decMin)/2.0))`
                  int(np.ceil(abs(decMax-decMin)/cellsize))]
        imsize = [int(aU.getOptimumSize(imsize[0])), int(aU.getOptimumSize(imsize[1]))]
    else:
        fovpixsize = int(np.ceil(fov/cellsize))
        imsize = int(aU.getOptimumSize(fovpixsize))
        imsize = [imsize, imsize]
    if maximsize is not None:
        maximsize = int(maximsize)
        if imsize[0] > maximsize:
            imsize[0] = maximsize
        if imsize[1] > maximsize:
            imsize[1] = maximsize
    if spw and (spw != ''):
        tclean_params['spw'] = spw
    tclean_params['cell'] = cell
    tclean_params['imsize'] = imsize
    tclean_params['specmode'] = 'mfs'
    tclean_params['gridder'] = 'mosaic'
    tclean_params['deconvolver'] = 'hogbom'
    return tclean_params



# 
# Selfcal continuum
# 
def selfcal_continuum_data(
        vis, # the input vis must be a continuum data
        outputvis, 
        field = '', 
        phasecenter = '', 
        intent = 'OBSERVE_TARGET#ON_SOURCE', 
        spw = None, 
        cell = None, 
        maximsize = None, 
        calmode = 'p', 
        solint = '120s', 
        gaintable = [], 
        verbose = True, 
        overwrite = False, 
        showgui = True, 
        cleanup = True, 
    ):
    # 
    # check input field and phasecenter
    if field == '' and phasecenter == '':
        raise Exception('Error! Please provide either field or phasecenter when calling selfcal_continuum_data!')
    # 
    # prepare result dict
    outputbasename = outputvis # re.sub(r'\.(ms|ms.split|ms.split.cal)$', r'', outputvis)
    caltable = outputbasename + '.pcal'
    result_dict = {
        'outputvis': outputvis, 
        'caltable': caltable, 
        'image_before_selfcal': outputbasename+'.clean.before.selfcal.image',
        'image_after_selfcal': outputbasename+'.clean.before.selfcal.image',
        'model_for_selfcal': outputbasename+'.clean.before.selfcal.model',
        'pb_image': outputbasename+'.clean.after.selfcal.pb',
    }
    # 
    # check existing outputvis
    if os.path.exists(outputvis): 
        if overwrite:
            print('Caution! Overwriting existing {!r}'.format(outputvis))
            if os.path.exists(outputvis+'.backup'):
                shutil.rmtree(outputvis+'.backup')
            shutil.move(outputvis, outputvis+'.backup')
        else:
            print('Found existing {!r} and overwrite is False. Skipping it.'.format(outputvis))
            return result_dict
    #if os.path.exists(outputvis+'.tmp'):
    #    shutil.rmtree(outputvis+'.tmp')
    os.system('rm -rf {}'.format(outputvis+'.tmp*'))
    if verbose:
        print('Copying {} -> {} for selfcal processing'.format(vis, outputvis+'.tmp'))
    shutil.copytree(vis, outputvis+'.tmp')
    original_vis = vis
    vis = outputvis+'.tmp'
    # 
    # find fields containing the phasecenter
    if phasecenter != '' and field == '':
        if verbose:
            print('Extracting field selections:')
        field = extract_field_selections(vis, phasecenter, pb_factor=1.3) # np.exp(-0.5*(1.3*2.35482/2)**2) = 0.3
        if verbose:
            print('  {}'.format(field))
    # 
    # make continuum dataset
    if spw is None:
        spw = aU.getScienceSpws(vis, intent=intent, returnString=True)
    mstransform_params = dict(field=field, intent=intent, spw=spw, 
                              chanaverage=True, chanbin=3840, # width=3840, 
                              timeaverage=True, timebin='30s', # must set timeaverage=True
                              datacolumn='data', keepflags=False, reindex=False)
    mstransform(vis, outputvis+'.tmp.cont', **mstransform_params)
    # 
    # copy FIELD/EPHEM0_*.tab
    if (len(glob.glob('{}/FIELD/EPHEM0_*.tab'.format(vis))) > 0) and \
       (len(glob.glob('{}/FIELD/EPHEM0_*.tab'.format(outputvis+'.tmp.cont'))) == 0):
        os.system('cp -r {}/FIELD/EPHEM0_*.tab {}/FIELD/'.format(vis, outputvis+'.tmp.cont'))
    # 
    vis = outputvis+'.tmp.cont'
    # 
    # make initial clean
    if verbose:
        print('Estimating tclean params:')
    tclean_params = estimate_tclean_params(vis, intent=intent, field=field, spw=spw, cell=cell, maximsize=maximsize) # cannot specify spw=spw because aU.pickCellSize issue ('meanfreq = mymsmd.meanfreq(int(spw))')
    if phasecenter != '':
        tclean_params['phasecenter'] = phasecenter
    if verbose:
        print(pprint.pformat(tclean_params, indent=4))
    dirty_image = outputvis+'.tmp.dirty'
    if verbose:
        print('Making dirty image:')
        print('tclean({})'.format(params2str(dict(vis=vis, imagename=dirty_image, niter=0, datacolumn='data', **tclean_params))))
    tclean(vis=vis, imagename=dirty_image, niter=0, datacolumn='data', **tclean_params)
    # 
    dirty_stats = imstat(dirty_image+'.image')
    with open(dirty_image+'.image.stats.json', 'w') as fp:
        json.dump(dirty_stats, fp, indent=4, cls=NpEncoder)
    if verbose:
        print('Dirty image statistics:')
        print(pprint.pformat(dirty_stats, indent=4))
    # 
    rms = dirty_stats['rms'][0]
    maxposf = dirty_stats['maxposf']
    phasecenterf = 'J2000 ' + ' '.join(maxposf.split(',')[0:2])
    tclean_params['threshold'] = 4.0*rms
    clean_image = outputvis+'.tmp.clean'
    if verbose:
        print('Making clean image:')
        print('tclean({})'.format(params2str(dict(vis=vis, imagename=clean_image, niter=100, datacolumn='data', 
                                                  savemodel='none', **tclean_params))))
    tclean(vis=vis, imagename=clean_image, niter=100, datacolumn='data', 
           savemodel='none', **tclean_params)
    if verbose:
        print('Writing clean model:')
        print('tclean({})'.format(params2str(dict(vis=vis, imagename=clean_image, niter=0, datacolumn='data', savemodel='modelcolumn', 
                                                  calcpsf=False, calcres=False, restoration=False, **tclean_params))))
    tclean(vis=vis, imagename=clean_image, niter=0, datacolumn='data', savemodel='modelcolumn', 
           calcpsf=False, calcres=False, restoration=False, **tclean_params)
    # 
    # gaincal
    #caltable = outputbasename + '.pcal'
    if gaintable is None or gaintable == '':
        gaintable = []
    elif isinstance(gaintable, str):
        gaintable = [gaintable]
    if verbose:
        print('Making calibration table:')
        print('gaincal({})'.format(params2str(dict(vis=vis, caltable=caltable, gaintable=gaintable, calmode=calmode, solint=solint))))
    gaincal(vis=vis, caltable=caltable, gaintable=gaintable, calmode=calmode, solint=solint)
    if not os.path.exists(caltable):
        raise Exception('Error! Failed to produce calibration table: {}'.format(caltable))
    # 
    # plotcal
    if verbose:
        print('Plotting calibration table')
    plotms(caltable, xaxis='time', yaxis='phase', plotrange=[0,0,-180,180], iteraxis='spw', gridrows=2, gridcols=2, 
           plotfile=caltable+'.plot.phase.vs.time.iter.spw.png', overwrite=True, highres=True, showgui=showgui)
    plotms(caltable, xaxis='time', yaxis='phase', plotrange=[0,0,-180,180], iteraxis='antenna', gridrows=3, gridcols=3, 
           plotfile=caltable+'.plot.phase.vs.time.iter.antenna.png', overwrite=True, highres=True, showgui=showgui)
    # 
    # applycal
    if verbose:
        print('Applying cal table to the extracted continuum dataset')
    applycal(vis=vis, gaintable=gaintable+[caltable], flagbackup=False)
    # 
    # reclean
    reclean_image = outputvis+'.tmp.reclean'
    if verbose:
        print('Remaking clean image:')
        print('tclean({})'.format(params2str(dict(vis=vis, imagename=reclean_image, niter=100, datacolumn='corrected', 
                                                  savemodel='none', pbcor=True, **tclean_params))))
    tclean(vis=vis, imagename=reclean_image, niter=100, datacolumn='corrected', 
           savemodel='none', pbcor=True, **tclean_params)
    # 
    # applycal to outputvis
    if verbose:
        print('Applying cal table to the output dataset')
    applycal(vis=outputvis+'.tmp', gaintable=gaintable+[caltable], flagbackup=False) # spwmap=spw
    # 
    # rename
    shutil.move(outputvis+'.tmp', outputvis)
    shutil.move(clean_image+'.model', outputbasename+'.clean.before.selfcal.model')
    shutil.move(clean_image+'.image', outputbasename+'.clean.before.selfcal.image')
    shutil.move(reclean_image+'.image', outputbasename+'.clean.after.selfcal.image')
    shutil.move(reclean_image+'.pb', outputbasename+'.clean.after.selfcal.pb')
    # 
    # cleanup
    if cleanup:
        print('Cleaning up {}'.format(outputvis+'.tmp*'))
        os.system('rm -rf {}'.format(outputvis+'.tmp*'))
    # 
    # return_dict
    result_dict = {
        'outputvis': outputvis, 
        'caltable': caltable, 
        'image_before_selfcal': outputbasename+'.clean.before.selfcal.image',
        'image_after_selfcal': outputbasename+'.clean.before.selfcal.image',
        'model_for_selfcal': outputbasename+'.clean.before.selfcal.model',
        'pb_image': outputbasename+'.clean.after.selfcal.pb',
    }
    return result_dict



# 
# timebin a uv data
# 
def timebin_uvdata(
        ms_file, 
        timebin='30s',
    ):
    old_size = get_dir_size(ms_file)
    os.system('touch "{}"'.format(ms_file+'.touch'))
    # get data column
    tb.open(ms_file)
    colnames = tb.colnames()
    tb.close()
    if 'CORRECTED_DATA' in colnames:
        datacolumn = 'corrected'
    else:
        datacolumn = 'data'

    # split with timebin
    split(vis=ms_file, outputvis=ms_file+'.timebinned', timebin=timebin, datacolumn=datacolumn, keepflags=False)
    if not os.path.exists(ms_file+'.original'):
        os.system('mv "{}" "{}"'.format(ms_file, ms_file+'.original'))
    else:
        raise Exception('Error! File exists: {}'.format(ms_file+'.original'))
    if not os.path.exists(ms_file):
        os.system('mv "{}" "{}"'.format(ms_file+'.timebinned', ms_file))
    else:
        raise Exception('Error! File exists: {}'.format(ms_file))
    if not os.path.exists(ms_file+'.timebin.log'):
        os.system('mv "{}" "{}"'.format(ms_file+'.touch', ms_file+'.timebin.log'))
        new_size = get_dir_size(ms_file)
        with open(ms_file+'.timebin.log', 'a') as fp:
            fp.write('timebin: 30s\n')
            fp.write('old_size: {:6g} GB\n'.format(float(old_size)/1024/1024/1024))
            fp.write('new_size: {:6g} GB\n'.format(float(new_size)/1024/1024/1024))
    else:
        raise Exception('Error! File exists: {}'.format(ms_file+'.timebin.log'))
    
    if os.path.exists(ms_file+'.touch'):
        print('Failed to split the ms data {}? Please remove the touch file and re-run the script.'.format(ms_file))
    else:
        print('Successfully time-binned the ms data {}'.format(ms_file))
        os.system('rm -rf "{}"'.format(ms_file+'.original'))






