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

import os, sys, re, json, copy, time, datetime, shutil
import numpy as np
import analysisUtils as aU
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
for iscope, scope in enumerate(inspect.stack()):
    print('inspect.stack()', iscope, scope.filename)
    if scope.filename.find('importlib')>=0 or scope.filename.find('<')>=0:
        continue
    if iscope == 0:
        continue
    caller_globals = scope.frame.f_globals
    break
#print('caller_globals', caller_globals.keys())
if check_casa(caller_globals):
    for key in caller_globals:
        if key not in globals(): 
            if key.find('casa')>=0 or key in ['tb', 'ia', 'flagdata', 'split', 'mstransform', 'tclean', 'imstat', 'gaincal', 'plotms', 'applycal']:
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
    'co21': 230.53800e9, 
    '13co21': 220.39868e9, 
    'c18o21': 219.56035e9, 
    'h30alpha': 231.900928e9, 
    'ci10': 492.16068e9, 
}


# 
# Find lines in science spw, return a dictionary of spw selections of the lines and continuum channels
# 
def extract_line_spw_selections(
        vis,
        vlsrk = 0.0, # km/s
        z = 0.0, 
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
        verbose=False, 
        overwrite=False,
    ):
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
    if verbose:
        print('flagdata with field {!r}, spw {!r}'.format(field, spw))
    flagdata(vis=outputvis+'.tmp', mode='manual', field=field, spw=spw, flagbackup=False)
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
    ):
    # 
    tclean_params = {}
    cell = aU.pickCellSize(vis, intent=intent, spw=spw, maxBaselinePercentile=maxBaselinePercentile, npix=npix, cellstring=True)
    cellsize = float(cell.replace('arcsec',''))
    # determine fov and imsize
    if fov is None:
        if field and (field != ''):
            result = aU.plotmosaic(vis, doplot=False, pblevel=pblevel, intent=intent, spw=spw, sourceid=field, verbose=False)
        else: # assuming that this vis data only contains science target fields
            result = aU.plotmosaic(vis, doplot=False, pblevel=pblevel, intent=intent, spw=spw, verbose=False)
        centralField, raMax, raMin, decMax, decMin = result
        imsize = [int(np.ceil(abs(raMax-raMin)*np.cos(np.deg2rad((decMax+decMin)/2.0))/cellsize)), 
                  int(np.ceil(abs(decMax-decMin)/cellsize))]
        imsize = [int(aU.getOptimumSize(imsize[0])), int(aU.getOptimumSize(imsize[1]))]
    else:
        fovpixsize = int(np.ceil(fov/cellsize))
        imsize = int(aU.getOptimumSize(fovpixsize))
        imsize = [imsize, imsize]
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
        spw = '', 
        gaintable = [], 
        verbose = False, 
        overwrite = False, 
        cleanup = True, 
    ):
    # 
    # check input field and phasecenter
    if field == '' and phasecenter == '':
        raise Exception('Error! Please provide either field or phasecenter when calling selfcal_continuum_data!')
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
            return
    #if os.path.exists(outputvis+'.tmp'):
    #    shutil.rmtree(outputvis+'.tmp')
    os.system('rm -rf {}'.format(outputvis+'.tmp*'))
    if verbose:
        print('Copying {} -> {} for selfcal processing'.format(vis, outputvis+'.tmp'))
    shutil.copytree(vis, outputvis+'.tmp')
    original_vis = vis
    vis = outputvis+'.tmp'
    outputbasename = re.sub(r'\.(ms|ms.split|ms.split.cal)$', r'', outputvis)
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
    spw = aU.getScienceSpws(vis, intent=intent, returnString=True)
    mstransform_params = dict(field=field, intent=intent, spw=spw, width=3840, timebin='30s', datacolumn='data', reindex=False)
    mstransform(vis, outputvis+'.tmp.cont', **mstransform_params)
    vis = outputvis+'.tmp.cont'
    # 
    # make initial clean
    if verbose:
        print('Estimating tclean params:')
    tclean_params = estimate_tclean_params(vis, intent=intent) # , field=field, spw=spw
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
    caltable = outputbasename + '.pcal'
    if gaintable is None or gaintable == '':
        gaintable = []
    elif isinstance(gaintable, str):
        gaintable = [gaintable]
    gaincal(vis=vis, caltable=caltable, gaintable=gaintable, calmode='p', solint='120s')
    # 
    # plotcal
    plotms(caltable, xaxis='time', yaxis='phase', plotrange=[0,0,-180,180], iteraxis='spw', gridrows=2, gridcols=2, 
           plotfile=caltable+'.plot.phase.vs.time.iter.spw.png', overwrite=True, highres=True)
    plotms(caltable, xaxis='time', yaxis='phase', plotrange=[0,0,-180,180], iteraxis='antenna', gridrows=3, gridcols=3, 
           plotfile=caltable+'.plot.phase.vs.time.iter.antenna.png', overwrite=True, highres=True)
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
        print('tclean({})'.format(params2str(dict(vis=vis, imagename=reclean_image, niter=100, datacolumn='corrected', savemodel='none', 
                                                  pbcor=True, **tclean_params))))
    tclean(vis=vis, imagename=reclean_image, niter=100, datacolumn='corrected', savemodel='none', 
           pbcor=True, **tclean_params)
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
        os.system('rm -rf {}'.format(outputvis+'.tmp*'))







