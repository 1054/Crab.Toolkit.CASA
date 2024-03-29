#!/usr/bin/env python
# 
# This needs to be run in CASA
# 
# CASA modules/functions used:
#     tb, casalog, mstransform, inp, saveinputs, exportfits, tclean
# 
# Example:
#     import a_dzliu_code_level_4_clean; reload(a_dzliu_code_level_4_clean); from a_dzliu_code_level_4_clean import dzliu_clean; dzliu_clean()
# 
from __future__ import print_function
import os, sys, re, json, copy, timeit, shutil
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
    from imstat_cli import imstat_cli_
    imstat = imstat_cli_()
else:
    # see CASA 6 updates here: https://alma-intweb.mtk.nao.ac.jp/~eaarc/UM2018/presentation/Nakazato.pdf
    from casatasks import tclean, mstransform, exportfits, imstat
    #from casatasks import sdbaseline
    #from casatools import ia



# 
# def print2
# 
def print2(message):
    print(message)
    casalog.post(message, 'INFO')







# 
# def get_optimized_imsize
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
# def get field phasecenters
# 
def get_field_phasecenters(vis, galaxy_name):
    # 
    # Requires CASA module/function tb.
    # 
    casalog.origin('get_field_phasecenters')
    # 
    tb.open(vis+os.sep+'FIELD')
    field_names = tb.getcol('NAME')
    field_phasecenters = [tb.getcell('DELAY_DIR', i) for i in range(tb.nrows())] # rad,rad
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
# def find_lab_line_name_and_freq
# 
def find_lab_line_name_and_freq(line_name):
    # 
    # No CASA module/function used here.
    # 
    lab_line_dict = {}
    lab_line_dict['HI21cm'] = {'rest-freq': 1.420405751}
    lab_line_dict['CO(1-0)'] = {'rest-freq': 115.2712018}
    lab_line_dict['CO(2-1)'] = {'rest-freq': 230.5380000}
    line_name_cleaned = re.sub(r'[^a-zA-Z0-9]', r'', line_name).lower()
    found_lab_line_name = None
    found_lab_line_freq = None
    for lab_line_name in lab_line_dict.keys():
        lab_line_name_cleaned = re.sub(r'[^a-zA-Z0-9]', r'', str(lab_line_name)).lower()
        #print('line_name_cleaned', line_name_cleaned, 'lab_line_name_cleaned', lab_line_name_cleaned)
        if line_name_cleaned == lab_line_name_cleaned:
            found_lab_line_name = lab_line_name
            found_lab_line_freq = lab_line_dict[lab_line_name]['rest-freq']
    if found_lab_line_name is None:
        raise Exception('Error! Could not find line name "%s" in our lab line dict! Please update the code!'%(line_name))
    return found_lab_line_name, found_lab_line_freq





# 
# def split_line_visibilities
# 
def split_line_visibilities(dataset_ms, output_ms, galaxy_name, line_name, line_velocity, line_velocity_width, line_velocity_resolution):
    # 
    # Requires CASA module/function tb, default, inp, saveinputs, mstransform.
    # 
    casalog.origin('split_line_visibilities')
    
    # 
    # check existing file
    # 
    if os.path.isdir(output_ms):
        print2('Found existing "%s"! Will not overwrite it! Skipping split_line_visibilities()!'%(output_ms))
        return
    
    # 
    # calc linefreq
    # 
    lab_line_name, lab_line_freq = find_lab_line_name_and_freq(line_name)
    #linefreq = lab_line_freq*(1.0-(line_velocity/2.99792458e5))*1e9 # Hz, for velocity with radio definition
    linefreq = lab_line_freq/(1.0+(line_velocity/2.99792458e5))*1e9 # Hz, for velocity with optical definition
    
    # 
    # check data column
    # 
    tb.open(dataset_ms)
    if 'CORRECTED_DATA' in tb.colnames():
        datacolumn = 'CORRECTED'
    else:
        datacolumn = 'DATA'
    tb.close()
    
    # 
    # get spectral_window and ref_freq
    # 
    tb.open(dataset_ms+os.sep+'SPECTRAL_WINDOW')
    spw_names = tb.getcol('NAME')
    spw_chan_freq_col = [tb.getcell('CHAN_FREQ', i) for i in range(tb.nrows())]
    spw_chan_width_col = [tb.getcell('CHAN_WIDTH', i) for i in range(tb.nrows())]
    spw_ref_freq_col = tb.getcol('REF_FREQUENCY')
    tb.close()
    
    #print('spw_names:', spw_names)
    valid_spw_indicies = np.argwhere([re.match(r'.*FULL_RES.*', t) is not None for t in spw_names]).flatten().tolist()
    if len(valid_spw_indicies) == 0:
        valid_spw_indicies = np.argwhere([re.match(r'.*WVR.*', t) is None for t in spw_names]).flatten().tolist()
    ref_freq_Hz = np.nan
    linespws = [] # the target line may be contained in multiple spws?
    linechanwidths = [] # a list of channel widths in Hz for each matched line spw, for all matched line spws
    #linechanfreq = []
    #linechanfreqs = []
    for i in valid_spw_indicies:
        spw_chan_freq_list = spw_chan_freq_col[i]
        spw_chan_width_list = spw_chan_width_col[i]
        spw_ref_freq = spw_ref_freq_col[i]
        print2('spw_%d, ref_freq %.3e Hz, chan_freq %.3e .. %.3e Hz (%d), chan_width %.3e Hz'%(i, spw_ref_freq, np.max(spw_chan_freq_list), np.min(spw_chan_freq_list), len(spw_chan_freq_list), np.min(spw_chan_width_list) ) )
        # find the target line in these spw
        if linefreq >= np.min(spw_chan_freq_list) and linefreq <= np.max(spw_chan_freq_list):
            # found our target line within this spw
            ref_freq_Hz = spw_ref_freq
            linespws.append(i)
            linechanwidths.append(spw_chan_width_list) # append all valid spw 
            #linechanfreq = spw_chan_freq_list
            #linechanfreqs.append(spw_chan_freq_list)
    
    if len(linespws) == 0:
        raise ValueError('Error! Target line %s at rest-frame %.3f GHz was not covered by the "SPECTRAL_WINDOW" of the input vis "%s"!'%(lab_line_name, lab_line_freq, dataset_ms))
    
    # if we have multiple spw and some linechanwidth do not match, then we only use the best linechanwidth spw
    valid_linespws_indicies = []
    linechanwidth = np.min([np.min(t) for t in linechanwidths])
    for i in range(len(linespws)):
        linechanwidth_i = np.min(linechanwidths[i])
        if np.isclose(linechanwidth_i, linechanwidth):
            valid_linespws_indicies.append(i)
        else:
            print2('Warning! Discared spw %d due to coarse channel width of %.3e Hz than the finest channel width of %.3e Hz.'%(i, linechanwidth_i, linechanwidth))
    
    linespws = [linespws[i] for i in valid_linespws_indicies]
    linechanwidths = [linechanwidths[i] for i in valid_linespws_indicies]
    #linechanfreq = linechanfreq[valid_linespws_indicies]
    #linechanfreqs = linechanfreqs[valid_linespws_indicies]
    
    # make numpy array
    linespws = np.array(linespws)
    
    # convert channel width from frequency to velocity
    linechanwidth = np.min([np.min(t) for t in linechanwidths])
    linechanwidth_kms = linechanwidth / ref_freq_Hz * 2.99792458e5 # km/s
    #linechanwidths_kms = linechanwidths / ref_freq_Hz * 2.99792458e5 # km/s
    
    print2('linespws = %s'%(linespws))
    print2('linechanwidth = %.3e Hz'%(linechanwidth))
    print2('linechanwidth_kms = %s'%(linechanwidth_kms))
    
    # 
    # chanbin
    # 
    #chanbin = line_velocity_resolution / np.min(linechanwidth_kms)
    #chanbin = int(np.round(chanbin))
    #print2('chanbin = %s'%(chanbin))
    #chanaverage = True
    #chanbin = chanbin
    
    # 
    # nchan and start and width
    # 
    regridms = True
    mode = 'frequency'
    if line_velocity_resolution > 0:
        # if the user has input a positive line_velocity_resolution, then we convert it to channel width factor
        width_channel_number = line_velocity_resolution / np.abs(linechanwidth_kms[0])
        width_channel_number = int(np.round(width_channel_number)) # in units of channel number, like a rebin factor
    else:
        width_channel_number = 1 # otherwise keep the original channel width
    width_freq_Hz = width_channel_number * linechanwidth
    width = '%.0fHz'%(width_freq_Hz)
    start_freq_Hz = linefreq - 0.5*(line_velocity_width/2.99792458e5)*ref_freq_Hz #<TODO># lowest freq (as document says) or left-most freq (depending on positive/negative chanwidth)?
    start = '%.0fHz'%(start_freq_Hz)
    restfreq = '%.0fHz'%(ref_freq_Hz)
    nchan = (line_velocity_width/2.99792458e5)*ref_freq_Hz / width_freq_Hz # the output number of channels, covering the full line_velocity_width
    nchan = int(np.round(nchan))
    
    # 
    # select galaxy name
    # 
    matched_field_name, matched_field_indices, matched_field_phasecenters = get_field_phasecenters(dataset_ms, galaxy_name)
    
    # 
    # mstransform
    # 
    mstransform_parameters = {}
    mstransform_parameters['vis'] = dataset_ms
    mstransform_parameters['outputvis'] = output_ms
    mstransform_parameters['field'] = ','.join(matched_field_indices.astype(str))
    mstransform_parameters['spw'] = ','.join(linespws.astype(str))
    mstransform_parameters['datacolumn'] = datacolumn
    mstransform_parameters['regridms'] = True
    mstransform_parameters['mode'] = mode
    mstransform_parameters['start'] = start
    mstransform_parameters['width'] = width
    mstransform_parameters['nchan'] = nchan
    mstransform_parameters['nspw'] = 1
    #mstransform_parameters['combinespws'] = True
    mstransform_parameters['outframe'] = 'LSRK'
    mstransform_parameters['veltype'] = 'radio'
    mstransform_parameters['timeaverage'] = True
    mstransform_parameters['timebin'] = '60s'
    # 
    # Reset tclean parameters
    # 
    #default(mstransform)
    # 
    # Apply to global variables
    # 
    #for t in mstransform_parameters:
    #    globals()[t] = mstransform_parameters[t]
    #    locals()[t] = mstransform_parameters[t]
    
    #inp(mstransform)
    
    #saveinputs('mstransform', os.path.dirname(output_ms)+os.sep+'saved_mstransform_inputs.txt')
    
    #mstransform()
    
    with open(os.path.dirname(output_ms)+os.sep+'saved_mstransform_inputs.json', 'w') as fp:
        json.dump(mstransform_parameters, fp, indent = 4)
    
    mstransform(**mstransform_parameters)
    
    if not os.path.isdir(output_ms):
        raise Exception('Error! Failed to run mstransform and produce "%s"!'%(output_ms))





# 
# def prepare_clean_parameters
# 
def prepare_clean_parameters(vis, imagename, imcell = None, imsize = None, niter = 30000, calcres = True, calcpsf = True, 
                             phasecenter = '', field = '', pbmask = 0.2, pblimit = 0.1, threshold = 0.0, specmode = 'cube'):
    # 
    # Requires CASA module/function tb.
    # 
    casalog.origin('prepare_clean_parameters')
    
    # 
    # check field, makes sure there is only one field -- not exactly true, because for mosaics there are multiple fields but they are the same galaxy.
    # 
    #tb.open(vis+os.sep+'FIELD')
    #field_count = tb.nrows()
    #tb.close()
    #if field_count > 1:
    #    raise Exception('Error! The input vis has multiple fields! Please split the target field before calling prepare_clean_parameters()!')
    
    # 
    # check field, makes sure there is only one spw
    # 
    tb.open(vis+os.sep+'SPECTRAL_WINDOW')
    spw_count = tb.nrows()
    spw_ref_freq_col = tb.getcol('REF_FREQUENCY')
    tb.close()
    # 
    if spw_count > 1:
        check_spw_ref_freq_consistency = True
        for ispw in range(len(spw_ref_freq_col)):
            if spw_ref_freq_col[ispw] != spw_ref_freq_col[0]:
                check_spw_ref_freq_consistency = False
                break
        if check_spw_ref_freq_consistency == False:
            raise Exception('Error! The input vis "%s" has multiple spws and they do not have the same REF_FREQUENCY! Please split the target line channels into one spw before calling prepare_clean_parameters()!'%(vis))
    # 
    ref_freq_Hz = spw_ref_freq_col[0]
    
    # 
    # check imcell
    # 
    if imcell is None:
        # 
        # get state and obs_mode
        # 
        tb.open(vis+os.sep+'STATE')
        state_ID_col = np.arange(tb.nrows())
        state_obs_mode_col = tb.getcol('OBS_MODE')
        tb.close()
        # 
        valid_state_indices = []
        for i,obs_mode in enumerate(state_obs_mode_col):
            if obs_mode.find('ON_SOURCE') >= 0 or obs_mode.find('OBSERVE_TARGET') >= 0:
                valid_state_indices.append(i)
        # 
        if 0 == len(valid_state_indices):
            raise ValueError('Error! No valid state "ON_SOURCE" in the "STATE" of the input vis "%s"! The "STATE" table contains following "OBS_MODE": %s'%(vis, str(state_obs_mode_col)))
        # 
        valid_state_indices = np.array(valid_state_indices)
        
        # 
        # get data_desc_id to spw_id mapper
        # 
        #tb.open(vis+os.sep+'DATA_DESCRIPTION')
        #data_desc_spw_col = tb.getcol('SPECTRAL_WINDOW_ID')
        #tb.close()
        #
        #valid_data_desc_indicies = []
        #for i,data_desc_spw in enumerate(data_desc_spw_col):
        #    if data_desc_spw in linespws:
        #        valid_data_desc_indicies.append(i)
        #
        #if 0 == len(valid_data_desc_indicies):
        #    raise ValueError('Error! No valid data desc spw in the "DATA_DESCRIPTION" of the input vis "%s"!'%(vis))
        #
        #valid_data_desc_indicies = np.array(valid_data_desc_indicies)
        
        # 
        # get uvdist
        # 
        tb.open(vis)
        uvw = tb.getcol('UVW') # shape (3, nrows), each row is a (u,v,w) array
        ##print2('tb.query(\'FIELD_ID in [%s] AND DATA_DESC_ID in [%s] AND STATE_ID in [%s]\', \'UVW\')'%(','.join(matched_field_indices.astype(str)), ','.join(valid_data_desc_indicies.astype(str)), ','.join(valid_state_indices.astype(str))))
        ##result = tb.query('FIELD_ID in [%s] AND DATA_DESC_ID in [%s] AND STATE_ID in [%s]'%(','.join(matched_field_indices.astype(str)), ','.join(valid_data_desc_indicies.astype(str)), ','.join(valid_state_indices.astype(str))), 'UVW')
        #print2('tb.query(\'FIELD_ID in [%s] AND STATE_ID in [%s]\', \'UVW\')'%(','.join(matched_field_indices.astype(str)), ','.join(valid_state_indices.astype(str))))
        #result = tb.query('FIELD_ID in [%s] AND STATE_ID in [%s]'%(','.join(matched_field_indices.astype(str)), ','.join(valid_state_indices.astype(str))), 'UVW')
        #uvw = result.getcol('UVW')
        tb.close()
        # 
        uvdist = np.sqrt(np.sum(np.square(uvw[0:2, :]), axis=0))
        maxuvdist = np.max(uvdist)
        print2('maxuvdist = %s [m]'%(maxuvdist))
        L80uvdist = np.percentile(uvdist, 80) # np.max(uvdist) # now I am using 90-th percentile of baselies, same as used by 'analysisUtils.py' pickCellSize() getBaselineStats(..., percentile=...)
        print2('L80uvdist = %s [m] (80-th percentile)'%(L80uvdist))
        # 
        synbeam = 2.99792458e8 / ref_freq_Hz / maxuvdist / np.pi * 180.0 * 3600.0 # arcsec
        synbeam = 0.574 * 2.99792458e8 / ref_freq_Hz / L80uvdist / np.pi * 180.0 * 3600.0 # arcsec # .574lambda/L80, see 'analysisUtils.py' estimateSynthesizedBeamFromASDM()
        synbeam_nprec = 2 # keep 2 valid 
        synbeam_ndigits = (synbeam_nprec-1) - int(np.floor(np.log10(synbeam))) # keep to these digits (precision) after point, e.g., 0.1523 -> nprec 2 -> round(0.1523*100)/100 = 0.15
        synbeam = (np.round(synbeam * 10**(synbeam_ndigits))) / 10**(synbeam_ndigits)
        oversampling = 5.0
        imcell_arcsec = synbeam / oversampling
        imcell = '%sarcsec'%(imcell_arcsec)
        print2('synbeam = %s [arcsec]'%(synbeam))
        print2('imcell_arcsec = %s'%(imcell_arcsec))
        # 
    else:
        imcell_arcsec = float(re.sub(r'arcsec$', r'', imcell))
    # 
    # check imsize
    # 
    if imsize is None:
        # 
        # get antdiam
        # 
        minantdiam = get_antenn_diameter(vis)
        # 
        # get field and phase centers
        # 
        matched_field_name, matched_field_indices, matched_field_phasecenters = get_field_phasecenters(vis, field)
        matched_field_min_RA_deg = np.min(matched_field_phasecenters[0, :])
        matched_field_max_RA_deg = np.max(matched_field_phasecenters[0, :])
        matched_field_min_Dec_deg = np.min(matched_field_phasecenters[1, :])
        matched_field_max_Dec_deg = np.max(matched_field_phasecenters[1, :])
        #print('matched_field_phasecenters.shape', matched_field_phasecenters.shape)
        #print('matched_field_phasecenters:', matched_field_phasecenters)
        #raise NotImplementedError()
        # 
        # calc primary beam
        # 
        pribeam = 1.13  * (2.99792458e8 / ref_freq_Hz / minantdiam / np.pi * 180.0 ) # in units of degrees
        print2('pribeam = %s [arcsec]'%(pribeam * 3600.0))
        print2('matched_field_min_RA_deg = %s'%(matched_field_min_RA_deg))
        print2('matched_field_max_RA_deg = %s'%(matched_field_max_RA_deg))
        print2('matched_field_min_Dec_deg = %s'%(matched_field_min_Dec_deg))
        print2('matched_field_max_Dec_deg = %s'%(matched_field_max_Dec_deg))
        imsize_RA_deg = (matched_field_max_RA_deg - matched_field_min_RA_deg) * np.cos(np.deg2rad((matched_field_max_Dec_deg+matched_field_min_Dec_deg)/2.0)) + 2.0*pribeam
        imsize_Dec_deg = (matched_field_max_Dec_deg - matched_field_min_Dec_deg) + 2.0*pribeam
        imsize_RA = imsize_RA_deg / (imcell_arcsec / 3600.0) # pixels
        imsize_Dec = imsize_Dec_deg / (imcell_arcsec / 3600.0) # pixels
        print2('imsize_RA = %s [arcsec]'%(imsize_RA_deg * 3600.0))
        print2('imsize_Dec = %s [arcsec]'%(imsize_Dec_deg * 3600.0))
        imsize = [get_optimized_imsize(imsize_RA), get_optimized_imsize(imsize_Dec)]
    # 
    print2('imsize = %s'%(imsize))
    # 
    # We can also use analysisUtils, but the results are very similar to my above implementation.
    # 
    #au_cellsize, au_imsize, au_centralField = au.pickCellSize(vis, imsize=True, npix=5)
    #au.plotmosaic(vis, 'NGC_4321')
    # 
    # Prepare tclean parameters
    # 
    clean_parameters = {}
    clean_parameters['vis'] = vis
    clean_parameters['selectdata'] = True
    clean_parameters['field'] = '' # ','.join(matched_field_indices.astype(str))
    clean_parameters['spw'] = ''
    clean_parameters['phasecenter'] = phasecenter
    clean_parameters['cell'] = imcell # tclean parameter name
    clean_parameters['imsize'] = imsize # tclean parameter name
    clean_parameters['imagename'] = imagename # output_dir+os.sep+'%s_%s_cleaned'%(galaxy_name_cleaned, linename)
    clean_parameters['gridder'] = 'mosaic' # 'standard'
    clean_parameters['specmode'] = specmode # 'cube' # for spectral line cube
    clean_parameters['outframe'] = 'LSRK'
    clean_parameters['deconvolver'] = 'hogbom'
    clean_parameters['usemask'] = 'pb' # construct a 1/0 mask at the 0.2 level
    clean_parameters['pbmask'] = pbmask # data outside this pbmask will not be fitted
    clean_parameters['threshold'] = threshold
    #clean_parameters['usemask'] = 'user' #<TODO># 
    #clean_parameters['mask'] = '' #<TODO># 
    clean_parameters['pblimit'] = pblimit # data outside this pblimit will be output as NaN
    clean_parameters['pbcor'] = True # create both pbcorrected and uncorrected images
    clean_parameters['restoration'] = True
    clean_parameters['restoringbeam'] = '%sarcsec'%(synbeam) # Automatically estimate a common beam shape/size appropriate for all planes.
    clean_parameters['nterms'] = 1 # nterms must be ==1 when deconvolver='hogbom' is chosen
    clean_parameters['chanchunks'] = -1 # This feature is experimental and may have restrictions on how chanchunks is to be chosen. For now, please pick chanchunks so that nchan/chanchunks is an integer. 
    clean_parameters['interactive'] = False
    clean_parameters['savemodel'] = 'virtual' # 'none', 'virtual', 'modelcolumn'. 'virtual' for simple gridding, 'modelcolumn' for gridder='awproject'.
    #niter = 30000
    #calcres = True # calculate initial residual image at the beginning of the first major cycle
    #calcpsf = True
    if calcres==False and niter==0:
        print2('Note: Only the PSF will be made and no data will be gridded in the first major cycle of cleaning.')
    elif calcres==False and niter>0:
        print2('Note: We will assume that a "*.residual" image already exists and that the minor cycle can begin without recomputing it.')
    elif calcres==True:
        if calcpsf==False and not (os.path.isfile(imagename+'.psf') and os.path.isfile(imagename+'.sumwt')):
            calcpsf = True # calcres=True requires that calcpsf=True or that the .psf and .sumwt images already exist on disk (for normalization purposes)
    clean_parameters['niter'] = niter
    clean_parameters['calcres'] = calcres
    clean_parameters['calcpsf'] = calcpsf
    
    # 
    # Check mpicasa
    if 'processor_origin' in globals():
        global processor_origin
        if processor_origin.find('MPIClient') >= 0:
            print2('going parallel!')
            clean_parameters['parallel'] = True
    
    # 
    # Reset tclean parameters
    # 
    #default(tclean)
    
    # 
    # Apply to global variables
    # 
    #for t in clean_parameters:
    #    globals()[t] = clean_parameters[t]
    
    # 
    # Return
    # 
    return clean_parameters




def run_tclean_with_clean_parameters(clean_parameters):
    # 
    # Requires CASA module/function default, inp, saveinputs, tclean.
    # 
    casalog.origin('run_tclean_with_clean_parameters')
    # 
    # Reset tclean parameters
    # 
    #default(tclean)
    # 
    # Apply to global variables
    # 
    #for t in clean_parameters:
    #    globals()[t] = clean_parameters[t]
    #    locals()[t] = clean_parameters[t]
    # 
    # Load vis, imagename
    # 
    vis = clean_parameters['vis']
    imagename = clean_parameters['imagename']
    # 
    # Check existing file
    # 
    if os.path.isdir(imagename+'.image'):
        raise Exception('Error! Found existing "%s"! Will not overwrite!'%(imagename+'.image'))
        return
    # 
    # Print and save inputs
    # 
    #inp(tclean)
    #saveinputs('tclean', os.path.dirname(imagename)+os.sep+'saved_tclean_inputs.txt')
    # 
    # Run tclean
    # 
    #tclean()
    # 
    # 
    # Print and save inputs (New method)
    # 
    with open(os.path.dirname(imagename)+os.sep+'saved_tclean_inputs.json', 'w') as fp:
        json.dump(clean_parameters, fp, indent = 4)
    # 
    # Run tcleans (New method)
    # 
    tclean(**clean_parameters)
    # 
    # Check outputs
    # 
    if os.path.isdir(imagename+'.image'):
        print2('Cleaning seems finished sucessfully.')
        exportfits(imagename+'.image', imagename+'.image.fits')
        exportfits(imagename+'.image.pbcor', imagename+'.image.pbcor.fits')
        exportfits(imagename+'.psf', imagename+'.psf.fits')
        exportfits(imagename+'.pb', imagename+'.pb.fits')
        exportfits(imagename+'.model', imagename+'.model.fits')
        exportfits(imagename+'.residual', imagename+'.residual.fits')
    else:
        raise Exception('Error! tclean failed to produce the output image "%s"!'%(imagename+'.image'))




def make_dirty_image(vis, imagename, **kwargs):
    # 
    # Check existing file
    # 
    if os.path.isdir(imagename+'.image'):
        print2('Found existing "%s"! Will not overwrite it! Skipping make_dirty_image()!'%(imagename+'.image'))
        return
    else:
        # 
        # prepare_clean_parameters
        clean_parameters = prepare_clean_parameters(vis, imagename, niter = 0, **kwargs)
        # 
        # run_tclean_with_clean_parameters
        run_tclean_with_clean_parameters(clean_parameters)
        # 
        # copy saved_tclean_inputs.txt
        #shutil.copy(os.path.dirname(imagename)+os.sep+'saved_tclean_inputs.txt', os.path.dirname(imagename)+os.sep+'saved_tclean_inputs_for_dirty_cube.txt')
        shutil.copy(os.path.dirname(imagename)+os.sep+'saved_tclean_inputs.json', os.path.dirname(imagename)+os.sep+'saved_tclean_inputs_for_dirty_cube.json')




def make_clean_image(vis, imagename, **kwargs):
    # 
    # Check existing file
    # 
    if os.path.isdir(imagename+'.image'):
        print2('Found existing "%s"! Will not overwrite it! Skipping make_clean_image()!'%(imagename+'.image'))
        return
    else:
        # 
        # prepare_clean_parameters
        clean_parameters = prepare_clean_parameters(vis, imagename, niter = 30000, **kwargs)
        # 
        # run_tclean_with_clean_parameters
        run_tclean_with_clean_parameters(clean_parameters)
        # 
        # copy saved_tclean_inputs.txt
        #shutil.copy(os.path.dirname(imagename)+os.sep+'saved_tclean_inputs.txt', os.path.dirname(imagename)+os.sep+'saved_tclean_inputs_for_clean_cube.txt')
        shutil.copy(os.path.dirname(imagename)+os.sep+'saved_tclean_inputs.json', os.path.dirname(imagename)+os.sep+'saved_tclean_inputs_for_clean_cube.json')




def make_clean_image_of_continuum(vis, imagename, **kwargs):
    # 
    # Check existing file
    # 
    if os.path.isdir(imagename+'.image'):
        print2('Found existing "%s"! Will not overwrite it! Skipping make_clean_image()!'%(imagename+'.image'))
        return
    else:
        # 
        # prepare_clean_parameters
        clean_parameters = prepare_clean_parameters(vis, imagename, niter = 30000, specmode = 'mfs', **kwargs)
        # 
        # run_tclean_with_clean_parameters
        run_tclean_with_clean_parameters(clean_parameters)
        # 
        # copy saved_tclean_inputs.txt
        #shutil.copy(os.path.dirname(imagename)+os.sep+'saved_tclean_inputs.txt', os.path.dirname(imagename)+os.sep+'saved_tclean_inputs_for_clean_cube_of_continuum.txt')
        shutil.copy(os.path.dirname(imagename)+os.sep+'saved_tclean_inputs.json', os.path.dirname(imagename)+os.sep+'saved_tclean_inputs_for_clean_cube_of_continuum.json')




# 
# def project_fits_cube
# 
def project_fits_cube(input_fits_cube, template_fits_cube, output_fits_cube, overwrite = False):
    # 
    # We will project the input_fits_cube to the pixel grid of the template_fits_cube_or_wcs
    # by matching the World Coordinate System (WCS). 
    # 
    # We can process both fits cubes and images. 
    # 
    # No CASA module required. 
    # 
    # Superseding part of the process_clean_mask() code
    # 
    import warnings
    from astropy.utils.exceptions import AstropyWarning
    warnings.simplefilter('ignore', category=AstropyWarning)
    from astropy.table import Table
    from astropy.io import fits
    from astropy import units as u
    from astropy import wcs
    from astropy.wcs import WCS
    from astropy.wcs.utils import proj_plane_pixel_scales
    from astropy.coordinates import SkyCoord, FK5
    from scipy.interpolate import griddata
    # 
    # Read template_fits_cube_or_wcs
    if type(template_fits_cube) is astropy.wcs.WCS:
        # If input is a fits wcs
        template_fits_header = template_fits_cube_or_wcs.to_header()
    elif type(template_fits_cube) is astropy.io.fits.Header:
        # If input is a fits wcs
        template_fits_header = template_fits_cube
    else:
        # If input is a fits file
        print('Reading "%s"'%(template_fits_cube))
        template_hdulist = fits.open(template_fits_cube)
        template_hdu = template_hdulist[0]
        template_fits_header = template_hdu.header
    # 
    # Read input_fits_cube
    print('Reading "%s"'%(input_fits_cube))
    input_hdulist = fits.open(input_fits_cube)
    input_hdu = input_hdulist[0]
    input_fits_header = input_hdu.header
    input_fits_data = input_hdu.data
    input_fits_data_shape = list(input_fits_data.shape)
    # 
    # Make sure the input is a 3D cube or at least a 2D image
    if int(input_fits_header['NAXIS']) < 2:
        raise Exception('Error! The input fits cube "%s" does not have more than 2 dimensions!'%(input_fits_cube))
    if int(template_fits_header['NAXIS']) < 2:
        raise Exception('Error! The input fits cube "%s" does not have more than 2 dimensions!'%(template_fits_cube))
    # 
    # Make sure the input fits data have at least equal dimensions as the template fits data. For extra higher dimensions, we will loop them over.
    #if int(input_fits_header['NAXIS']) < int(template_fits_header['NAXIS']):
    #    raise Exception('Error! The input fits cube "%s" has %d dimensions but the template "%s" has %d dimensions!'%(input_fits_cube, input_fits_header['NAXIS'], template_fits_cube, template_fits_header['NAXIS']))
    # 
    # Take the template fits dimension.
    naxis = min(int(input_fits_header['NAXIS']), int(template_fits_header['NAXIS']))
    # 
    # Do velocity-to-frequency conversion before checking header consistency
    if naxis >= 3:
        if input_fits_header['CTYPE3'].strip().upper() == 'VRAD' and template_fits_header['CTYPE3'].strip().upper() == 'FREQ':
            ctype3 = input_fits_header['CTYPE3']
            cunit3 = input_fits_header['CUNIT3'].strip().replace(' ','').lower()
            crpix3 = input_fits_header['CRPIX3']
            crval3 = input_fits_header['CRVAL3']
            cdelt3 = input_fits_header['CDELT3']
            if cunit3 == 'km/s' or cunit3 == 'kms-1':
                c30 = 2.99792458e5
            else:
                c30 = 2.99792458e8
            input_fits_header['CRVAL3'] = (1.0-(crval3/c30))*input_fits_header['RESTFRQ'] # Hz, (nu0-nu)/nu0 = (v/c), so nu = (1-(v/c))*nu0
            input_fits_header['CDELT3'] = (-(cdelt3/c30))*input_fits_header['RESTFRQ'] # Hz, reversed order
            input_fits_header['CTYPE3'] = 'FREQ'
            input_fits_header['CUNIT3'] = 'Hz'
    # 
    # Check header consistency and record NAXISi
    inaxis = []
    tnaxis = []
    for i in range(1, naxis+1):
        if input_fits_header['CTYPE%d'%(i)].strip().upper() != template_fits_header['CTYPE%d'%(i)].strip().upper():
            raise Exception('Error! The input fits cube "%s" CTYPE%d is %s but the template "%s" CTYPE%d is %s!'%(\
                                input_fits_cube, i, input_fits_header['CTYPE%d'%(i)].strip(), \
                                template_fits_cube, i, template_fits_header['CTYPE%d'%(i)].strip() ) )
        inaxis.append(int(input_fits_header['NAXIS%d'%(i)]))
        tnaxis.append(int(template_fits_header['NAXIS%d'%(i)]))
    inaxis = np.array(inaxis[::-1]) # [nx, ny, nchan, ....], it is inverted to the Python array dimension order
    tnaxis = np.array(tnaxis[::-1]) # [nx, ny, nchan, ....], it is inverted to the Python array dimension order
    inaxis_str = 'x'.join(inaxis.astype(str))
    tnaxis_str = 'x'.join(tnaxis.astype(str))
    # 
    # If input fits data have extra higher dimensions than the template data, we will condense the additional dimensions into one extra dimension, and later we will loop them over and do the interpolation slice-by-slice.
    idataslice = 0
    if int(input_fits_header['NAXIS']) > naxis:
        print2('Warning! The input fits cube "%s" has %d dimensions while the template "%s" has only %d dimensions! We will loop the extra higher dimensions of the input data slice-by-slice and project to template pixel grid.'%(input_fits_cube, input_fits_header['NAXIS'], template_fits_cube, template_fits_header['NAXIS']))
        idatashape = copy.copy(input_fits_data_shape[::-1][0:naxis])
        idataslice = np.product(input_fits_data_shape[::-1][naxis:])
        idatashape.append(idataslice)
        idatashape = idatashape[::-1]
        input_fits_data.shape = idatashape
        output_fits_data_shape = copy.copy(tnaxis[::-1])
        # here we reshape the input data to shrink all higher dimensions into one extra dimension, 
        # e.g., from (m,n,a,b,c) to (x,a,b,c), when the template data shape is (a,b,c), and x = m*n. 
    # 
    # Otherwise if the template fits data have extra higher dimensions, we will only use the common dimensions, while pad extra dimensions with np.newaxis
    elif int(template_fits_header['NAXIS']) > naxis:
        print2('Warning! The input fits cube "%s" has %d dimensions while the template "%s" has %d dimensions! We will ignore the extra higher dimensions in the template for computation and simply reshape the output data to the template dimensions.'%(input_fits_cube, input_fits_header['NAXIS'], template_fits_cube, template_fits_header['NAXIS']))
        odatashape = []
        idataslice = -(int(template_fits_header['NAXIS']) - naxis)
        for i in range(naxis):
            odatashape.append(tnaxis[i]) # template_fits_header['NAXIS%d'(i+1)]
        for i in range(naxis, naxis-idataslice):
            odatashape.append(1) # template_fits_header['NAXIS%d'(i+1)]
        odatashape = odatashape[::-1]
        output_fits_data_shape = odatashape
    # 
    # Otherwise the output fits data shape is the same as the template fits data shape
    else:
        output_fits_data_shape = copy.copy(tnaxis[::-1])
    # 
    # Get input fits WCS
    iwcs = WCS(input_fits_header, naxis=naxis)
    idata = input_fits_data
    # 
    # Get template fits WCS
    twcs = WCS(template_fits_header, naxis=naxis)
    tdatashape = tnaxis[::-1] # this will also be the shape of (each slice of) the output data array
    # 
    # Make pixel mgrid
    timestart = timeit.default_timer()
    if naxis == 2:
        # 
        print2('Generating pixel mgrid with %s pixels'%(inaxis_str))
        iy, ix = np.mgrid[0:inaxis[1], 0:inaxis[0]]
        ipixcoords = np.column_stack([ix.flatten(), iy.flatten()])
        # 
        print2('Generating pixel mgrid with %s pixels'%(tnaxis_str))
        ty, tx = np.mgrid[0:tnaxis[1], 0:tnaxis[0]]
        tpixcoords = np.column_stack([tx.flatten(), ty.flatten()])
        # 
    elif naxis == 3:
        # 
        print2('Generating pixel mgrid with %s pixels'%(inaxis_str))
        ichan, iy, ix = np.mgrid[0:inaxis[2], 0:inaxis[1], 0:inaxis[0]]
        ipixcoords = np.column_stack([ix.flatten(), iy.flatten(), ichan.flatten()])
        # 
        print2('Generating pixel mgrid with %s pixels'%(tnaxis_str))
        tchan, ty, tx = np.mgrid[0:tnaxis[2], 0:tnaxis[1], 0:tnaxis[0]]
        tpixcoords = np.column_stack([tx.flatten(), ty.flatten(), tchan.flatten()])
        # 
    else:
        raise NotImplementedError('Error! The cube projection and interpolation have not been implemented for an NAXIS of %d!'%(naxis))
    # 
    timestop = timeit.default_timer()
    print('Used %s seconds'%(timestop-timestart))
    # 
    #raise NotImplementedError() # debug point
    # 
    timestart = timeit.default_timer()
    # 
    # Convert each pixel coordinate to skycoordinate for the template pixel grid which is also the output pixel grid.
    print2('Computing wcs_pix2world for %s pixels'%(tnaxis_str))
    oskycoords = twcs.wcs_pix2world(tpixcoords, 0)
    print2('oskycoords.shape = %s'%(list(oskycoords.shape)))
    #tra, tdec, tfreq = oskycoords.T
    # 
    # Convert each pixel skycoordinate to the coordinate in the input mask cube, so that we can do interpolation. 
    print2('Computing wcs_world2pix for %s pixels'%(tnaxis_str))
    opixcoords = iwcs.wcs_world2pix(oskycoords, 0)
    print2('opixcoords.shape = %s'%(list(opixcoords.shape)))
    # 
    timestop = timeit.default_timer()
    print('Used %s seconds'%(timestop-timestart))
    # 
    # Loop each input data slice
    odata = []
    odataarray = None
    for i in range(max(1,idataslice)):
        # 
        # Do interpolation with scipy.interpolate.griddata
        timestart = timeit.default_timer()
        print2('Interpolating griddata...')
        if idataslice > 0:
            idataarray = idata[i].flatten()
        else:
            idataarray = idata.flatten()
        idatamask = ~np.isnan(idataarray)
        odataarray = griddata(ipixcoords[idatamask], \
                              idataarray[idatamask], \
                              opixcoords, \
                              method = 'nearest', \
                              fill_value = 0 )
        timestop = timeit.default_timer()
        print('Used %s seconds'%(timestop-timestart))
        # 
        # The interpolation is done with serialized arrays, so we reshape the output interpolated aray to 3D cube 
        odataarray = odataarray.reshape(tdata.shape).astype(idata.dtype)
        if idataslice > 0:
            odata.append(odataarray)
        else:
            odata = odataarray
    #odatashape = list(odataarray.shape)[::-1]
    #odatashape.append(input_fits_data_shape[::-1][naxis:])
    #odatashape = odatashape[::-1]
    output_fits_data_shape = output_fits_data_shape[::-1]
    output_fits_data_shape[0:naxis] = list(odata.shape)[::-1][0:naxis]
    output_fits_data_shape = output_fits_data_shape[::-1]
    print('input_fits_data_shape = %s'%(input_fits_data_shape))
    print('output_fits_data_shape = %s'%(output_fits_data_shape))
    print('odata.shape = %s'%(list(odata.shape)))
    odata.shape = output_fits_data_shape
    # 
    # Prepare output fits header
    output_fits_header = twcs.to_header()
    if idataslice > 0:
        output_fits_header['NAXIS'] = len(input_fits_data_shape)
        for i in range(naxis, len(input_fits_data_shape)):
            for keybase in ['NAXIS', 'CTYPE', 'CUNIT', 'CRPIX', 'CRVAL', 'CDELT']:
                key = keybase+'%d'%(i+1)
                if key in input_fits_header:
                    output_fits_header[key] = input_fits_header[key]
    elif idataslice < 0:
        output_fits_header['NAXIS'] = len(output_fits_data_shape)
        for i in range(naxis, len(output_fits_data_shape)):
            for keybase in ['NAXIS', 'CTYPE', 'CUNIT', 'CRPIX', 'CRVAL', 'CDELT']:
                key = keybase+'%d'%(i+1)
                if key in template_fits_header:
                    output_fits_header[key] = template_fits_header[key]
    # 
    # Output the interpolated mask cube as a fits file
    print2('Writing fits file...')
    output_hdu = fits.PrimaryHDU(data = odata, header = output_fits_header)
    if not re.match(r'.*\.fits$', output_fits_cube, re.IGNORECASE):
        output_fits_cube += '.fits'
    output_hdu.writeto(output_mask_cube, overwrite = overwrite)
    print2('Output to "%s"!'%(output_fits_cube))

# 
# def test_project_fits_cube
# 
def test_project_fits_cube():
    project_fits_cube(input_fits_cube = 'ngc4321_co21_clean_mask.fits', 
                      template_fits_cube = 'run_tclean_2015.1.00956.S._.12m._.1/ngc4321_co21_dirty.image.pbcor.fits', 
                      output_fits_cube = 'test_project_fits_cube.fits', 
                      overwrite = True)

# 
# test here!
# 
#test_project_fits_cube()
#raise NotImplementedError()




# 
# def process_clean_mask
# 
def process_clean_mask(input_mask_cube, template_image_cube, output_mask_cube):
    # 
    # We will read the input_mask_cube, register the cube to the template_image_cube WCS, 
    # and output the mask cube. 
    # 
    # No CASA module required
    # 
    # TODO: The user can supply some parameters to add some mask conditions using the template_image_cube
    # 
    import warnings
    from astropy.utils.exceptions import AstropyWarning
    warnings.simplefilter('ignore', category=AstropyWarning)
    from astropy.table import Table
    from astropy.io import fits
    from astropy import units as u
    from astropy import wcs
    from astropy.wcs import WCS
    from astropy.wcs.utils import proj_plane_pixel_scales
    from astropy.coordinates import SkyCoord, FK5
    from scipy.interpolate import griddata
    # 
    # Read template_image_cube
    template_image_cube_file = template_image_cube
    for tempfile in [template_image_cube+'.image.pbcor.fits', template_image_cube+'.image.fits', template_image_cube+'.fits']:
        if os.path.isfile(tempfile):
            template_image_cube_file = tempfile
            break
    template_hdulist = fits.open(template_image_cube_file)
    template_hdu = template_hdulist[0]
    # 
    # Get WCS
    tnstokes, tnchan, tny, tnx = template_hdu.data.shape
    twcs = WCS(template_hdu.header, naxis=3)
    tdata = template_hdu.data
    # 
    # Read input_mask_cube
    input_hdulist = fits.open(input_mask_cube)
    input_hdu = input_hdulist[0]
    if input_hdu.header['CTYPE3'].strip() == 'VRAD':
        ctype3 = input_hdu.header['CTYPE3']
        cunit3 = input_hdu.header['CUNIT3'].strip().replace(' ','').lower()
        crpix3 = input_hdu.header['CRPIX3']
        crval3 = input_hdu.header['CRVAL3']
        cdelt3 = input_hdu.header['CDELT3']
        if cunit3 == 'km/s' or cunit3 == 'kms-1':
            c30 = 2.99792458e5
        else:
            c30 = 2.99792458e8
        input_hdu.header['CRVAL3'] = (1.0-(crval3/c30))*input_hdu.header['RESTFRQ'] # Hz, (nu0-nu)/nu0 = (v/c), so nu = (1-(v/c))*nu0
        input_hdu.header['CDELT3'] = (-(cdelt3/c30))*input_hdu.header['RESTFRQ'] # Hz, reversed order
        input_hdu.header['CTYPE3'] = 'FREQ'
        input_hdu.header['CUNIT3'] = 'Hz'
    inchan, iny, inx = input_hdu.data.shape
    iwcs = WCS(input_hdu.header, naxis=3)
    idata = input_hdu.data
    # 
    # Prepare pixel mgrid
    print2('generating pixel mgrid with %dx%dx%d pixels'%(inx, iny, inchan))
    igchan, igy, igx = np.mgrid[0:inchan, 0:iny, 0:inx]
    # 
    print2('generating pixel mgrid with %dx%dx%d pixels'%(tnx, tny, tnchan))
    tgchan, tgy, tgx = np.mgrid[0:tnchan, 0:tny, 0:tnx]
    # 
    #raise NotImplementedError() # debug point
    # 
    # Convert each pixel coordinate to skycoordinate for the template pixel grid which is also the output pixel grid.
    print2('computing wcs_pix2world for %dx%dx%d pixels'%(tnx, tny, tnchan))
    oskycoords = twcs.wcs_pix2world(np.column_stack([tgx.flatten(), tgy.flatten(), tgchan.flatten()]), 0)
    print2('oskycoords.shape = %s'%(list(oskycoords.shape)))
    #tra, tdec, tfreq = oskycoords.T
    # 
    # Convert each pixel skycoordinate to the coordinate in the input mask cube, so that we can do interpolation. 
    print2('computing wcs_world2pix for %dx%dx%d pixels'%(tnx, tny, tnchan))
    opixcoords = iwcs.wcs_world2pix(oskycoords, 0)
    print2('opixcoords.shape = %s'%(list(opixcoords.shape)))
    ogx, ogy, ogchan = opixcoords.T
    # 
    # Do interpolation with scipy.interpolate.griddata
    #imask = ~np.isnan(idata) <TODO> Is NaN a problem?
    print2('interpolating griddata...')
    odata = griddata(np.column_stack([igchan.flatten(), igy.flatten(), igx.flatten()]), idata.flatten(), \
                     np.column_stack([ogchan.flatten(), ogy.flatten(), ogx.flatten()]), \
                     method='nearest', fill_value=0)
    # 
    # The interpolation is done with serialized arrays, so we reshape the output interpolated aray to 3D cube 
    odata = odata.reshape(tdata.shape).astype(int)
    # 
    # <TODO> To implement some mask operations according to the input template_image_cube.
    #        for example, do an AND operation to combine current mask and all pixels with S/N>4 in the template image cube. 
    # <TODO> 
    # 
    # Output the interpolated mask cube as a fits file
    print2('writing fits file...')
    output_hdu = fits.PrimaryHDU(data = odata, header = twcs.to_header())
    if not re.match(r'.*\.fits$', output_mask_cube, re.IGNORECASE):
        output_mask_cube += '.fits'
    output_hdu.writeto(output_mask_cube, overwrite = True)
    print2('Output to "%s"!'%(output_mask_cube))















def dzliu_clean():
    # 
    # Reset tclean and prepare user input paraameters
    # 
    dataset_ms = 'Level_3_Concat/pointing6.ms'
    dataset_name = re.sub(r'\.ms$', r'', os.path.basename(dataset_ms))
    galaxy_name = 'pointing6' # 
    galaxy_center = ''
    phasecenter = '' # galaxy_center
    galaxy_name_cleaned = re.sub(r'[^a-zA-Z0-9]', r'', galaxy_name).lower()
    #lab_line_name = 'HI21cm'
    #lab_line_freq = 1.420405751 # GHz
    linename = 'HI21cm'
    #linefreq = lab_line_freq*1e9 / (1.0+redshift) # Hz
    linevcen = 10230.0 # km/s, optical definition
    linevwid = 1000.0 # km/s
    linevres = -1 # 50.0 # km/s
    redshift = linevcen / 2.99792458e5
    
    # 
    # Prepare output dir and name
    # 
    output_dir = 'Level_4_Clean/run_tclean_%s'%(dataset_name)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    line_ms = output_dir+os.sep+'%s_%s.ms'%(galaxy_name_cleaned, linename)
    line_dirty_cube = output_dir+os.sep+'%s_%s_dirty'%(galaxy_name_cleaned, linename)
    line_clean_cube = output_dir+os.sep+'%s_%s_clean'%(galaxy_name_cleaned, linename)
    line_clean_cube_of_continuum = output_dir+os.sep+'%s_%s_clean'%(galaxy_name_cleaned, 'cont')
    
    # 
    # Split line data and make channel averaging
    # 
    split_line_visibilities(dataset_ms, line_ms, galaxy_name, linename, linevcen, linevwid, linevres)
    
    # 
    # Make dirty image
    # 
    make_dirty_image(line_ms, line_dirty_cube, phasecenter = phasecenter)
    
    #
    # Compute rms in the dirty image
    # 
    result_imstat_dict = imstat(line_dirty_cube+'.image')
    threshold = result_imstat_dict['rms'][0] * 3.0 #<TODO># 3-sigma
    
    # 
    # Make clean image
    # 
    make_clean_image(line_ms, line_clean_cube, phasecenter = phasecenter, threshold = threshold, pblimit = 0.05, pbmask = 0.05)
    
    # 
    # Make clean image of the rough continuum 
    # 
    make_clean_image_of_continuum(line_ms, line_clean_cube_of_continuum, phasecenter = phasecenter, threshold = threshold, pblimit = 0.05, pbmask = 0.05)






############
#   main   #
############

dzliu_main_func_name = 'dzliu_clean' # make sure this is the right main function in this script file

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




