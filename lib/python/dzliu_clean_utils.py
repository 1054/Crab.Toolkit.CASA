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
        #print2('matched_field_phasecenters = %s'%(matched_field_phasecenters))
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


def cleanup_tclean_products(imagename, cleanup_mask=True, cleanup_fits=True, exit_on_error=True):
    if imagename.endswith('.image'):
        imagename = re.sub(r'\.image$', r'', imagename)
    suffix_list = ['.image', '.image.pbcor', '.model', '.pb', '.psf', '.residual', '.sumwt', '.weight', '.tt0', '.alpha', '.beta'] #<TODO># depends on CASA version and tclean cube type
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


def apply_pbcor_to_tclean_image(imagename, overwrite=True, exit_on_error=True):
    if imagename.endswith('.image'):
        imagename = re.sub(r'\.image$', r'', imagename)
    infile = imagename+'.image'
    pbimage = imagename+'.pb'
    outfile = imagename+'.image.pbcor'
    if not os.path.isdir(infile) or not os.path.isdir(pbimage):
        raise Exception('Error! Data not found: "%s" or "%s"'%(infile, pbimage))
    if os.path.isdir(outfile):
        if not overwrite:
            raise Exception('Found existing data "%s"! Please clean it up first!'%(outfile))
        else:
            print2('Found existing data "%s", overwriting it.'%(outfile))
            shutil.rmtree(outfile)
    print2('Running CASA task: impbcor(imagename=%r, pbimage=%r, outfile=%r, mode=%r, cutoff=%s)'%(infile, pbimage, outfile, 'divide', 0.1))
    impbcor(imagename=infile, pbimage=pbimage, outfile=outfile, mode='divide', cutoff=0.1)
    if os.path.isdir(outfile):
        print2('Output to "%s"'%(outfile))
    else:
        if exit_on_error:
            raise Exception('Error! Failed to run CASA impbcor and output "%s"'%(outfile))
    # 


def export_tclean_products_as_fits_files(imagename, dropstokes=True, overwrite=True, exit_on_error=True):
    if imagename.endswith('.image'):
        imagename = re.sub(r'\.image$', r'', imagename)
    suffix_list = ['.image', '.image.pbcor', '.mask', '.model', '.pb', '.psf', '.residual', '.tt0', '.alpha', '.beta'] #<TODO># depends on CASA version and tclean cube type
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




