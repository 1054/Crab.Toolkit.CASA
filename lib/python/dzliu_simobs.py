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
import pyfits
import glob
from distutils.version import LooseVersion
from collections import OrderedDict
from taskinit import casalog, tb #, ms, iatool
#casalog.filter('DEBUG2')
# 
#import resource
#resource.setrlimit(resource.RLIMIT_AS, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
#resource.setrlimit(resource.RLIMIT_NOFILE, (10240, 40960)) # ulimit -Sn # ulimit -Hn
# 
#from taskinit import casac
#tb = casac.table
#from __casac__.table import table as tb
try:
    fn
except:
    from __casac__.functional import functional as fn
try:
    cs
except:
    from __casac__.coordsys import coordsys as cs
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
    from importfits_cli import importfits_cli_
    importfits = importfits_cli_()
    from concat_cli import concat_cli_
    concat = concat_cli_()
    from split_cli import split_cli_
    split = split_cli_()
    from imstat_cli import imstat_cli_
    imstat = imstat_cli_()
    from simobserve_cli import simobserve_cli_
    simobserve = simobserve_cli_()
else:
    # see CASA 6 updates here: https://alma-intweb.mtk.nao.ac.jp/~eaarc/UM2018/presentation/Nakazato.pdf
    from casatasks import tclean, mstransform, exportfits, importfits, concat, split, imstat
    #from casatasks import sdbaseline
    #from casatools import ia
# 
#if not os.path.isfile(os.getenv('CASAPATH').split(' ')[0]+'/data/alma/simmos/vla.a.cfg'):
#    raise Exception('Error! VLA config file not found from "$CASAPATH/data/alma/simmos/vla.a.cfg"!')
global casa_simmos_dir
casa_simmos_dir = os.getenv('CASAPATH').split(' ')[0]+'/data/alma/simmos'
if not os.path.isdir(casa_simmos_dir):
    raise Exception('Error! CASA simulation config file directory is not found: %r'%(casa_simmos_dir))



# 
# def print2
# 
def print2(message):
    print(message)
    casalog.post(message, 'INFO')




# 
# def arcsec2float
# 
def arcsec2float(arcsec_str):
    arcsec_value = np.nan
    if re.match(r'^.*arcsec$', arcsec_str.strip().lower()):
        arcsec_value = float(re.sub(r'arcsec$', r'', arcsec_str.strip().lower()))
    elif re.match(r'^.*arcmin$', arcsec_str.strip().lower()):
        arcsec_value = float(re.sub(r'arcmin$', r'', arcsec_str.strip().lower())) * 60.0
    elif re.match(r'^.*(deg|degree)$', arcsec_str.strip().lower()):
        arcsec_value = float(re.sub(r'(deg|degree)$', r'', arcsec_str.strip().lower())) * 3600.0
    else:
        try:
            arcsec_value = float(arcsec_str)
        except:
            raise ValueError('Error! Could not convert the input string "%s" to a float value in units of arcsec!'%(arcsec_str))
    return arcsec_value



# 
def parse_velocity_or_frequency(input_str, only_velocity=False, only_frequency=False):
    # 
    input_str = str(input_str)
    trimmed_str = re.sub(r'[\s\^_]', r'', input_str)
    # 
    regex_velocity = re.compile(r'^([0-9eE.+-]+)(km/s|kms-1|m/s|ms-1)$', re.IGNORECASE)
    regex_frequency = re.compile(r'^([0-9eE.+-]+)(GHz|MHz|kHz|Hz)$', re.IGNORECASE)
    # 
    output_kms = np.nan
    output_GHz = np.nan
    # 
    if regex_velocity.match(trimmed_str):
        var_value = regex_velocity.sub(r'\1', trimmed_str)
        var_unit = regex_velocity.sub(r'\2', trimmed_str)
        if var_unit.lower().startswith('km'):
            output_kms = float(var_value)
        elif var_unit.lower().startswith('m'):
            output_kms = float(var_value) / 1e3
        else:
            raise Exception('Wrong input velocity unit? %r'%(input_str))
    # 
    elif regex_frequency.match(trimmed_str):
        var_value = regex_frequency.sub(r'\1', trimmed_str)
        var_unit = regex_frequency.sub(r'\2', trimmed_str)
        if var_unit.lower() == 'GHz'.lower():
            output_GHz = float(var_value)
        elif var_unit.lower() == 'MHz'.lower():
            output_GHz = float(var_value) / 1e3
        elif var_unit.lower() == 'kHz'.lower():
            output_GHz = float(var_value) / 1e6
        elif var_unit.lower() == 'Hz'.lower():
            output_GHz = float(var_value) / 1e9
        else:
            raise Exception('Wrong input frequency unit? %r'%(input_str))
    # 
    if (only_frequency and np.isnan(output_GHz)):
        raise Exception('Wrong input frequency unit? %r'%(input_str))
    elif (only_velocity and np.isnan(output_kms)):
        raise Exception('Wrong input velocity unit? %r'%(input_str))
    elif (np.isnan(output_kms) and np.isnan(output_GHz)):
        raise Exception('Wrong input velocity or frequency unit? %r'%(input_str))
    # 
    if only_frequency:
        return output_GHz
    elif only_velocity:
        return output_kms
    else:
        return output_kms, output_GHz



# 
def parse_flux(input_str):
    # 
    input_str = str(input_str)
    trimmed_str = re.sub(r'[\s\^_]', r'', input_str)
    # 
    regex_flux_with_unit = re.compile(r'^([0-9eE.+-]+)(Jy|mJy|uJy)$', re.IGNORECASE)
    regex_flux_no_unit = re.compile(r'^([0-9eE.+-]+)$', re.IGNORECASE)
    # 
    output_Jy = np.nan
    # 
    if regex_flux_with_unit.match(trimmed_str):
        var_value = regex_flux_with_unit.sub(r'\1', trimmed_str)
        var_unit = regex_flux_with_unit.sub(r'\2', trimmed_str)
        if var_unit.lower() == 'Jy'.lower():
            output_Jy = float(var_value)
        elif var_unit.lower() == 'mJy'.lower():
            output_Jy = float(var_value) / 1e3
        elif var_unit.lower() == 'uJy'.lower():
            output_Jy = float(var_value) / 1e6
        else:
            raise Exception('Wrong input flux unit? %r'%(input_str))
    elif regex_flux_no_unit.match(trimmed_str):
        var_value = regex_flux_no_unit.sub(r'\1', trimmed_str)
        output_Jy = float(var_value)
    # 
    if np.isnan(output_Jy):
        raise Exception('Wrong input flux unit? %r'%(input_str))
    # 
    return output_Jy



# 
def parse_angle(input_str):
    # 
    input_str = str(input_str)
    trimmed_str = re.sub(r'[\s\^_]', r'', input_str)
    # 
    regex_angle_with_unit = re.compile(r'^([0-9eE.+-]+)(arcsec|arcmin|degree|deg)$', re.IGNORECASE)
    regex_angle_no_unit = re.compile(r'^([0-9eE.+-]+)$', re.IGNORECASE)
    # 
    output_arcsec = np.nan
    # 
    if regex_angle_with_unit.match(trimmed_str):
        var_value = regex_angle_with_unit.sub(r'\1', trimmed_str)
        var_unit = regex_angle_with_unit.sub(r'\2', trimmed_str)
        if var_unit.lower().startswith('arcsec'):
            output_arcsec = float(var_value)
        elif var_unit.lower().startswith('arcmin'):
            output_arcsec = float(var_value) * 60.
        elif var_unit.lower().startswith('deg'):
            output_arcsec = float(var_value) * 3600.
        else:
            raise Exception('Wrong input angle unit? %r'%(input_str))
    elif regex_angle_no_unit.match(trimmed_str):
        var_value = regex_angle_no_unit.sub(r'\1', trimmed_str)
        output_arcsec = float(var_value)
    # 
    if np.isnan(output_arcsec):
        raise Exception('Wrong input angle unit? %r'%(input_str))
    # 
    return output_arcsec



# 
def parse_skycoord(input_str):
    # 
    input_str = str(input_str).strip()
    # 
    regex_skycoord = re.compile(r'^(J2000)\s+([^\s]+)\s+([^\s]+)$', re.IGNORECASE)
    # 
    output_RA_deg = np.nan
    output_Dec_deg = np.nan
    # 
    if regex_skycoord.match(input_str):
        var_1 = regex_skycoord.sub(r'\1', input_str)
        var_2 = regex_skycoord.sub(r'\2 \3', input_str)
        cs2 = cs()
        this_skycoord = cs2.newcoordsys(direction=True)
        this_skycoord.setdirection(refcode=var_1, proj='SIN', projpar=[0,0], refpix=[0,0], refval=var_2)
        output_RA_rad, output_Dec_rad = this_skycoord.referencevalue(format='n')['numeric']
        this_skycoord.done()
        # see -- https://casa.nrao.edu/docs/CasaRef/coordsys-Tool.html
        output_RA_deg = np.rad2deg(output_RA_rad)
        output_Dec_deg = np.rad2deg(output_Dec_rad)
    # 
    if np.isnan(output_RA_deg) or np.isnan(output_Dec_deg):
        raise Exception('Wrong input coordsys format? %r'%(input_str))
    # 
    return output_RA_deg, output_Dec_deg



# 
def parse_config(input_str, return_max_baseline = True):
    # 
    input_str = str(input_str)
    # 
    output_config_file = ''
    # 
    regex_config_1 = re.compile(r'^C([0-9\*]*)-([0-9]+)$')
    regex_config_2 = re.compile(r'^alma\.cycle([0-9\*]+)\.([0-9]+)\.cfg$')
    # 
    if regex_config_1.match(input_str):
        #antenna_number = regex_config_1.sub(r'\1', input_str)
        config_number = int(regex_config_1.sub(r'\2', input_str))
        output_config_file = 'alma.cycle*.%d.cfg'%(config_number)
    # 
    if regex_config_2.match(input_str):
        cycle_number = int(regex_config_2.sub(r'\1', input_str))
        config_number = int(regex_config_2.sub(r'\2', input_str))
        output_config_file = 'alma.cycle*.%d.cfg'%(config_number)
    # 
    if output_config_file == '':
        raise Exception('Wrong input config format? %r'%(input_str))
    # 
    # find the latest cycle config file
    global casa_simmos_dir
    found_configs = glob.glob(casa_simmos_dir+os.sep+output_config_file)
    if found_configs is None or len(found_configs) == 0:
        raise Exception('Error! Cound not find '+casa_simmos_dir+os.sep+output_config_file)
    found_configs.sort(key=LooseVersion)
    output_config_file = os.path.basename(found_configs[-1])
    # 
    if return_max_baseline:
        config_baseline_dict = {1: 0.16, 2: 0.31, 3: 0.50, 4: 0.78, 5: 1.4, 6: 2.5, 7: 3.6, 8: 8.5, 9: 13.9, 10: 16.2}
        max_baseline = config_baseline_dict[config_number]
        return output_config_file, max_baseline
    # 
    return output_config_file



# 
def calc_beamsize_for_antennalist(antennalist, skyfreq_GHz):
    global casa_simmos_dir
    x_list = []
    y_list = []
    diam_list = []
    with open(casa_simmos_dir+os.sep+antennalist, 'r') as fp:
        for textline in fp:
            if not textline.startswith('#') and not textline.strip() == '':
                x, y, z, diam, antname = textline.split()[0:5]
                x_list.append(float(x))
                y_list.append(float(y))
                diam_list.append(float(diam))
    # 
    x = np.array(x_list)
    y = np.array(y_list)
    diam = np.mean(diam_list)
    # 
    uvdist = []
    for i in range(len(x)):
        for j in range(len(y)):
            if j != i:
                uvdist.append(np.sqrt((x[i]-x[j])**2 + (y[i]-y[j])**2))
    # 
    ref_freq_Hz = skyfreq_GHz * 1e9
    print2('calc_beamsize_for_antennalist: sky freq. = %s [GHz]'%(skyfreq_GHz))
    # 
    print2('calc_beamsize_for_antennalist: ant. num. = %d'%(len(x)))
    maxuvdist = np.max(uvdist)
    print2('calc_beamsize_for_antennalist: maxuvdist = %s [m]'%(maxuvdist))
    L80uvdist = np.percentile(uvdist, 80) # np.max(uvdist) # now I am using 90-th percentile of baselies, same as used by 'analysisUtils.py' pickCellSize() getBaselineStats(..., percentile=...)
    print2('calc_beamsize_for_antennalist: L80uvdist = %s [m] (80-th percentile)'%(L80uvdist))
    synbeam = 2.99792458e8 / ref_freq_Hz / maxuvdist / np.pi * 180.0 * 3600.0 # arcsec
    print2('calc_beamsize_for_antennalist: syn. beam = %s [arcsec] (max baseline)'%(synbeam))
    synbeam = 0.574 * 2.99792458e8 / ref_freq_Hz / L80uvdist / np.pi * 180.0 * 3600.0 # arcsec # .574lambda/L80, see 'analysisUtils.py' estimateSynthesizedBeamFromASDM()
    print2('calc_beamsize_for_antennalist: syn. beam = %s [arcsec] (80-th percentile baselines)'%(synbeam))
    # 
    return synbeam




# 
# def dzliu_simobs
# 
def dzliu_simobs(
        project = 'simobs', 
        totaltime = '30min', 
        integration = '30s', 
        source_shape = 'gaussian', 
        source_size = None, 
        source_total_flux = None, 
        source_center = None, 
        line_shape = 'gaussian', 
        line_FWHM = None, 
        line_center = None, 
        channel_width = None, 
        pixel_size = None, 
        beam_size = None, 
        config = None, 
        skymodel = None, 
        skyfreq = None, 
        setpointings = True, 
        ptgfile = None, 
        phasecenter = None, 
        thermalnoise = 'tsys-atm', 
        verbose = True, 
        overwrite = True, 
    ):
    # 
    casalog.origin('dzliu_simobs')
    # 
    # set default values
    if skyfreq is None:
        skyfreq = '230GHz'
        print2('Setting skyfreq to default %r.'%(skyfreq))
    elif type(skyfreq) is not str:
        raise Exception('Error! The type of skyfreq %s must be str.'%(type(skyfreq)))
    # 
    skyfreq_GHz = parse_velocity_or_frequency(skyfreq, only_frequency=True)
    # 
    if source_size is None:
        source_size = [0.5, 0.5, 0.0] # arcsec, arcsec, deg
        print2('Setting source_size to default %s x %s arcsec (PA %s deg).'%(source_size[0],source_size[1],source_size[2]))
    elif type(source_size) not in [list, tuple, np.ndarray]:
        raise Exception('Error! The type of source_size %s must be [list, tuple, np.ndarray].'%(type(source_size)))
    elif len(source_size) != 3:
        raise Exception('Error! The length of source_size %s must be 3 [major_FWHM_in_arcsec, minor_FWHM_in_arcsec, PA_in_degree].'%(str(source_size)))
    # 
    if source_total_flux is None:
        source_total_flux = '1.0Jy'
        print2('Setting source_total_flux to default %r.'%(source_total_flux))
    elif type(source_total_flux) is not str:
        raise Exception('Error! The type of source_total_flux %s must be str.'%(type(source_total_flux)))
    # 
    source_total_flux_Jy = parse_flux(source_total_flux)
    # 
    if config is None:
        if beam_size is not None:
            raise NotImplementedError('TODO: calculate baseline length from beam_size_arcsec and skyfreq_GHz then find the config.')
        config = 'C-5'
        print2('Setting config to default %r.'%(config))
    # 
    antennalist, max_baseline = parse_config(config, return_max_baseline=True)
    # 
    beam_size_arcsec = calc_beamsize_for_antennalist(antennalist, skyfreq_GHz)
    # 
    #if beam_size is None:
    #    beam_size = [0.1, 0.1, 0.0] # arcsec, arcsec, deg
    #    print2('Setting beam_size to default %s x %s arcsec (PA %s deg).'%(beam_size[0],beam_size[1],beam_size[2]))
    #elif type(beam_size) not in [list, tuple, np.ndarray]:
    #    raise Exception('Error! The type of beam_size %s must be [list, tuple, np.ndarray].'%(type(beam_size)))
    #elif len(beam_size) != 3:
    #    raise Exception('Error! The length of beam_size %s must be 3 [major_FWHM_in_arcsec, minor_FWHM_in_arcsec, PA_in_degree].'%(str(beam_size)))
    # 
    if pixel_size is None:
        oversampling = 5.0
        pixel_size_arcsec = beam_size_arcsec/oversampling # arcsec
        pixel_size_arcsec = int(pixel_size_arcsec*1000.0)/1000.0 # arcsec, round to 1mas precision
        pixel_size = '%sarcsec'%(pixel_size_arcsec)
        print2('Setting pixel_size to default 1/%s beam_size %s.'%(oversampling, pixel_size))
    elif type(pixel_size) not in [float, np.float32, np.float64]:
        raise Exception('Error! The type of pixel_size %s must be [float, np.float32, np.float64].'%(type(pixel_size)))
    # 
    pixel_size_arcsec = parse_angle(pixel_size)
    # 
    if line_FWHM is None:
        line_FWHM = '400km/s'
        print2('Setting line_FWHM to default %r.'%(line_FWHM))
    elif type(line_FWHM) is not str:
        raise Exception('Error! The type of line_FWHM %s must be str.'%(type(line_FWHM)))
    # 
    line_FWHM_kms, line_FWHM_GHz = parse_velocity_or_frequency(line_FWHM)
    # 
    if phasecenter is None:
        phasecenter = 'J2000 150.0deg 2.0deg'
        print2('Setting phasecenter to default %r.'%(phasecenter))
    elif type(phasecenter) is not str:
        raise Exception('Error! The type of phasecenter %s must be str.'%(type(phasecenter)))
    # 
    if source_center is None:
        source_center = phasecenter
        print2('Setting source_center to default phasecenter %r.'%(source_center))
    elif type(source_center) is not str:
        raise Exception('Error! The type of source_center %s must be str.'%(type(source_center)))
    # 
    source_center_RA_deg, source_center_Dec_deg = parse_skycoord(source_center)
    # 
    if line_center is None:
        line_center = skyfreq
        print2('Setting line_center to default skyfreq %r.'%(line_center))
    elif type(line_center) is not str:
        raise Exception('Error! The type of line_center %s must be str.'%(type(line_center)))
    # 
    line_center_kms, line_center_GHz = parse_velocity_or_frequency(line_center)
    # 
    if channel_width is None:
        channel_width = '8.0MHz'
        print2('Setting channel_width to default %r.'%(channel_width))
    elif type(channel_width) is not str:
        raise Exception('Error! The type of channel_width %s must be str.'%(type(channel_width)))
    # 
    channel_width_kms, channel_width_GHz = parse_velocity_or_frequency(channel_width)
    # 
    # 
    # check velocity frequency definition
    if np.isnan(line_FWHM_GHz):
        line_FWHM_GHz = line_FWHM_kms/2.99792458e5 * skyfreq_GHz
    if np.isnan(line_center_GHz):
        line_center_GHz = line_center_kms/2.99792458e5 * skyfreq_GHz
    if np.isnan(channel_width_GHz):
        channel_width_GHz = channel_width_kms/2.99792458e5 * skyfreq_GHz
    # 
    # 
    # make skymodel
    if skymodel is not None:
        print2('Using user input skymodel %r.'%(skymodel))
        if skymodel.endswith('.fits'):
            skymodel_fitsfile = skymodel
        else:
            skymodel_fitsfile = skymodel+'.fits'
            if not os.path.isfile(skymodel_fitsfile):
                exportfits(skymodel, skymodel_fitsfile)
                if not os.path.isfile(skymodel_fitsfile):
                    raise Exception('Error! Could not run exportfits!')
        skymodel_data = pyfits.getdata(skymodel_fitsfile)
    else:
        print2('Making skymodel with source_shape %r and source_size %s'%(source_shape, source_size))
        skymodel = project+'.skymodel.image'
        skymodel_fitsfile = project+'.skymodel.image.fits'
        print2('os.getcwd(): %r'%(os.getcwd()))
        if os.path.isfile(skymodel_fitsfile):
            print('Reading %r'%(skymodel_fitsfile))
            skymodel_data = pyfits.getdata(skymodel_fitsfile)
        else:
            fwhm = [source_size[0]/pixel_size_arcsec, source_size[1]/pixel_size_arcsec]
            imsize = int(np.ceil(15.0*np.max(fwhm))) #<TODO># 15*FWHM
            ygrid, xgrid = np.mgrid[1:imsize+1, 1:imsize+1] # 1-based
            xyarray = np.column_stack([xgrid.ravel(), ygrid.ravel()])
            fn2 = fn()
            skymodel_data_2D_fn = fn2.gaussian2d(amplitude=1.0, center=[(imsize+1)/2.0, (imsize+1)/2.0], 
                                                 fwhm=fwhm, pa='%sdeg'%(source_size[2]))
            # -- see https://casadocs.readthedocs.io/en/stable/_modules/casatools/functional.html#functional.gaussian2d
            skymodel_data_2D = skymodel_data_2D_fn.f(xyarray.ravel())
            skymodel_data_2D = skymodel_data_2D / np.nansum(skymodel_data_2D) * source_total_flux_Jy
            skymodel_data_2D = skymodel_data_2D.reshape(imsize, imsize)
            # 
            spectral_FWHM = line_FWHM_GHz / channel_width_GHz
            spectral_size = int(np.ceil(7.0*spectral_FWHM)) #<TODO># 15*FWHM
            spectral_grid = np.arange(spectral_size)+1.0 # 1-based
            spectral_data_1D_fn = fn2.gaussian1d(amplitude=1.0, center=(spectral_size+1)/2.0, 
                                                 fwhm=spectral_FWHM)
            # -- see https://casadocs.readthedocs.io/en/stable/_modules/casatools/functional.html#functional.gaussian1d
            spectral_data_1D = spectral_data_1D_fn.f(spectral_grid.ravel())
            # 
            skymodel_data = np.full((spectral_size, imsize, imsize), fill_value = 0.0)
            for ichan in range(spectral_size):
                skymodel_data[ichan, :, :] = skymodel_data_2D[:, :] * spectral_data_1D[ichan]
            # 
            # savefits
            hdu = pyfits.PrimaryHDU(skymodel_data)
            hdu.header['BUNIT'] = 'Jy/pixel'
            hdu.header['EQUINOX'] = 2000.0
            hdu.header['RADESYS'] = 'FK5'
            hdu.header['SPECSYS'] = 'LSRK'
            hdu.header['CTYPE1'] = 'RA---SIN'
            hdu.header['CTYPE2'] = 'DEC--SIN'
            hdu.header['CTYPE3'] = 'FREQ' # 'VRAD'
            hdu.header['CUNIT1'] = 'deg'
            hdu.header['CUNIT2'] = 'deg'
            hdu.header['CUNIT3'] = 'Hz'
            hdu.header['CDELT1'] = -pixel_size_arcsec/3600.0
            hdu.header['CDELT2'] = pixel_size_arcsec/3600.0
            hdu.header['CDELT3'] = channel_width_GHz*1e9
            hdu.header['CRVAL1'] = source_center_RA_deg
            hdu.header['CRVAL2'] = source_center_Dec_deg
            hdu.header['CRVAL3'] = line_center_GHz*1e9
            hdu.header['CRPIX1'] = (imsize+1.0)/2.0
            hdu.header['CRPIX2'] = (imsize+1.0)/2.0
            hdu.header['CRPIX3'] = (spectral_size+1.0)/2.0
            hdu.header['RESTFRQ'] = (skyfreq_GHz*1e9, 'Hz')
            hdu.header['SKYFREQ'] = (skyfreq_GHz, 'GHz')
            hdu.header['PIXSC'] = (pixel_size_arcsec, 'arcsec')
            hdu.header['CHANWID'] = (channel_width_GHz*1e3, 'MHz')
            hdu.header['LINEFWHM'] = (line_FWHM_GHz*1e3, 'MHz')
            hdu.header['SFLUX'] = (source_total_flux_Jy, 'Jy, source total flux')
            hdu.header['SMAJ'] = (source_size[0], 'arcsec, source major FWHM')
            hdu.header['SMIN'] = (source_size[1], 'arcsec, source minor FWHM')
            hdu.header['SPA'] = (source_size[2], 'deg, source PA FWHM')
            hdu.header['BMAJ'] = (0.0001/3600.0, 'arcsec, beam major FWHM')
            hdu.header['BMIN'] = (0.0001/3600.0, 'arcsec, beam minor FWHM')
            hdu.header['BPA'] = (0.0, 'deg, beam PA FWHM')
            #hdulist = pyfits.HDUList([hdu])
            #print2(str(hdu.header))
            #print2(str(hdu.data))
            #print2('hdulist: %s'%(hdulist))
            if skymodel_fitsfile.find(os.sep)>=0:
                if not os.path.isdir(os.path.dirname(skymodel_fitsfile)):
                    os.makedirs(os.path.dirname(skymodel_fitsfile))
            if os.path.isfile(skymodel_fitsfile):
                shutil.move(skymodel_fitsfile, skymodel_fitsfile+'.backup')
            #with open(skymodel_fitsfile, 'ab+') as ofp:
            #    hdulist.writeto(ofp, clobber=True, checksum=True)
            hdu.writeto(skymodel_fitsfile, clobber=True, checksum=True)
            #hdulist.close()
            print('Written to %r'%(skymodel_fitsfile))
        # 
        # import fits
        #if os.path.isdir(skymodel):
        #    shutil.move(skymodel, skymodel+'.backup')
        #importfits(skymodel_fitsfile, skymodel)
        # 
        skymodel = skymodel_fitsfile
    # 
    # 
    # run simobserve
    output_noisy_ms = project+'.'+antennalist.replace('.txt','')+'.noisy.ms'
    if not os.path.isdir(project+os.sep+output_noisy_ms):
        
        print('Running simobserve ...')
        
        if os.path.isdir(project):
            shutil.rmtree(project)
        
        simobserve_params = OrderedDict()
        
        # project
        simobserve_params['project'] = project
        
        # skymodel
        simobserve_params['skymodel'] = skymodel
        
        # inbright
        simobserve_params['inbright'] = '' # '%sJy/pixel'%(source_total_flux / np.nansum(skymodel_data) * np.nanmax(skymodel_data)) # 'Jy/pixel' of the brightest pixel
        
        # indirection
        simobserve_params['indirection'] = source_center # source_center
        
        # incenter
        simobserve_params['incenter'] = ''
        
        # inwidth
        simobserve_params['inwidth'] = '%sMHz'%(channel_width_GHz*1e3)
        
        # phasecenter
        simobserve_params['direction'] = phasecenter
        
        # integration
        simobserve_params['integration'] = integration
        
        # totaltime
        simobserve_params['totaltime'] = totaltime
        
        # maptype
        simobserve_params['maptype'] = 'hexagonal'
        
        # obsmode
        simobserve_params['obsmode'] = 'int'
        
        # antennalist
        simobserve_params['antennalist'] = antennalist # 'alma;%s'%(beam_size)
        
        # sdant and sdantlist
        simobserve_params['sdant'] = 0
        simobserve_params['sdantlist'] = '' # 'aca.tp.cfg'
        
        # thermalnoise
        simobserve_params['thermalnoise'] = thermalnoise
        
        # others
        simobserve_params['overwrite'] = overwrite
        simobserve_params['verbose'] = verbose
        
        # print and write simobserve_params to file
        simobserve_params_str = 'simobserve('+', '.join("{!s}={!r}".format(k, simobserve_params[k]) for k in simobserve_params.keys())+')'
        print2(simobserve_params_str)
        
        simobserve(**simobserve_params)
        
        #listobs('casasim_Project/casasim_Project.vla.a.ms')
        
        #if os.path.isdir('casasim_Project/casasim_Project.vla.a.noisy.ms'):
        #    os.system('rm -rf casasim_Project/casasim_Project.vla.a.noisy.ms')
    
        # 
        # Manually add noise and "corrupt" the visibilities, 
        # see https://casa.nrao.edu/docs/CasaRef/simulator-Tool.html#simulator.setnoise.html
        # 
        #if not os.path.isdir('casasim_Project/casasim_Project.vla.a.noisy.ms'):
        #    
        #    print('')
        #    print('Adding noise ...')
        #    
        #    os.system('rm -rf casasim_Project/casasim_Project.vla.a.noisy.ms')
        #    
        #    os.system('cp -r casasim_Project/casasim_Project.vla.a.ms casasim_Project/casasim_Project.vla.a.noisy.ms')
        #    
        #    sm.openfromms('casasim_Project/casasim_Project.vla.a.noisy.ms')
        #    #sm.predict('')
        #    sm.setnoise(mode='simplenoise', simplenoise='0.78e-3 Jy')
        #    #sm.setpa(mode='calculate')
        #    sm.corrupt()
        #    sm.done()
        #    
        #    if os.path.isdir('output_dirty_image_3000.image'):
        #        os.system('rm -rf output_dirty_image_3000.*')
        
    else:
        
        print2('Found existing %r. Will not overwrite it.'%(project+os.sep+output_noisy_ms))
    
    
    # 
    # Run tclean
    # 
    













############
#   main   #
############

dzliu_main_func_name = 'dzliu_simobs' # make sure this is the right main function in this script file

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




