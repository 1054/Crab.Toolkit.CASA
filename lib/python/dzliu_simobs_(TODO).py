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
#casalog.filter('DEBUG2')
# 
#import resource
#resource.setrlimit(resource.RLIMIT_AS, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
#resource.setrlimit(resource.RLIMIT_NOFILE, (10240, 40960)) # ulimit -Sn # ulimit -Hn
# 
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
else:
    # see CASA 6 updates here: https://alma-intweb.mtk.nao.ac.jp/~eaarc/UM2018/presentation/Nakazato.pdf
    from casatasks import tclean, mstransform, exportfits, concat, split, imstat
    #from casatasks import sdbaseline
    #from casatools import ia
# 
if not os.path.isfile(os.getenv('CASAPATH').split(' ')[0]+'/data/alma/simmos/vla.a.cfg'):
    raise Exception('Error! VLA config file not found from "$CASAPATH/data/alma/simmos/vla.a.cfg"!')




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
# def dzliu_simobs
# 
def dzliu_simobs(project = 'simobs', skymodel = '', beam = '0.1arcsec', skyfreq = '230GHz'):
    # 
    casalog.origin('dzliu_simobs')
    # 
    if not os.path.isdir('casasim_Project/casasim_Project.vla.a.ms'):
        
        print('Running simobserve ...')
        
        if os.path.isdir(project):
            shutil.rmtree(project)
        
        #project = ''
        #totaltime = '2' # repeating twice 'integration' for each pointing # '120h'
        #beam = '0.75arcsec'
        #skyfreq = '3GHz'
        
        #skymodel = ''  # flux units should be Jy/pixel
        #inbright = '0.1e-6' # Jy/pixel
        #indirection = 'J2000 10h00m28.5867s +02d12m21.200s' # reset 
        #incenter = '3GHz' # reset
        #inwidth = '2048MHz' # reset
        #incell = '0.1arcsec' # (arcsec2float(beamsize)/5.0)
        
        setpointings = False # set to False so that our own ptfile is used
        ptgfile = 'casasim_Project.ptg.txt'
        direction = 'J2000 10h00m28.5867s +02d12m21.200s' # mosaic center
        mapsize = ''
        integration = '10min' # time per pointing. 
        totaltime = '2' # repeating 'integration' for each pointing
        maptype = 'hexagonal'
        obsmode = 'int'
        #antennalist = 'vla.a.cfg'
        antennalist = 'alma;%s'%(beamsize)
        hourangle = 'transit'
        
        thermalnoise = ''
        
        overwrite = True
        
        inp(simobserve)
        
        simobserve()
        
        listobs('casasim_Project/casasim_Project.vla.a.ms')
        
        if os.path.isdir('casasim_Project/casasim_Project.vla.a.noisy.ms'):
            os.system('rm -rf casasim_Project/casasim_Project.vla.a.noisy.ms')
    
    # 
    # Manually add noise and "corrupt" the visibilities, 
    # see https://casa.nrao.edu/docs/CasaRef/simulator-Tool.html#simulator.setnoise.html
    # 
    if not os.path.isdir('casasim_Project/casasim_Project.vla.a.noisy.ms'):
        
        print('')
        print('Adding noise ...')
        
        os.system('rm -rf casasim_Project/casasim_Project.vla.a.noisy.ms')
        
        os.system('cp -r casasim_Project/casasim_Project.vla.a.ms casasim_Project/casasim_Project.vla.a.noisy.ms')
        
        sm.openfromms('casasim_Project/casasim_Project.vla.a.noisy.ms')
        #sm.predict('')
        sm.setnoise(mode='simplenoise', simplenoise='0.78e-3 Jy')
        #sm.setpa(mode='calculate')
        sm.corrupt()
        sm.done()
        
        if os.path.isdir('output_dirty_image_3000.image'):
            os.system('rm -rf output_dirty_image_3000.*')













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




