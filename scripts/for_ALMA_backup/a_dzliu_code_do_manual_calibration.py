

import os, sys, re, json, shutil, glob
sys.path.append('/software/casa/analysis_scripts/')
sys.path.append(os.path.expanduser('~/Cloud/Github/Crab.Toolkit.CASA/lib/python'))
import analysisUtils as aU
import numpy as np
import matplotlib as mpl
#mpl.use('agg')
import matplotlib.pyplot as plt
from dzliu_uvdata_utils import (check_casa, extract_field_selections, extract_line_spw_selections, flag_line_data, selfcal_continuum_data)



if not os.path.exists('uid___A002_Xb44b49_X8393.ms'):
    importasdm('../raw/uid___A002_Xb44b49_X8393.asdm.sdm', 'uid___A002_Xb44b49_X8393.ms')
if not os.path.exists('uid___A002_Xb58564_X67c.ms'):
    importasdm('../raw/uid___A002_Xb58564_X67c.asdm.sdm',  'uid___A002_Xb58564_X67c.ms')
if not os.path.exists('uid___A002_Xb5aa7c_X4c2d.ms'):
    importasdm('../raw/uid___A002_Xb5aa7c_X4c2d.asdm.sdm', 'uid___A002_Xb5aa7c_X4c2d.ms')
if not os.path.exists('uid___A002_Xb68dbd_X7818.ms'):
    importasdm('../raw/uid___A002_Xb68dbd_X7818.asdm.sdm', 'uid___A002_Xb68dbd_X7818.ms')


vis = 'uid___A002_Xb44b49_X8393.ms'

# get bandpass and phase calibrators info
calibrators = aU.getCalibrators(vis, returnDict=True)
bandpasscalname = calibrators['BANDPASS']
print('bandpasscalname: {}'.format(bandpasscalname))
phasecalname = calibrators['PHASE']
print('phasecalname: {}'.format(phasecalname))
bandpass_cal_field_id = aU.getCalibrators(vis, intent=['BANDPASS'], returnNames=False)
print('bandpass_cal_field_id: {}'.format(bandpass_cal_field_id))
phase_cal_field_id = aU.getCalibrators(vis, intent=['PHASE'], returnNames=False)
print('phase_cal_field_id: {}'.format(phase_cal_field_id))
bandpass_cal_field_id = str(bandpass_cal_field_id[0]) # '0'
phase_cal_field_id = str(phase_cal_field_id[0]) # '4'

# get science spw
spw_dict = aU.getScienceSpws(vis, intent='OBSERVE_TARGET#ON_SOURCE', returnFreqRanges=True)
spw = ','.join([str(key) for key in spw_dict])

# do bandpass calibration
## do bandpass calibrator phase calibration
gaincal(vis, caltable=vis+'.bandpass.pcal', field=bandpass_cal_field_id, solint='int', calmode='p', gaintype='G')
plotms(vis+'.bandpass.pcal', xaxis='time', yaxis='phase', gridcols=3, gridrows=3, iteraxis='antenna', 
       avgchannel='128', plotrange=[0,0,-180,180], coloraxis='spw')
# plotms(vis+'.bandpass.pcal', xaxis='time', yaxis='amp', gridcols=3, gridrows=3, iteraxis='antenna', 
#        avgchannel='128', coloraxis='corr')
# plotms(vis+'.bandpass.pcal', xaxis='time', yaxis='amp', gridcols=3, gridrows=3, iteraxis='antenna', 
#        avgchannel='128', coloraxis='spw')
# plotms(vis+'.bandpass.pcal', xaxis='time', yaxis='amp', gridcols=3, gridrows=3, iteraxis='antenna', 
#        avgchannel='128', coloraxis='spw')
## do bandpass calibration
bandpass(vis=vis, caltable=vis+'.bpcal', field=bandpass_cal_field_id, 
         solint='inf,64MHz', combine='scan', solnorm=True,
         gaintable=[vis+'.bandpass.pcal']) # solint='inf,10chan', 
plotbandpass(caltable=vis+'.bpcal',
        xaxis='chan',
        yaxis='both',
        subplot=32)
plotms(vis+'.bpcal', xaxis='chan', yaxis='amp', gridcols=3, gridrows=3, iteraxis='antenna', 
       coloraxis='spw')
plotms(vis+'.bpcal', xaxis='chan', yaxis='phase', gridcols=3, gridrows=3, iteraxis='antenna', 
       plotrange=[0,0,-180,180], coloraxis='spw')
## apply bandpass calibration
applycal(vis=vis, field='', gaintable=[vis+'.bandpass.pcal'], interp=['linear'], 
         gainfield=[bandpass_cal_field_id], applymode='calonly')
## split out bandpass calibrated ms
outputvis = re.sub(r'\.ms$', r'', vis) + '_bpcal.ms'
split(vis=vis, outputvis=outputvis, datacolumn='corrected', keepflags=False)

vis = outputvis



# do flux calibration
## get flux calibrator info
flux_cal_name = calibrators['FLUX']
print('flux_cal_name: {}'.format(flux_cal_name))
flux_cal_field_id = aU.getCalibrators(vis, intent=['FLUX'], returnNames=False)
print('flux_cal_field_id: {}'.format(flux_cal_field_id))
flux_cal_field_id = str(flux_cal_field_id[0]) # '0'
## set model for planets
setjy(vis=vis, field=flux_cal_field_id, standard='Butler-JPL-Horizons 2012', usescratch=True)
## compute phase solution for BANDPASS+FLUX+(FIRST)PHASE
## 2025-07-06 05:30:10 INFO listobs      0    none J1924-2914          19:24:51.055957 -29.14.30.12103 ICRS    0        5088960
## 2025-07-06 05:30:10 INFO listobs      1    none J2229-0832          22:29:40.084340 -08.32.54.43552 ICRS    1         609120
## 2025-07-06 05:30:10 INFO listobs      2    none Neptune             22:54:11.128660 -07.54.40.41093 ICRS    2         780840
## 2025-07-06 05:30:10 INFO listobs      3    none J1745-2900          17:45:40.036070 -29.00.28.16710 ICRS    3         609120
## 2025-07-06 05:30:10 INFO listobs      4    none J1733-3722          17:33:15.193142 -37.22.32.39451 ICRS    4        1684800
## 0: BANDPASS, 2: FLUX, 1+3: POINT, 4: PHASE
gaincal(vis=vis, caltable=vis+'.pcal.first', field='0,1,2,3,4', # '{},{},{}'.format(bandpass_cal_field_id, flux_cal_field_id, phase_cal_field_id),
        solint='inf', calmode='p', gaintype='G')

plotcal(caltable=vis+'.pcal.first', xaxis='time', yaxis='phase', subplot=331, iteration='antenna', plotrange=[0,0,-180,180])




# do phase calibration

targetname = 'sgr_a'

phasecenter = 'J2000 17h45m40.0335s -29d00m28.2302s'

field_selection = extract_field_selections(vis, phasecenter, verbose=True) # '3,73,74'



phasecalvis = vis + '.phasecal.{}.selfcal'.format(phasecalname)
result = selfcal_continuum_data(vis, phasecalvis, field=phasecalname, spw=spw, intent='CALIBRATE_PHASE#ON_SOURCE', maximsize=1000)
phasecaltable = result['caltable']

outputvis = vis + '.selfcal'
result = selfcal_continuum_data(vis, outputvis, phasecenter=phasecenter, gaintable=[phasecaltable], maximsize=1000)
caltable = result['caltable']
