#!/usr/bin/env python
# 
from __future__ import print_function
import os, sys, re, json, copy
import astropy
from astropy.io import fits
import numpy as np

# 
# Read user input
# 
source_position = [100, 100] # 0-based
source_morphology = [50.0, 50.0, 0.0] # major axis radius, minor axis radius and position angle in units of pixel, pixel and degree, respectively. 
#pixel_scale = [0.01, 0.01] # arcsec
image_size = [301, 301]
naxis1, naxis2 = image_size
output_image = 'skymodel_flat_disc.fits'
output_json = 'skymodel_flat_disc.json'

image_data = np.full((naxis2, naxis1), 0.0)
gridy, gridx = np.mgrid[0:naxis2, 0:naxis1]
discx, discy = source_position
theta_angle = np.arctan2((gridy-discy), (gridx-discx)) # in radian
theta_rota = theta_angle - np.deg2rad(source_morphology[2]+90.0)
print('theta_rota[0, 0]', np.rad2deg(theta_rota[0, 0]))
source_mask = ( ((gridx-discx)**2 + (gridy-discy)**2) < ((source_morphology[0]*np.cos(theta_rota))**2 + (source_morphology[1]*np.sin(theta_rota))**2) )
image_data[source_mask] = 1.0

hdu = fits.PrimaryHDU(image_data)
#hdu.header['PIXSCALE'] = pixel_scale[1]
#hdu.header.comments['PIXSCALE'] = 'in units of arcsec'
hdu.writeto(output_image, overwrite=True)
print('Output to "%s"!'%(output_image))

output_info = {}
output_info['source_position'] = source_position
output_info['source_morphology'] = source_morphology
#output_info['pixel_scale'] = pixel_scale
output_info['image_size'] = image_size
with open(output_json, 'w') as jfp:
    json.dump(output_info, jfp, indent=4)
print('Output to "%s"!'%(output_json))


