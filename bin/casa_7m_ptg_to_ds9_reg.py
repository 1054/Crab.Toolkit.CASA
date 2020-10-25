#!/usr/bin/env python
# 
from __future__ import print_function
import os, sys, re, json, copy
import astropy
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle, FK5
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from regions import CircleSkyRegion, read_ds9, write_ds9, PixCoord
import numpy as np
#import matplotlib as mpl
#import matplotlib.pyplot as plt
sqrt = np.sqrt
ln = np.log
pi = np.pi

# 
# Read user input
# 
if len(sys.argv) <= 3:
    print('Usage:')
    print('    casa_7m_ptg_to_ds9_reg.py "scan_pointings_7m.ptg.txt" "output_ds9.reg" sky_frequency_GHz')
    print('    # Note that sky_frequency_GHz is a float value in units of GHz')
    sys.exit()
casa_ptg_file = sys.argv[1]
output_ds9_regions_file = sys.argv[2]
sky_frequency = float(sys.argv[3])


# Calc ALMA primary beam
sky_wavelength = 2.99792458e5/sky_frequency # um
primary_beam_diam_7m = 7.0 # ALMA 12m
primary_beam_diam_12m = 12.0 # ALMA 12m
#primary_beam_tape = 10.0 # https://safe.nrao.edu/wiki/bin/view/ALMA/AlmaPrimaryBeamCorrection
#primary_beam_bpar = 1.243 - 0.343 * primary_beam_tape + 0.12 * primary_beam_tape**2 # http://legacy.nrao.edu/alma/memos/html-memos/alma456/memo456.pdf -- Eq(18)
primary_beam_bpar = 1.13
primary_beam_fwhm_7m = primary_beam_bpar * sky_wavelength / (primary_beam_diam_7m*1e6) / pi * 180.0 * 3600.0 # arcsec
primary_beam_fwhm_12m = primary_beam_bpar * sky_wavelength / (primary_beam_diam_12m*1e6) / pi * 180.0 * 3600.0 # arcsec

primary_beam_fwhm = primary_beam_fwhm_7m


pointing_table = Table.read(casa_ptg_file, format='ascii.commented_header')
#print(pointing_table)

# Fix CASA Dec format
if re.match(r'^([^\.]+)\.([^\.]+)\.([^\.]+)\.([^\.]+)$', pointing_table['DEC'][0]):
    pointing_table['DEC'] = [re.sub(r'^([^\.]+)\.([^\.]+)\.([^\.]+)\.([^\.]+)$', r'\1:\2:\3.\4', t) for t in pointing_table['DEC']]
    # CASA Dec format is like [+-]XX.XX.XX.XXXX

scan_grid_ra_dec = ['%s %s'%(str(pointing_table['RA'][i]), str(pointing_table['DEC'][i])) for i in range(len(pointing_table))]
#print(scan_grid_ra_dec[0])
scan_grid_skycoords = SkyCoord(scan_grid_ra_dec, frame=FK5, unit=(u.hourangle, u.deg), obstime="J2000")

radius = Angle(primary_beam_fwhm/3600.0/2.0, 'deg')
region = [CircleSkyRegion(scan_skycoord, radius) for scan_skycoord in scan_grid_skycoords]
write_ds9(region, output_ds9_regions_file)
print('Output to "%s"!'%(output_ds9_regions_file))



