#!/usr/bin/env python
# 
# Run this code with python3
# 
# Note that the beams are different among images!
# 

import os, sys, re, glob, shutil, copy
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.coordinates import SkyCoord, FK5
from regions import DS9Parser, read_ds9, write_ds9, RectangleSkyRegion, CircleSkyRegion, PixelRegion, LinePixelRegion, ds9_objects_to_string
from reproject import reproject_interp
from scipy.interpolate import griddata


# User-defined parameters
Region_files = glob.glob('ds9_regions_of_mosaic_pointings_in_*_boxes_unpadded.reg')
DataSet_names = [re.sub(r'ds9_regions_of_mosaic_pointings_in_(.*)_boxes_unpadded.reg', r'\1', os.path.basename(t)) for t in Region_files] # [os.path.basename(os.path.dirname(os.path.dirname((t)))) for t in DataSet_dirs]


# Check input
if len(Region_files) == 0:
    print('Error! File "ds9_regions_of_mosaic_pointings_in_*_boxes_unpadded.reg" not found! Please run \'a_dzliu_code_print_phasecenters_in_casa.py\' or \'a_dzliu_code_run_dividing_mosaic_and_imaging_in_casa.py\' first!')
    sys.exit(255)


# define functions
def pop_out_obsgeo(header):
    header_obsgeo_x = None
    if 'OBSGEO-X' in header:
        header_obsgeo_x = header['OBSGEO-X']
        del header['OBSGEO-X']
    header_obsgeo_y = None
    if 'OBSGEO-Y' in header:
        header_obsgeo_y = header['OBSGEO-Y']
        del header['OBSGEO-Y']
    header_obsgeo_z = None
    if 'OBSGEO-Z' in header:
        header_obsgeo_z = header['OBSGEO-Z']
        del header['OBSGEO-Z']
    return header, header_obsgeo_x, header_obsgeo_y, header_obsgeo_z

def push_in_obsgeo(header, header_obsgeo_x, header_obsgeo_y, header_obsgeo_z):
    if header_obsgeo_x is not None:
        header['OBSGEO-X'] = header_obsgeo_x
    if header_obsgeo_y is not None:
        header['OBSGEO-Y'] = header_obsgeo_y
    if header_obsgeo_z is not None:
        header['OBSGEO-Z'] = header_obsgeo_z
    return header


# Start processing
print('getcwd: %s'%(os.getcwd()))


# loop datasets
for Region_file, DataSet_name in list(zip(Region_files, DataSet_names)):
    
    print('DataSet: '+DataSet_name)

    # read ds9 regions
    list_of_regions = read_ds9(Region_file)
    
    # prepare lists
    list_of_icol = []
    list_of_irow = []
    list_of_image_file = [] # science image
    list_of_pb_file = [] # primary beam image
    list_of_corner_coords = [] # science image x0y0 (0,0) x1y1 (naxis1-1,naxis2-1)
    
    # prepare output header
    output_image = None
    output_weight = None
    output_coverage = None
    output_header = None
    output_pixsc = None
    
    # loop regions to find the overall extent
    for region_idx, region_obj in enumerate(list_of_regions):
        if isinstance(region_obj, RectangleSkyRegion):
            #print(region_obj)
            #print(region_obj.meta['text'])
            icol = int(re.sub(r'divided mosaic ([0-9]+) ([0-9]+)', r'\1', region_obj.meta['text']))
            irow = int(re.sub(r'divided mosaic ([0-9]+) ([0-9]+)', r'\2', region_obj.meta['text']))
            print('i %d, icol irow: %d %d'%(region_idx, icol, irow))
            image_files = glob.glob('Level_4_Data_Images_Divide_Mosaic/%s_Mosaic_%d_%d/output_*_dirty.image.pbcor.fits'%(DataSet_name, icol, irow))
            if len(image_files) > 0:
                image_file = image_files[0]
                pb_file = image_file.replace('.image.pbcor.fits', '.pb.fits')
                # 
                #image, header = fits.getdata(image_file, header=True)
                #pb = fits.getdata(pb_file)
                header = fits.getheader(image_file)
                wcs = WCS(header, naxis=2)
                pixsc = proj_plane_pixel_scales(wcs)[-1] * 3600.0 # arcsec
                print('image: %r, NAXIS1 NAXIS2: %d %d'%(image_file, header['NAXIS1'], header['NAXIS2']))
                # 
                #RA_corners, Dec_corners = wcs.wcs_pix2world([0, header['NAXIS1']-1], [0, header['NAXIS2']-1], 0)
                # 
                #print(region_obj)
                RA_corners = [ region_obj.center.ra.deg + region_obj.width.to(u.deg).value / 2.0 / np.cos(np.deg2rad(region_obj.center.dec.deg)), 
                               region_obj.center.ra.deg - region_obj.width.to(u.deg).value / 2.0 / np.cos(np.deg2rad(region_obj.center.dec.deg)) ]
                Dec_corners = [ region_obj.center.dec.deg - region_obj.height.to(u.deg).value / 2.0, 
                                region_obj.center.dec.deg + region_obj.height.to(u.deg).value / 2.0 ]
                # 
                list_of_icol.append(icol)
                list_of_irow.append(irow)
                list_of_image_file.append(image_file)
                list_of_pb_file.append(pb_file)
                list_of_corner_coords.append([SkyCoord(RA_corners[0], Dec_corners[0], frame=FK5, unit=(u.deg, u.deg)),\
                                              SkyCoord(RA_corners[1], Dec_corners[1], frame=FK5, unit=(u.deg, u.deg))])
                # 
                if output_header is None:
                    output_header = copy.deepcopy(header)
                    output_pixsc = pixsc
    #raise NotImplementedError()
    
    # compute the full extent of the output image
    max_RA = np.max([np.max([t[0].ra.deg, t[1].ra.deg]) for t in list_of_corner_coords])
    min_RA = np.min([np.min([t[0].ra.deg, t[1].ra.deg]) for t in list_of_corner_coords])
    max_Dec = np.max([np.max([t[0].dec.deg, t[1].dec.deg]) for t in list_of_corner_coords])
    min_Dec = np.min([np.min([t[0].dec.deg, t[1].dec.deg]) for t in list_of_corner_coords])
    dRA = (max_RA-min_RA)*np.cos(np.deg2rad((max_Dec+min_Dec)/2.0))*3600.0 # arcsec
    dDec = (max_Dec-min_Dec)*3600.0 # arcsec
    
    # create output image
    output_header['NAXIS1'] = int(np.ceil(dRA/pixsc))
    output_header['NAXIS2'] = int(np.ceil(dDec/pixsc))
    output_header['CRPIX1'] = (output_header['NAXIS1']+1.0)/2.0
    output_header['CRPIX2'] = (output_header['NAXIS2']+1.0)/2.0
    output_header['CRVAL1'] = (max_RA+min_RA)/2.0
    output_header['CRVAL2'] = (max_Dec+min_Dec)/2.0
    output_wcs = WCS(output_header, naxis=2)
    print('Mosaic image extent dRA dDec: %s, %s, size: %d x %d, pixel size: %s'%(dRA, dDec, output_header['NAXIS1'], output_header['NAXIS2'], pixsc))
    #print('output_header', output_header)
    
    # initialize output image weight coverage arrays <TODO> assuming 3D!
    out_nchan = output_header['NAXIS3']
    out_ny = output_header['NAXIS2']
    out_nx = output_header['NAXIS1']
    output_image = np.full((out_nchan, out_ny, out_nx), fill_value=0.0)
    output_weight = np.full((out_nchan, out_ny, out_nx), fill_value=0.0)
    output_coverage = np.full((out_nchan, out_ny, out_nx), fill_value=0)
    
    # fix observer issue which makes reproject fails: SpectralCoord instance has no attribute 'in_observer_velocity_frame'
    output_header, geox, geoy, geoz = pop_out_obsgeo(output_header)
    
    # prepare a duplicate header for reproject in segments
    reproject_header = copy.deepcopy(output_header)
    
    # reproject
    for i in range(len(list_of_image_file)):
        
        icol, irow = list_of_icol[i], list_of_irow[i]
        
        if not (icol == 4 and irow == 2):
            continue
        
        image_data, image_header = fits.getdata(list_of_image_file[i], header=True)
        pb_data, pb_header = fits.getdata(list_of_pb_file[i], header=True)
        
        image_header, geox_not_used, geoy_not_used, geoz_not_used = pop_out_obsgeo(image_header)
        pb_header, geox_not_used, geoy_not_used, geoz_not_used = pop_out_obsgeo(pb_header)
        
        if len(image_data.shape) != 3: raise Exception('Error! Input image should be a 3D data!')
        if len(pb_data.shape) != 3: raise Exception('Error! Input pb should be a 3D data!')
        
        #image_data[np.isnan(image_data)] = 0.0 #<TODO># pad with zero?
        #pb_data[np.isnan(pb_data)] = 1e-10 #<TODO># pad with zero?
        ichan = 0
        image_data = image_data[ichan]
        pb_data = pb_data[ichan]
        image_header['NAXIS'] = 2
        pb_header['NAXIS'] = 2
        
        corner00_sky, corner11_sky = list_of_corner_coords[i]
        corner_pix_x, corner_pix_y = output_wcs.wcs_world2pix([corner00_sky.ra.deg, corner11_sky.ra.deg], [corner00_sky.dec.deg, corner11_sky.dec.deg], 0)
        corner00_pix_x, corner11_pix_x = corner_pix_x
        corner00_pix_y, corner11_pix_y = corner_pix_y
        x0 = int(corner00_pix_x)
        y0 = int(corner00_pix_y)
        if x0 < 0: x0 = 0
        if y0 < 0: y0 = 0
        x1 = int(corner11_pix_x)
        y1 = int(corner11_pix_y)
        if x1 > output_header['NAXIS1']-1: x1 = output_header['NAXIS1']-1
        if y1 > output_header['NAXIS2']-1: y1 = output_header['NAXIS2']-1
        nx = x1 - x0 + 1
        ny = y1 - y0 + 1
        
        
        # use reproject_interp
        reproject_header['NAXIS1'] = nx
        reproject_header['NAXIS2'] = ny
        reproject_header['CRPIX1'] = output_header['CRPIX1'] - x0
        reproject_header['CRPIX2'] = output_header['CRPIX2'] - y0
        
        image_hdu = fits.PrimaryHDU(data = image_data, header = image_header)
        pb_hdu = fits.PrimaryHDU(data = pb_data, header = pb_header)
        
        print('reprojecting image %d [%dx%d] to (%d,%d)+[%dx%d] at icol irow %d %d (progress %d/%d)'%(\
            i, image_header['NAXIS1'], image_header['NAXIS2'], x0, y0, nx, ny, icol, irow, i+1, len(list_of_image_file)))
        new_image = reproject_interp(image_hdu, reproject_header, return_footprint=False, order='nearest-neighbor')
        
        print('reprojecting pb %d [%dx%d] to (%d,%d)+[%dx%d] at icol irow %d %d (progress %d/%d)'%(\
            i, pb_header['NAXIS1'], pb_header['NAXIS2'], x0, y0, nx, ny, icol, irow, i+1, len(list_of_image_file)))
        new_pb = reproject_interp(pb_hdu, reproject_header, return_footprint=False, order='nearest-neighbor')
        
        
        # use griddata
        
        
        #print('new_image.shape', new_image.shape)
        #print('new_pb.shape', new_pb.shape)
        
        if icol == 4 and irow == 2:
            output_hdu = fits.PrimaryHDU(data = new_image, header = reproject_header)
            output_file = 'debug_new_image.fits'
            if os.path.isfile(output_file):
                shutil.move(output_file, output_file+'.backup')
            output_hdu.writeto(output_file)
            print('Output to "%s"'%(output_file))
            # 
            output_hdu = image_hdu
            output_file = 'debug_input_image.fits'
            if os.path.isfile(output_file):
                shutil.move(output_file, output_file+'.backup')
            output_hdu.writeto(output_file)
            print('Output to "%s"'%(output_file))
        
        output_image[:, y0:y0+ny, x0:x0+nx] += new_image * (new_pb**2)
        output_weight[:, y0:y0+ny, x0:x0+nx] += (new_pb**2)
        output_coverage[:, y0:y0+ny, x0:x0+nx] += 1

    # weight by pb**2
    mask = (output_weight>0)
    output_image[mask] = output_image[mask] / output_weight[mask]
    
    # put back OBSGEO
    output_header = push_in_obsgeo(output_header, geox, geoy, geoz)
    
    # write to fits
    output_hdu = fits.PrimaryHDU(data = output_image, header = output_header)
    output_file = 'Level_4_Data_Images_Divide_Mosaic/%s_Mosaicked.fits'%(DataSet_name)
    if os.path.isfile(output_file):
        shutil.move(output_file, output_file+'.backup')
    output_hdu.writeto(output_file)
    print('Output to "%s"'%(output_file))
    
    # write weight to fits
    output_hdu = fits.PrimaryHDU(data = output_weight, header = output_header)
    output_file = 'Level_4_Data_Images_Divide_Mosaic/%s_Mosaicked.weight.fits'%(DataSet_name)
    if os.path.isfile(output_file):
        shutil.move(output_file, output_file+'.backup')
    output_hdu.writeto(output_file)
    print('Output to "%s"'%(output_file))
    
    # write weight to fits
    output_hdu = fits.PrimaryHDU(data = output_coverage, header = output_header)
    output_file = 'Level_4_Data_Images_Divide_Mosaic/%s_Mosaicked.coverage.fits'%(DataSet_name)
    if os.path.isfile(output_file):
        shutil.move(output_file, output_file+'.backup')
    output_hdu.writeto(output_file)
    print('Output to "%s"'%(output_file))



# done
print('Done!')








