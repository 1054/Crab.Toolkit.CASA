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
#from reproject import reproject_interp # always return nan in some frame
from scipy.interpolate import griddata


# User-defined parameters
Region_files = glob.glob('ds9_regions_of_mosaic_pointings_in_*_boxes_unpadded.reg')
DataSet_names = [re.sub(r'ds9_regions_of_mosaic_pointings_in_(.*)_boxes_unpadded.reg', r'\1', os.path.basename(t)) for t in Region_files] # [os.path.basename(os.path.dirname(os.path.dirname((t)))) for t in DataSet_dirs]


# Check input
if len(Region_files) == 0:
    print('Error! File "ds9_regions_of_mosaic_pointings_in_*_boxes_unpadded.reg" not found! Please run \'a_dzliu_code_print_phasecenters_in_casa.py\' or \'a_dzliu_code_run_dividing_mosaic_and_imaging_in_casa.py\' first!')
    sys.exit(255)


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
    list_of_image_corners_coords = [] # science image lower left and upper right corner, i.e., x0y0 (0,0) x1y1 (naxis1-1,naxis2-1)
    list_of_region_corners_coords = [] # unpadded region lower left and upper right corner
    list_of_corner_coords = [] # if edge region then use image corner otherwise region corner
    
    # prepare output header
    output_image = None
    output_weight = None
    output_coverage = None
    output_header = None
    output_pixsc = None
    
    # prepare to check all input image consistency, they must have the same pixsc and nchan
    pixsc = None
    nchan = None
    
    # loop regions to find the overall extent
    for region_idx, region_obj in enumerate(list_of_regions):
        if isinstance(region_obj, RectangleSkyRegion):
            #print(region_obj)
            #print(region_obj.meta['text'])
            icol = int(re.sub(r'divided mosaic ([0-9]+) ([0-9]+)', r'\1', region_obj.meta['text']))
            irow = int(re.sub(r'divided mosaic ([0-9]+) ([0-9]+)', r'\2', region_obj.meta['text']))
            print('i %d, icol irow: %d %d'%(region_idx, icol, irow))
            image_files = glob.glob('Level_4_Data_Images_Divide_Mosaic/%s_Mosaic_%d_%d/output_*_clean.image.pbcor.fits'%(DataSet_name, icol, irow))
            if len(image_files) > 0:
                image_file = image_files[0]
                pb_file = image_file.replace('.image.pbcor.fits', '.pb.fits')
                # 
                header = fits.getheader(image_file)
                print('image: %r, NAXIS1 NAXIS2: %d %d'%(image_file, header['NAXIS1'], header['NAXIS2']))
                # 
                wcs = WCS(header, naxis=2)
                pixsc1 = proj_plane_pixel_scales(wcs)[-1] * 3600.0 # arcsec
                if pixsc is None:
                    pixsc = pixsc1
                elif not np.isclose(pixsc, pixsc1, atol=1e-3, rtol=1e-3):
                    print('Error! The input image does not have the same pixsc as others (%s vs %s)!'%(pixsc1, pixsc))
                    raise Exception('Error! The input image does not have the same pixsc as others (%s vs %s)!'%(pixsc1, pixsc))
                nchan1 = int(header['NAXIS3'])
                if nchan is None:
                    nchan = nchan1
                elif nchan != nchan1:
                    print('Error! The input image does not have the same nchan as others (%s vs %s)!'%(nchan1, nchan))
                    raise Exception('Error! The input image does not have the same nchan as others (%s vs %s)!'%(nchan1, nchan))
                # 
                image_corners_RA, image_corners_Dec = wcs.wcs_pix2world([0, header['NAXIS1']-1], [0, header['NAXIS2']-1], 0)
                # 
                region_corners_RA = [ region_obj.center.ra.deg + region_obj.width.to(u.deg).value / 2.0 / np.cos(np.deg2rad(region_obj.center.dec.deg)), 
                                      region_obj.center.ra.deg - region_obj.width.to(u.deg).value / 2.0 / np.cos(np.deg2rad(region_obj.center.dec.deg)) ]
                region_corners_Dec = [ region_obj.center.dec.deg - region_obj.height.to(u.deg).value / 2.0, 
                                       region_obj.center.dec.deg + region_obj.height.to(u.deg).value / 2.0 ]
                # 
                list_of_icol.append(icol)
                list_of_irow.append(irow)
                list_of_image_file.append(image_file)
                list_of_pb_file.append(pb_file)
                list_of_image_corners_coords.append([SkyCoord(image_corners_RA[0], image_corners_Dec[0], frame=FK5, unit=(u.deg, u.deg)),\
                                                     SkyCoord(image_corners_RA[1], image_corners_Dec[1], frame=FK5, unit=(u.deg, u.deg))])
                list_of_region_corners_coords.append([SkyCoord(region_corners_RA[0], region_corners_Dec[0], frame=FK5, unit=(u.deg, u.deg)),\
                                                      SkyCoord(region_corners_RA[1], region_corners_Dec[1], frame=FK5, unit=(u.deg, u.deg))])
                # 
                if output_header is None:
                    output_header = copy.deepcopy(header)
                    output_pixsc = pixsc
    #raise NotImplementedError()
    
    # compute the full extent of the output image
    ncol = np.max(list_of_icol)+1
    nrow = np.max(list_of_irow)+1
    max_RA = np.max([np.max([t[0].ra.deg, t[1].ra.deg]) for t in list_of_image_corners_coords])
    min_RA = np.min([np.min([t[0].ra.deg, t[1].ra.deg]) for t in list_of_image_corners_coords])
    max_Dec = np.max([np.max([t[0].dec.deg, t[1].dec.deg]) for t in list_of_image_corners_coords])
    min_Dec = np.min([np.min([t[0].dec.deg, t[1].dec.deg]) for t in list_of_image_corners_coords])
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
    #output_header, geox, geoy, geoz = pop_out_obsgeo(output_header)
    
    # prepare a duplicate header for reproject in segments
    reproject_header = copy.deepcopy(output_header)
    
    # reproject
    for i in range(len(list_of_image_file)):
        
        icol, irow = list_of_icol[i], list_of_irow[i]
        print('icol, irow = %d, %d (ncol, nrow = %d, %d)'%(icol, irow, ncol, nrow))
        
        #if not (icol == 4 and irow == 2):
        #    continue
        
        image_data, image_header = fits.getdata(list_of_image_file[i], header=True)
        pb_data, pb_header = fits.getdata(list_of_pb_file[i], header=True)
        this_wcs = WCS(image_header, naxis=2)
        this_nx = image_header['NAXIS1']
        this_ny = image_header['NAXIS2']
        
        #image_header, geox_not_used, geoy_not_used, geoz_not_used = pop_out_obsgeo(image_header)
        #pb_header, geox_not_used, geoy_not_used, geoz_not_used = pop_out_obsgeo(pb_header)
        
        if len(image_data.shape) != 3: raise Exception('Error! Input image should be a 3D data!')
        if len(pb_data.shape) != 3: raise Exception('Error! Input pb should be a 3D data!')
        
        #image_data[np.isnan(image_data)] = 0.0 #<TODO># pad with zero?
        #pb_data[np.isnan(pb_data)] = 1e-10 #<TODO># pad with zero?
        
        # use reproject_interp
        #reproject_header['NAXIS1'] = nx
        #reproject_header['NAXIS2'] = ny
        #reproject_header['CRPIX1'] = output_header['CRPIX1'] - x0
        #reproject_header['CRPIX2'] = output_header['CRPIX2'] - y0
        #
        #image_hdu = fits.PrimaryHDU(data = image_data, header = image_header)
        #pb_hdu = fits.PrimaryHDU(data = pb_data, header = pb_header)
        #
        #print('reprojecting image %d [%dx%d] to (%d,%d)+[%dx%d] at icol irow %d %d (progress %d/%d)'%(\
        #    i, image_header['NAXIS1'], image_header['NAXIS2'], x0, y0, nx, ny, icol, irow, i+1, len(list_of_image_file)))
        #new_image = reproject_interp(image_hdu, reproject_header, return_footprint=False, order='nearest-neighbor')
        #
        #print('reprojecting pb %d [%dx%d] to (%d,%d)+[%dx%d] at icol irow %d %d (progress %d/%d)'%(\
        #    i, pb_header['NAXIS1'], pb_header['NAXIS2'], x0, y0, nx, ny, icol, irow, i+1, len(list_of_image_file)))
        #new_pb = reproject_interp(pb_hdu, reproject_header, return_footprint=False, order='nearest-neighbor')
        # 
        #if icol == 4 and irow == 2:
        #    output_hdu = fits.PrimaryHDU(data = new_image, header = reproject_header)
        #    output_file = 'debug_new_image.fits'
        #    if os.path.isfile(output_file):
        #        shutil.move(output_file, output_file+'.backup')
        #    output_hdu.writeto(output_file)
        #    print('Output to "%s"'%(output_file))
        #    # 
        #    output_hdu = image_hdu
        #    output_file = 'debug_input_image.fits'
        #    if os.path.isfile(output_file):
        #        shutil.move(output_file, output_file+'.backup')
        #    output_hdu.writeto(output_file)
        #    print('Output to "%s"'%(output_file))
        
        
        # note that the output and the input images have the same pixsc, 
        # so we do interpolation only for correcting for fractional pixel shift.
        # 
        # first set corner pixels, for the left/right-most regions, we use the image boundary.
        if icol == 0:
            corner00_ra = list_of_image_corners_coords[i][0].ra.deg
        else:
            corner00_ra = list_of_region_corners_coords[i][0].ra.deg
        if icol == ncol-1:
            corner11_ra = list_of_image_corners_coords[i][1].ra.deg
        else:
            corner11_ra = list_of_region_corners_coords[i][1].ra.deg
        if irow == 0:
            corner00_dec = list_of_image_corners_coords[i][0].dec.deg
        else:
            corner00_dec = list_of_region_corners_coords[i][0].dec.deg
        if irow == nrow-1:
            corner11_dec = list_of_image_corners_coords[i][1].dec.deg
        else:
            corner11_dec = list_of_region_corners_coords[i][1].dec.deg
        # 
        # convert ra dec to pixel
        out_corners_x, out_corners_y = output_wcs.wcs_world2pix([corner00_ra, corner11_ra], [corner00_dec, corner11_dec], 0)
        this_corners_x, this_corners_y = this_wcs.wcs_world2pix([corner00_ra, corner11_ra], [corner00_dec, corner11_dec], 0)
        out_corners_x = np.round(out_corners_x, 3).tolist()
        out_corners_y = np.round(out_corners_y, 3).tolist()
        this_corners_x = np.round(this_corners_x, 3).tolist()
        this_corners_y = np.round(this_corners_y, 3).tolist()
        print('out_corners_x, out_corners_y = %s, %s'%(out_corners_x, out_corners_y))
        print('this_corners_x, this_corners_y = %s, %s'%(this_corners_x, this_corners_y))
        ox0 = int(np.round(out_corners_x[0]))
        ox1 = int(np.floor(out_corners_x[1]))
        oy0 = int(np.round(out_corners_y[0]))
        oy1 = int(np.floor(out_corners_y[1]))
        ix0 = int(np.round(this_corners_x[0]))
        ix1 = int(np.floor(this_corners_x[1]))
        iy0 = int(np.round(this_corners_y[0]))
        iy1 = int(np.floor(this_corners_y[1]))
        print('ox0, ox1, oy0, oy1 = %s, %s, %s, %s'%(ox0, ox1, oy0, oy1))
        print('ix0, ix1, iy0, iy1 = %s, %s, %s, %s'%(ix0, ix1, iy0, iy1))
        nx = ox1-ox0+1
        ny = oy1-oy0+1
        print('nx = %s, ix1-ix0+1 = %s, ox1-ox0+1 = %s'%(nx, ix1-ix0+1, ox1-ox0+1))
        print('ny = %s, iy1-iy0+1 = %s, oy1-oy0+1 = %s'%(ny, iy1-iy0+1, oy1-oy0+1))
        # 
        # check image boundary
        tuned = False
        if ox0 < 0:
            dx0 = -ox0
            ox0 = ox0 + dx0
            ix0 = ix0 + dx0
            tuned = True
        if oy0 < 0:
            dy0 = -oy0
            oy0 = oy0 + dy0
            iy0 = iy0 + dy0
            tuned = True
        if ox1 > out_nx-1:
            dx1 = ox1-out_nx+1
            ox1 = ox1 - dx1
            ix1 = ix1 - dx1
            tuned = True
        if oy1 > out_ny-1:
            dy1 = oy1-out_ny+1
            oy1 = oy1 - dy1
            iy1 = iy1 - dy1
            tuned = True
        if ox0 > out_nx-1:
            continue
        if oy0 > out_ny-1:
            continue
        # 
        nx = ox1-ox0+1
        ny = oy1-oy0+1
        # 
        if nx != ix1-ix0+1:
            ix1 = ix0+nx-1
            tuned = True
        if ny != iy1-iy0+1:
            iy1 = iy0+ny-1
            tuned = True
        # 
        if ix1 > this_nx-1:
            dx1 = ix1-this_nx+1
            ix1 = ix1 - dx1
            ox1 = ox1 - dx1
            nx = nx - dx1
            tuned = True
        if iy1 > this_ny-1:
            dy1 = iy1-this_ny+1
            iy1 = iy1 - dy1
            oy1 = oy1 - dy1
            ny = ny - dy1
            tuned = True
        # 
        if tuned:
            print('ox0, ox1, oy0, oy1 = %s, %s, %s, %s (tuned)'%(ox0, ox1, oy0, oy1))
            print('ix0, ix1, iy0, iy1 = %s, %s, %s, %s (tuned)'%(ix0, ix1, iy0, iy1))
            print('nx = %s, ix1-ix0+1 = %s, ox1-ox0+1 = %s (tuned)'%(nx, ix1-ix0+1, ox1-ox0+1))
            print('ny = %s, iy1-iy0+1 = %s, oy1-oy0+1 = %s (tuned)'%(ny, iy1-iy0+1, oy1-oy0+1))
        # 
        this_gridy, this_gridx = np.mgrid[0:this_ny, 0:this_nx] # this full image grid
        tpixcoords = np.column_stack([this_gridx.flatten(), this_gridy.flatten()])
        fgridy, fgridx = np.mgrid[iy0:iy1+1, ix0:ix1+1] # dimension nx x ny for output
        fgridy = fgridy.astype(float) + (this_corners_y[0]-iy0) # fractional pixel coordinates within rectangle corners, but with dimension of nx x ny
        fgridx = fgridx.astype(float) + (this_corners_x[0]-ix0) # fractional pixel coordinates within rectangle corners, but with dimension of nx x ny
        fpixcoords = np.column_stack([fgridx.flatten(), fgridy.flatten()])
        # 
        # run griddata channel by channel
        for j in range(nchan):
            
            new_image_intpix = image_data[j, iy0:iy1+1, ix0:ix1+1]
            new_pb_intpix = pb_data[j, iy0:iy1+1, ix0:ix1+1]
            
            tdataarray = image_data[j, :, :].ravel()
            tdatamask = ~np.isnan(tdataarray)
            if np.count_nonzero(tdatamask) > 0:
                print('scipy.interpolate.griddata image_data[%d, %d:%d, %d:%d] (%d valid points) to fractioal pixel coordinates fgridx %.3f-%.3f fgridy %.3f-%.3f (%d points)'%(\
                    j, 0, this_ny, 0, this_nx, np.count_nonzero(tdatamask), 
                    fpixcoords[0][0], fpixcoords[-1][0], fpixcoords[0][1], fpixcoords[-1][1], nx*ny))
                odataarray = griddata(tpixcoords[tdatamask], \
                                      tdataarray[tdatamask], \
                                      fpixcoords, \
                                      method = 'cubic', \
                                      fill_value = np.nan ) # 2D cubic
            else:
                continue
            new_image = odataarray.reshape(new_image_intpix.shape)
            if np.count_nonzero(np.isnan(new_image_intpix))>0:
                new_image[np.isnan(new_image_intpix)] = np.nan
            
            tdataarray = pb_data[j, :, :].ravel()
            tdatamask = ~np.isnan(tdataarray)
            if np.count_nonzero(tdatamask) > 0:
                print('scipy.interpolate.griddata pb_data[%d, %d:%d, %d:%d] (%d valid points) to fractioal pixel coordinates fgridx %.3f-%.3f fgridy %.3f-%.3f (%d points)'%(\
                    j, 0, this_ny, 0, this_nx, np.count_nonzero(tdatamask), 
                    fpixcoords[0][0], fpixcoords[-1][0], fpixcoords[0][1], fpixcoords[-1][1], nx*ny))
                odataarray = griddata(tpixcoords[tdatamask], \
                                      tdataarray[tdatamask], \
                                      fpixcoords, \
                                      method = 'cubic', \
                                      fill_value = np.nan ) # 2D cubic
            else:
                continue
            new_pb = odataarray.reshape(new_pb_intpix.shape)
            if np.count_nonzero(np.isnan(new_pb_intpix))>0:
                new_pb[np.isnan(new_pb_intpix)] = np.nan
            
            print('channel: %s/%s, new_image.shape: %s, new_pb.shape: %s'%(j+1, nchan, new_image.shape, new_pb.shape))
            
            mask = np.logical_and(~np.isnan(new_image), ~np.isnan(new_pb))
            output_image_slice = output_image[j, oy0:oy1+1, ox0:ox1+1]
            output_weight_slice = output_weight[j, oy0:oy1+1, ox0:ox1+1]
            output_coverage_slice = output_coverage[j, oy0:oy1+1, ox0:ox1+1]
            output_image_slice[mask] += new_image[mask] * (new_pb[mask]**2)
            output_weight_slice[mask] += (new_pb[mask]**2)
            output_coverage_slice[mask] += 1
            output_image[j, oy0:oy1+1, ox0:ox1+1] += output_image_slice[:, :]
            output_weight[j, oy0:oy1+1, ox0:ox1+1] += output_weight_slice[:, :]
            output_coverage[j, oy0:oy1+1, ox0:ox1+1] += output_coverage_slice[:, :]
    
    # weight by pb**2
    mask = (output_coverage>0)
    output_image[mask] = output_image[mask] / output_weight[mask]
    output_image[~mask] = np.nan
    
    # put back OBSGEO
    #output_header = push_in_obsgeo(output_header, geox, geoy, geoz)
    
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








