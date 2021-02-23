#!/usr/bin/env python
# 
# This needs to be run in CASA
# 
# CASA modules/functions used:
#     tb, casalog, mstransform, inp, saveinputs, exportfits
# 
# Example:
#     import dzliu_linear_mosaic; reload(dzliu_linear_mosaic); from a_dzliu_code_level_4_clean import dzliu_clean; dzliu_clean()
# 
# For old CASA 4.7.2
#     pip-2.7 install --target=~/Applications/CASA-472.app/Contents/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages astropy
# 
from __future__ import print_function
import os, sys, re, json, copy, timeit, time, datetime, shutil
import numpy as np
if '__casac__' in globals() or 'start_casa.py' in sys.argv[0]:
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
    if 'casalog' in globals():
        casalog.post(message, 'INFO')




# 
# def velo2freq
# 
def velo2freq(input_fits_header):
    # 
    if input_fits_header['NAXIS'] >= 3:
        if input_fits_header['CTYPE3'].strip().upper() == 'VRAD':
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
    return input_fits_header




# 
# my_fft2
# 
def my_fft2(data2D):
    naxis2, naxis1 = data2D.shape
    xc = (naxis1-1)/2.0
    yc = (naxis2-1)/2.0
    fft_plane_size = 2**math.ceil(math.log2(max(naxis1,naxis2))) + 1
    ygrid, xgrid = np.mgrid[0:naxis2, 0:naxis1]
    vgrid, ugrid = np.mgrid[0:fft_plane_size, 0:fft_plane_size]
    xgrid = xgrid - xc
    ygrid = ygrid - yc
    ugrid = ugrid - (fft_plane_size-1)/2.0
    vgrid = vgrid - (fft_plane_size-1)/2.0
    data2Dflat = data2D
    print('Calculating my_fft2 for %dx%d data points'%(fft_plane_size, fft_plane_size))
    print('data2D[120,128]', data2D[120,128], (xgrid[120,128], ygrid[120,128]))
    print('data2D[120,129]', data2D[120,129], (xgrid[120,129], ygrid[120,129]))
    print(                  [(ugrid[j,i],vgrid[j,i])
                                for j in range(120,125) \
                                    for i in [128] \
                            ])
    print(                  [ np.sum(data2Dflat * np.exp(-2.0*np.pi*1j*(ugrid[j,i]*xgrid + vgrid[j,i]*ygrid)/100.0)) \
                                for j in range(120,125) \
                                    for i in [128] \
                            ])
    sys.exit()
    fftdiagram2D = np.array([ np.sum(data2Dflat * np.exp(-2.0*np.pi*1j*(ugrid[j,i]*xgrid + vgrid[j,i]*ygrid))) \
                                for j in range(fft_plane_size) \
                                    for i in range(fft_plane_size) \
                            ]).reshape([fft_plane_size, fft_plane_size])
    return fftdiagram2D
    




# 
# def convolve_fits_cube_data
# 
def convolve_fits_cube_data(input_fits_cube_data, input_pixscale, input_BMAJ_BMIN_BPA, output_BMAJ_BMIN_BPA):
    # 
    # load module
    from astropy.convolution import convolve, Gaussian2DKernel
    from astropy.modeling.models import Gaussian2D
    #import numpy.fft.fft2 as fft2
    #import numpy.fft.ifft2 as ifft2
    #import numpy.fft.fftshift as fftshift
    #import numpy.fft.ifftshift as ifftshift
    import scipy.fftpack
    from scipy.fftpack import fft2
    from scipy.fftpack import ifft2
    from scipy.fftpack import fftshift
    from scipy.fftpack import ifftshift
    # 
    # construct convolution kernel which is dividing iFFT(FFT(output_Beam)/FFT(input_Beam))
    # 
    # check user input
    if type(input_fits_cube_data) is not np.ndarray:
        raise ValueError('Error! Please input a numpy.ndarray as the input_fits_cube_data parameter of convolve_fits_cube_data()! type(input_fits_cube_data) = %s'%(type(input_fits_cube_data)))
    if np.isscalar(input_BMAJ_BMIN_BPA) or len(input_BMAJ_BMIN_BPA) != 3:
        raise ValueError('Error! Please input a 3-element list or tuple as the input_BMAJ_BMIN_BPA parameter of convolve_fits_cube_data()! type(input_BMAJ_BMIN_BPA) = %s'%(type(input_BMAJ_BMIN_BPA)))
    if np.isscalar(output_BMAJ_BMIN_BPA) or len(output_BMAJ_BMIN_BPA) != 3:
        raise ValueError('Error! Please input a 3-element list or tuple as the output_BMAJ_BMIN_BPA parameter of convolve_fits_cube_data()! type(output_BMAJ_BMIN_BPA) = %s'%(type(output_BMAJ_BMIN_BPA)))
    if not np.isscalar(input_pixscale):
        input_pixscale = np.abs(input_pixscale[-1])
    # 
    print('convolve_fits_cube_data: input_BMAJ_BMIN_BPA = %s'%(input_BMAJ_BMIN_BPA))
    print('convolve_fits_cube_data: output_BMAJ_BMIN_BPA = %s'%(output_BMAJ_BMIN_BPA))
    # 
    if len(input_fits_cube_data.shape) == 3:
        naxis = 3
        naxis1 = input_fits_cube_data.shape[-1] # x
        naxis2 = input_fits_cube_data.shape[-2] # y
        nchan = input_fits_cube_data.shape[0]
    elif len(input_fits_cube_data.shape) == 2:
        naxis = 2
        naxis1 = input_fits_cube_data.shape[-1] # x
        naxis2 = input_fits_cube_data.shape[-2] # y
        nchan = 1
        input_fits_cube_data = input_fits_cube_data[np.newaxis, :, :]
    else:
        raise ValueError('Error! The input_fits_cube_data of convolve_fits_cube_data() should be either 2D or 3D!')
    # 
    # resample input fits cube
    #fft_plane_size = 2**math.ceil(math.log2(max(naxis1,naxis2))) + 1
    #print('fft_plane_size = %s'%(fft_plane_size))
    #if fft_plane_size < 512:
    #    fft_plane_size = 512
    # 
    # generate two Gaussian2D
    fwhm2stddev = 1.0/(2.0*np.sqrt(2.0*np.log(2.0)))
    print('fwhm2stddev = %s'%(fwhm2stddev))
    input_beam_Gaussian2D = Gaussian2D(1.0, \
                    naxis1/2.0, \
                    naxis2/2.0, \
                    input_BMAJ_BMIN_BPA[0]/input_pixscale*fwhm2stddev, \
                    input_BMAJ_BMIN_BPA[1]/input_pixscale*fwhm2stddev, \
                    np.deg2rad(input_BMAJ_BMIN_BPA[2]+90.0))
    output_beam_Gaussian2D = Gaussian2D(1.0, \
                    naxis1/2.0, \
                    naxis2/2.0, \
                    output_BMAJ_BMIN_BPA[0]/input_pixscale*fwhm2stddev, \
                    output_BMAJ_BMIN_BPA[1]/input_pixscale*fwhm2stddev, \
                    np.deg2rad(output_BMAJ_BMIN_BPA[2]+90.0))
    # 
    # mesh
    ygrid, xgrid = np.mgrid[0:naxis2, 0:naxis1]
    input_beam_Gaussian2D_image = input_beam_Gaussian2D(xgrid, ygrid)
    output_beam_Gaussian2D_image = output_beam_Gaussian2D(xgrid, ygrid)
    set_dump_debug_image = True
    if set_dump_debug_image:
        print('Dumping \'dump_input_beam_Gaussian2D_image.fits\'')
        hdu = fits.PrimaryHDU(input_beam_Gaussian2D_image)
        hdu.writeto('dump_input_beam_Gaussian2D_image.fits', overwrite=True)
        # 
        print('Dumping \'dump_output_beam_Gaussian2D_image.fits\'')
        hdu = fits.PrimaryHDU(output_beam_Gaussian2D_image)
        hdu.writeto('dump_output_beam_Gaussian2D_image.fits', overwrite=True)
    # 
    # prepare kernel
    output_fits_cube_data = np.full([nchan, naxis2, naxis1], np.nan)
    #input_beam_fft_plane_diagram = fftshift(fft2(input_beam_Gaussian2D_image, )) # [fft_plane_size, fft_plane_size]
    #output_beam_fft_plane_diagram = fftshift(fft2(output_beam_Gaussian2D_image, )) # [fft_plane_size, fft_plane_size]
    input_beam_fft_plane_diagram = my_fft2(input_beam_Gaussian2D_image)
    # 
    if set_dump_debug_image:
        print('Dumping \'dump_input_beam_fft_plane_diagram.fits\'')
        hdu = fits.PrimaryHDU(np.absolute(input_beam_fft_plane_diagram))
        hdu.writeto('dump_input_beam_fft_plane_diagram.fits', overwrite=True)
        hdu = fits.PrimaryHDU(input_beam_fft_plane_diagram.real)
        hdu.writeto('dump_input_beam_fft_plane_diagram.real.fits', overwrite=True)
        hdu = fits.PrimaryHDU(input_beam_fft_plane_diagram.imag)
        hdu.writeto('dump_input_beam_fft_plane_diagram.imag.fits', overwrite=True)
    # 
    output_beam_fft_plane_diagram = my_fft2(output_beam_Gaussian2D_image)
    #print(np.argwhere(np.isclose(np.absolute(input_beam_fft_plane_diagram), 0.0, atol=1e-60, rtol=1e-60)))
    # 
    if set_dump_debug_image:
        print('Dumping \'dump_output_beam_fft_plane_diagram.fits\'')
        hdu = fits.PrimaryHDU(np.absolute(output_beam_fft_plane_diagram))
        hdu.writeto('dump_output_beam_fft_plane_diagram.fits', overwrite=True)
        hdu = fits.PrimaryHDU(output_beam_fft_plane_diagram.real)
        hdu.writeto('dump_output_beam_fft_plane_diagram.real.fits', overwrite=True)
        hdu = fits.PrimaryHDU(output_beam_fft_plane_diagram.imag)
        hdu.writeto('dump_output_beam_fft_plane_diagram.imag.fits', overwrite=True)
    # 
    convolution_kernel_fft_plane_diagram = np.divide(output_beam_fft_plane_diagram, input_beam_fft_plane_diagram)
    # 
    if set_dump_debug_image:
        print('Dumping \'dump_convolution_kernel_fft_plane_diagram.fits\'')
        hdu = fits.PrimaryHDU(np.absolute(convolution_kernel_fft_plane_diagram))
        hdu.writeto('dump_convolution_kernel_fft_plane_diagram.fits', overwrite=True)
    # 
    # loop 3D
    print('Looping %d channels'%(nchan))
    for i in range(nchan):
        input_image_data = input_fits_cube_data[i, :, :]
        #input_image_fft_plane_diagram = fftshift(fft2(input_image_data, )) # [fft_plane_size, fft_plane_size]
        input_image_fft_plane_diagram = my_fft2(input_image_data)
        output_image_fft_plane_diagram = input_image_fft_plane_diagram * convolution_kernel_fft_plane_diagram
        output_fits_cube_data[i, :, :] = np.absolute(ifft2(ifftshift(output_image_fft_plane_diagram), [naxis2, naxis1] ) )
        if i == nchan-1 and set_dump_debug_image:
            # 
            print('Dumping \'dump_output_image_fft_plane_diagram.fits\'')
            hdu = fits.PrimaryHDU(np.absolute(output_image_fft_plane_diagram))
            hdu.writeto('dump_output_image_fft_plane_diagram.fits', overwrite=True)
    # 
    if naxis == 2:
        output_fits_cube_data = output_fits_cube_data[0, :, :]
    # 
    if set_dump_debug_image:
        print('Dumping \'dump_output_fits_cube_data.fits\'')
        hdu = fits.PrimaryHDU(output_fits_cube_data)
        hdu.writeto('dump_output_fits_cube_data.fits', overwrite=True)
    # 
    raise NotImplementedError('TODO: complex divided by zero error')
    # 
    return output_fits_cube_data


def test_convolve_fits_cube_data():
    hdul = fits.open('ngc3627_1_7m_co21.fits')
    input_fits_cube_data = hdul[0].data
    input_fits_cube_header = hdul[0].header
    input_fits_cube_wcs = WCS(input_fits_cube_header, naxis=2)
    input_pixscale = proj_plane_pixel_scales(input_fits_cube_wcs)
    input_BMAJ_BMIN_BPA = [input_fits_cube_header['BMAJ'], input_fits_cube_header['BMIN'], input_fits_cube_header['BPA']]
    output_BMAJ_BMIN_BPA = np.array([15.0, 15.0, 0.0])/3600.0
    convolve_fits_cube_data(input_fits_cube_data, input_pixscale, input_BMAJ_BMIN_BPA, output_BMAJ_BMIN_BPA)








# 
# def project_fits_cube_data
# 
def project_fits_cube_data(input_fits_cube_data, input_fits_header, template_fits_header):
    # 
    # We will project the input_fits_cube to the pixel grid of the template_fits_cube_wcs
    # by matching the World Coordinate System (WCS). 
    # 
    # We can process both fits cubes and images. 
    # 
    # No CASA module required. 
    # 
    import warnings
    from astropy.utils.exceptions import AstropyWarning
    warnings.simplefilter('ignore', category=AstropyWarning)
    from astropy.io import fits
    from astropy import units as u
    from astropy import wcs
    from astropy.wcs import WCS
    from astropy.wcs.utils import proj_plane_pixel_scales
    from astropy.coordinates import SkyCoord, FK5
    from scipy.interpolate import griddata
    # 
    # Check input
    #if type(input_fits_cube_wcs) is not WCS:
    #    raise ValueError('Error! The input input_fits_cube_wcs is not a astropy.wcs.WCS!')
    #if type(template_fits_cube_wcs) is not WCS:
    #    raise ValueError('Error! The input template_fits_cube_wcs is not a astropy.wcs.WCS!')
    #if len(list(input_fits_cube_data.shape)) != input_fits_cube_wcs.naxis:
    #    raise ValueError('Error! The input input_fits_cube_data has a shape of %s which is inconsistent with the input WCS with NAXIS %d!'%(list(input_fits_cube_data.shape), input_fits_cube_wcs.naxis))
    if type(input_fits_header) is not fits.Header:
        raise ValueError('Error! The input input_fits_header is not a astropy.io.fits.Header!')
    if type(template_fits_header) is not fits.Header:
        raise ValueError('Error! The input template_fits_header is not a astropy.io.fits.Header!')
    if len(list(input_fits_cube_data.shape)) != input_fits_header['NAXIS']:
        raise ValueError('Error! The input input_fits_cube_data has a shape of %s which is inconsistent with the input fits header with NAXIS %d!'%(list(input_fits_cube_data.shape), input_fits_header['NAXIS']))
    # 
    #template_fits_header = template_fits_cube_wcs.to_header()
    #input_fits_header = input_fits_cube_wcs.to_header()
    idata = input_fits_cube_data
    # 
    # Make sure the input is a 3D cube or at least a 2D image
    if int(input_fits_header['NAXIS']) < 2:
        raise Exception('Error! The input fits header does not have more than 2 dimensions!')
    if int(template_fits_header['NAXIS']) < 2:
        raise Exception('Error! The template fits header does not have more than 2 dimensions!')
    # 
    # Do velocity-to-frequency conversion before checking header consistency
    if input_fits_header['CTYPE3'].strip().upper() == 'VRAD' and template_fits_header['CTYPE3'].strip().upper() == 'FREQ':
        input_fits_header = velo2freq(input_fits_header)
    # 
    # Take the minimum dimension of the input and the template fits dimension.
    naxis = min(int(input_fits_header['NAXIS']), int(template_fits_header['NAXIS']))
    # 
    # Check header consistency
    for i in range(1, naxis+1):
        if input_fits_header['CTYPE%d'%(i)].strip().upper() != template_fits_header['CTYPE%d'%(i)].strip().upper():
            raise Exception('Error! The input fits cube CTYPE%d is %s but the template CTYPE%d is %s!'%(\
                                i, input_fits_header['CTYPE%d'%(i)].strip(), \
                                i, template_fits_header['CTYPE%d'%(i)].strip() ) )
    # 
    # Store dimension arrays 'NAXISi'
    inaxis = np.array([int(input_fits_header['NAXIS%d'%(i)]) for i in range(1, input_fits_header['NAXIS']+1)]) # [nx, ny, nchan, ....], it is inverted to the Python array dimension order
    onaxis = np.array([int(template_fits_header['NAXIS%d'%(i)]) for i in range(1, template_fits_header['NAXIS']+1)]) # [nx, ny, nchan, ....], it is inverted to the Python array dimension order
    inaxis_str = 'x'.join(inaxis[0:naxis].astype(str))
    onaxis_str = 'x'.join(onaxis[0:naxis].astype(str)) # output naxis, same as template
    idatashape = copy.copy(inaxis[0:naxis][::-1])
    odatashape = copy.copy(onaxis[0:naxis][::-1])
    # 
    # If input fits data have extra higher dimensions than the template data, we will condense the additional dimensions into one extra dimension, and later we will loop them over and do the interpolation slice-by-slice.
    idataslice = 0
    if len(inaxis) > naxis:
        print2('Warning! The input fits cube has %d dimensions while the template has only %d dimensions! We will loop the extra higher dimensions of the input data slice-by-slice and project to template pixel grid.'%(input_fits_header['NAXIS'], template_fits_header['NAXIS']))
        idataslice = np.product(inaxis[naxis:])
        idatashape = np.concatenate([[idataslice], idatashape])
        idata.shape = idatashape # this will reshape the idata array.
        odatashape = np.concatenate([inaxis[naxis:][::-1], odatashape])
        # here we reshape the input data to shrink all higher dimensions into one extra dimension, 
        # e.g., from (m,n,a,b,c) to (x,a,b,c), when the template data shape is (a,b,c), and x = m*n. 
    # 
    # Otherwise if the template fits data have extra higher dimensions, we will only use the common dimensions, while pad extra dimensions with np.newaxis
    elif len(onaxis) > naxis:
        print2('Warning! The input fits cube has %d dimensions while the template has %d dimensions! We will ignore the extra higher dimensions in the template for computation and simply reshape the output data to the template dimensions.'%(input_fits_header['NAXIS'], template_fits_header['NAXIS']))
        odatashape = []
        idataslice = -(len(onaxis) - naxis) # negative value means template has more dimension than the input fits cube.
        odatashape = np.concatenate([[1]*idataslice, odatashape]) # just fill extra dimensions with 1.
        output_fits_data_shape = odatashape
    # 
    # Otherwise the output fits data shape is the same as the template fits data shape
    else:
        output_fits_data_shape = copy.copy(onaxis[::-1])
    # 
    # Make pixel mgrid
    timestart = timeit.default_timer()
    set_do_3D_interpolation = False
    if naxis == 2:
        set_do_3D_interpolation = False
    elif naxis == 3:
        # 
        # Check channel number consistency
        if naxis == 3:
            if onaxis[2] != inaxis[2]:
                print2('Warning! The input fits cube and template fits cube have different channel number! Will do 3D griddata interpolation which is very slow!')
                set_do_3D_interpolation = True
    else:
        raise NotImplementedError('Error! The cube projection and interpolation have not been implemented for an NAXIS of %d!'%(naxis))
    # 
    # if set_do_3D_interpolation
    if set_do_3D_interpolation:
        # 
        # Get input fits WCS with naxis=naxis
        iwcs = WCS(input_fits_header, naxis=3)
        # 
        # Get template fits WCS with naxis=naxis
        twcs = WCS(template_fits_header, naxis=3)
        # 
        # input grid 3D
        print2('Generating pixel mgrid with %s pixels'%(inaxis_str))
        ichan, iy, ix = np.mgrid[0:inaxis[2], 0:inaxis[1], 0:inaxis[0]]
        ipixcoords = np.column_stack([ix.flatten(), iy.flatten(), ichan.flatten()])
        # 
        # output grid = template grid
        print2('Generating pixel mgrid with %s pixels'%(onaxis_str))
        tchan, ty, tx = np.mgrid[0:onaxis[2], 0:onaxis[1], 0:onaxis[0]]
        tpixcoords = np.column_stack([tx.flatten(), ty.flatten(), tchan.flatten()])
    else:
        # 
        # Get input fits WCS with naxis=naxis
        iwcs = WCS(input_fits_header, naxis=2)
        # 
        # Get template fits WCS with naxis=naxis
        twcs = WCS(template_fits_header, naxis=2)
        # 
        # input grid 2D
        print2('Generating pixel mgrid with %s pixels'%(inaxis_str[0:2]))
        iy, ix = np.mgrid[0:inaxis[1], 0:inaxis[0]]
        ipixcoords = np.column_stack([ix.flatten(), iy.flatten()])
        # 
        # output grid = template grid
        print2('Generating pixel mgrid with %s pixels'%(onaxis_str[0:2]))
        ty, tx = np.mgrid[0:onaxis[1], 0:onaxis[0]]
        tpixcoords = np.column_stack([tx.flatten(), ty.flatten()])
    # 
    # Timer stops
    timestop = timeit.default_timer()
    print2('Used %s seconds'%(timestop-timestart))
    # 
    # Timer starts
    timestart = timeit.default_timer()
    # 
    # Convert each pixel coordinate to skycoordinate for the template pixel grid which is also the output pixel grid.
    print2('Computing wcs_pix2world for %s pixels'%(onaxis_str))
    oskycoords = twcs.wcs_pix2world(tpixcoords, 0)
    print2('oskycoords.shape = %s'%(str(list(oskycoords.shape))))
    # 
    # Convert each pixel skycoordinate to the coordinate in the input mask cube, so that we can do interpolation. 
    print2('Computing wcs_world2pix for %s pixels'%(onaxis_str))
    opixcoords = iwcs.wcs_world2pix(oskycoords, 0)
    print2('opixcoords.shape = %s'%(str(list(opixcoords.shape))))
    # 
    # Timer stops
    timestop = timeit.default_timer()
    print2('Used %s seconds'%(timestop-timestart))
    # 
    # Loop each input data slice
    print2('Looping %d data slices...'%(max(1,idataslice)))
    timestart = timeit.default_timer()
    odata = []
    odataarray = None
    for i in range(max(1,idataslice)):
        # 
        # 
        if idataslice > 0:
            idataarray = idata[i]
        else:
            idataarray = idata
        # 
        # 
        if naxis == 2:
            idataarray = idataarray[np.newaxis, :, :]
            odataarray = np.full([1, onaxis[1], onaxis[0]], np.nan) # 1, ny, nx
        else:
            odataarray = np.full([onaxis[2], onaxis[1], onaxis[0]], np.nan) # nchan, ny, nx
        # 
        # 
        if set_do_3D_interpolation:
            # 
            # Do interpolation with scipy.interpolate.griddata in 3D
            print2('Interpolating griddata in 3D for data slice %d/%d ...'%(i+1, max(1,idataslice)))
            print2('idataarray.shape = %s'%(list(idataarray.shape)))
            print2('ipixcoords.shape = %s'%(list(ipixcoords.shape)))
            print2('opixcoords.shape = %s'%(list(opixcoords.shape)))
            # note that if the last dimension has a size of 1, griddata will fail. So we have to wrap around.
            if naxis == 3 and idata.shape[0] == 1:
                ipixcoords_2 = ipixcoords[:, 0:2]
                opixcoords_2 = opixcoords[:, 0:2]
                idataarray_2 = idataarray
                idatamask_2 = ~np.isnan(idataarray_2)
                if np.count_nonzero(idatamask_2) > 0:
                    odataarray = griddata(ipixcoords_2[idatamask_2], \
                                          idataarray_2[idatamask_2], \
                                          opixcoords_2, \
                                          method = 'cubic', \
                                          fill_value = np.nan ) # 2D cubic
            else:
                # real 3D interpolation, very slow....
                idatamask = ~np.isnan(idataarray)
                if np.count_nonzero(idatamask_2) > 0:
                    odataarray = griddata(ipixcoords[idatamask], \
                                          idataarray[idatamask], \
                                          opixcoords, \
                                          method = 'linear', \
                                          fill_value = np.nan ) # 3D linear
            # 
            # The interpolation is done with serialized arrays, so we reshape the output interpolated aray to 3D cube 
            odataarray = odataarray.reshape(onaxis[0:naxis][::-1]).astype(idata.dtype)
            # 
        else:
            # 
            # simple channel-by-channel 2D interpolation
            for j in range(idataarray.shape[0]):
                print2('Interpolating 2D image at channel %d/%d for data slice %d/%d'%(j+1, idataarray.shape[0], i+1, max(1,idataslice)))
                ipixcoords_2 = ipixcoords # list of [x,y] coordinate pairs
                opixcoords_2 = opixcoords # list of [x,y] coordinate pairs
                idataarray_2 = idataarray[j, :, :].flatten() # data at this slice at this slice at this channel
                idatamask_2 = ~np.isnan(idataarray_2)
                if np.count_nonzero(idatamask_2) > 0:
                    odataarray_2 = griddata(ipixcoords_2[idatamask_2], \
                                            idataarray_2[idatamask_2], \
                                            opixcoords_2, \
                                            method = 'cubic', \
                                            fill_value = np.nan ) # 2D linear
                                            # method = 'linear'
                    odataarray[j, :, :] = odataarray_2.reshape([onaxis[1], onaxis[0]]) # ny, nx
            # 
        # 
        # 
        if idataslice > 0:
            odata.append(odataarray)
        else:
            odata = odataarray
        # 
        # Timer stops
        timestop = timeit.default_timer()
        print2('Used %s seconds'%(timestop-timestart))
    # 
    odata = np.array(odata)
    # 
    print2('input_fits_data_shape = %s'%(idatashape.tolist()))
    print2('output_fits_data_shape = %s'%(odatashape.tolist()))
    odata.shape = odatashape
    # 
    return odata






# 
# def examine_overlap_pixels
# 
def examine_overlap_pixels(input_image_1, input_image_2):
    # 
    # Input two data cubes with the same 3rd dimension.
    # 
    import warnings
    from astropy.utils.exceptions import AstropyWarning
    warnings.simplefilter('ignore', category=AstropyWarning)
    from astropy.io import fits
    from astropy import units as u
    from astropy import wcs
    from astropy.wcs import WCS
    from astropy.wcs.utils import proj_plane_pixel_scales
    from astropy.coordinates import SkyCoord, FK5
    raise NotImplementedError('Sorry, dzliu_linear_mosaic.examine_overlap_pixels() not implemented! Use phangs_alma_QA_for_multipart_combined_cube.py instead!')
    




# 
# def dzliu_linear_mosaic
# 
def dzliu_linear_mosaic(input_image_name_list, output_fits_cube):
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
    from astropy.io import fits
    from astropy import units as u
    from astropy import wcs
    from astropy.wcs import WCS
    from astropy.wcs.utils import proj_plane_pixel_scales
    from astropy.coordinates import SkyCoord, FK5
    # 
    # check input
    if np.isscalar(input_image_name_list):
        input_image_name_list = [input_image_name_list]
    # 
    # find image fits file and pb file
    input_fits_cube_list = []
    input_fits_pb_list = []
    for i in range(len(input_image_name_list)):
        input_fits_cube_list.append(input_image_name_list[i])
        if re.match(r'.*\.image\.pbcor\.fits$', input_image_name_list[i], re.IGNORECASE):
            input_image_name_list[i] = re.sub(r'\.image\.pbcor\.fits$', r'', input_image_name_list[i], re.IGNORECASE)
            is_pbcorreced = True
        elif re.match(r'.*\.image\.fits$', input_image_name_list[i], re.IGNORECASE):
            input_image_name_list[i] = re.sub(r'\.image\.fits$', r'', input_image_name_list[i], re.IGNORECASE)
        elif re.match(r'.*\.pb\.fits$', input_image_name_list[i], re.IGNORECASE):
            input_image_name_list[i] = re.sub(r'\.pb\.fits$', r'', input_image_name_list[i], re.IGNORECASE)
        elif re.match(r'.*\.fits$', input_image_name_list[i], re.IGNORECASE):
            input_image_name_list[i] = re.sub(r'\.fits$', r'', input_image_name_list[i], re.IGNORECASE)
        elif re.match(r'.*\.image$', input_image_name_list[i], re.IGNORECASE):
            input_image_name_list[i] = re.sub(r'\.image$', r'', input_image_name_list[i], re.IGNORECASE)
        elif re.match(r'.*_pb\.fits$', input_image_name_list[i], re.IGNORECASE):
            input_image_name_list[i] = re.sub(r'_pb\.fits$', r'', input_image_name_list[i], re.IGNORECASE)
        else:
            raise ValueError('Error! The input fits image "%s" should ends with ".image.pbcor.fits", ".image.fits" or ".fits" or ".image"!')
        # check fits cube file existence
        if not os.path.isfile(input_fits_cube_list[i]):
            raise Exception('Error! The fits image "%s" was not found! input_image_name_list[%d]: %s'%(input_fits_cube_list[i], i, input_image_name_list[i]))
        # check fits pb file existence and set 'input_fits_pb_list'
        if os.path.isfile(input_image_name_list[i]+'_pb.fits'):
            input_fits_pb_list.append(input_image_name_list[i]+'_pb.fits')
        elif os.path.isfile(input_image_name_list[i]+'.pb.fits'):
            input_fits_pb_list.append(input_image_name_list[i]+'.pb.fits')
        else:
            raise Exception('Error! The fits pb "%s" was not found! input_image_name_list[%d]: %s'%(input_image_name_list[i]+'.pb.fits', i, input_image_name_list[i]))
    # 
    if len(input_fits_pb_list) != len(input_fits_cube_list):
        raise Exception('Error! The input fits cube and pb are inconsistent! input_fits_cube_list: %s; input_fits_pb_list: %s'%(input_fits_cube_list, input_fits_pb_list))
    # 
    ninput = len(input_fits_cube_list)
    # 
    # 
    # check output, remove suffix
    #if output_fits_cube is None:
    #    output_fits_cube = 'run_linear_mosaic_%s/linear_mosaic'%(os.path.basename(input_fits_cube_list[0]))
    # 
    if re.match(r'.*\.fits$', output_fits_cube, re.IGNORECASE):
        output_fits_cube = re.sub(r'\.fits$', r'', output_fits_cube, re.IGNORECASE)
    # 
    if output_fits_cube.find(os.sep)>=0:
        if not os.path.isdir(os.path.dirname(output_fits_cube)):
            os.makedirs(os.path.dirname(output_fits_cube))
        if not os.path.isdir(os.path.dirname(output_fits_cube)):
            raise Exception('Error! Could not create output directory "%s"!'%(os.path.dirname(output_fits_cube)))
    # 
    # read fits header and wcs info
    fits_header_list = []
    fits_wcs_2D_list = []
    fits_pixscale_list = []
    fits_dimension_list = []
    fits_corner_00_RA_list = []
    fits_corner_00_Dec_list = []
    fits_corner_11_RA_list = []
    fits_corner_11_Dec_list = []
    fits_corner_01_RA_list = []
    fits_corner_01_Dec_list = []
    fits_corner_10_RA_list = []
    fits_corner_10_Dec_list = []
    for i in range(ninput):
        if not re.match(r'.*\.fits$', input_fits_cube_list[i], re.IGNORECASE):
            raise Exception('Error! Please input fits files!')
        if not re.match(r'.*\.fits$', input_fits_pb_list[i], re.IGNORECASE):
            raise Exception('Error! Please input fits files!')
        fits_header = None
        with fits.open(input_fits_cube_list[i]) as hdulist:
            hdu = hdulist[0]
            fits_header = copy.copy(hdu.header)
        fits_wcs_2D = WCS(fits_header, naxis=2)
        fits_corner_coords = fits_wcs_2D.wcs_pix2world([[0.5, 0.5], 
                                                        [fits_header['NAXIS1']+0.5, fits_header['NAXIS2']+0.5],
                                                        [0.5, fits_header['NAXIS2']+0.5],
                                                        [fits_header['NAXIS1']+0.5, 0.5]], 
                                                       1) # lower-left, and upper-right.
                                                          # pixel corners are 0.5, 0.5! see -- https://docs.astropy.org/en/stable/_modules/astropy/wcs/wcs.html
        print2('fits_corner_coords[%d][00][RA,Dec] = %.10f, %.10f, image = %r'%(i, fits_corner_coords[0][0], fits_corner_coords[0][1], input_fits_cube_list[i]))
        print2('fits_corner_coords[%d][11][RA,Dec] = %.10f, %.10f, image = %r'%(i, fits_corner_coords[1][0], fits_corner_coords[1][1], input_fits_cube_list[i]))
        print2('fits_corner_coords[%d][01][RA,Dec] = %.10f, %.10f, image = %r'%(i, fits_corner_coords[2][0], fits_corner_coords[2][1], input_fits_cube_list[i]))
        print2('fits_corner_coords[%d][10][RA,Dec] = %.10f, %.10f, image = %r'%(i, fits_corner_coords[3][0], fits_corner_coords[3][1], input_fits_cube_list[i]))
        fits_header_list.append(copy.copy(fits_header))
        fits_wcs_2D_list.append(copy.copy(fits_wcs_2D))
        fits_pixscale_list.append(proj_plane_pixel_scales(fits_wcs_2D)[1]*3600.0) # arcsec
        fits_dimension_list.append(fits_header['NAXIS'])
        fits_corner_00_RA_list.append(fits_corner_coords[0][0]) # note that RA increases to the left, so this corner has the largest RA. 
        fits_corner_00_Dec_list.append(fits_corner_coords[0][1])
        fits_corner_11_RA_list.append(fits_corner_coords[1][0])
        fits_corner_11_Dec_list.append(fits_corner_coords[1][1])
        fits_corner_01_RA_list.append(fits_corner_coords[2][0]) # note that RA increases to the left, so this corner has the largest RA. 
        fits_corner_01_Dec_list.append(fits_corner_coords[2][1])
        fits_corner_10_RA_list.append(fits_corner_coords[3][0])
        fits_corner_10_Dec_list.append(fits_corner_coords[3][1])
    # 
    # make sure pixel scales are the same, otherwise <TODO>
    for i in range(ninput):
        if not np.isclose(fits_pixscale_list[i], fits_pixscale_list[0]):
            raise Exception('Error! Pixel scales are inconsistent! Values are %s for input fits cubes %s.'%(fits_pixscale_list, input_fits_cube_list))
    # 
    # make sure all are images or cubes, and if cubes all channels are consistent
    naxis = 0
    if np.min(fits_dimension_list) < 2:
        raise Exception('Error! Data dimensions are smaller than 2! Dimensions are %s for input fits cubes %s.'%(fits_dimension_list, input_fits_cube_list))
    elif np.min(fits_dimension_list) == 2:
        naxis = 2
        for i in range(ninput):
            if fits_dimension_list[i] != 2:
                raise Exception('Error! Data dimensions are inconsistent! Dimensions are %s for input fits cubes %s.'%(fits_dimension_list, input_fits_cube_list))
    else:
        naxis = 3
        nchan = fits_header_list[0]['NAXIS3']
        for i in range(ninput):
            if fits_header_list[i]['NAXIS3'] != nchan:
                raise Exception('Error! The third dimension of the data cubes are inconsistent! NAXIS3 are %s for input fits cubes %s.'%(str([t['NAXIS3'] for t in fits_header_list]), input_fits_cube_list))
    print2('naxis = %d'%(naxis))
    # 
    # compute output image sizes and prepare output fits header and wcs
    pixscale = np.max(np.abs(fits_pixscale_list))
    argmin_RA = np.argmin(fits_corner_11_RA_list).flatten()[0] # upper right
    argmin_Dec = np.argmin(fits_corner_00_Dec_list).flatten()[0]
    argmax_RA = np.argmax(fits_corner_01_RA_list).flatten()[0] # upper left
    argmax_Dec = np.argmax(fits_corner_11_Dec_list).flatten()[0]
    min_RA = fits_corner_11_RA_list[argmin_RA]
    min_Dec = fits_corner_00_Dec_list[argmin_Dec]
    max_RA = fits_corner_01_RA_list[argmax_RA]
    max_Dec = fits_corner_11_Dec_list[argmax_Dec]
    print2('argmin_RA = %s'%(argmin_RA))
    print2('argmin_Dec = %s'%(argmin_Dec))
    print2('argmax_RA = %s'%(argmax_RA))
    print2('argmax_Dec = %s'%(argmax_Dec))
    print2('min_RA = %s'%(min_RA))
    print2('min_Dec = %s'%(min_Dec))
    print2('max_RA = %s'%(max_RA))
    print2('max_Dec = %s'%(max_Dec))
    nx = (max_RA - min_RA) * np.cos(np.deg2rad(max_Dec)) * 3600.0 / pixscale
    ny = (max_Dec - min_Dec) * 3600.0 / pixscale
    nx = int(np.ceil(nx))
    ny = int(np.ceil(ny))
    if naxis == 2:
        output_slices = np.full([ninput, ny, nx], 0.0)
        output_weights = np.full([ninput, ny, nx], 0.0)
        output_data = np.full([ny, nx], np.nan)
        output_mean = np.full([ny, nx], 0.0)
        output_cov = np.full([ny, nx], 0)
    elif naxis == 3:
        output_slices = np.full([ninput, nchan, ny, nx], 0.0)
        output_weights = np.full([ninput, nchan, ny, nx], 0.0)
        output_data = np.full([nchan, ny, nx], np.nan)
        output_mean = np.full([nchan, ny, nx], 0.0)
        output_cov = np.full([nchan, ny, nx], 0)
    else:
        raise NotImplementedError('NAXIS %d not implemented!'%(naxis))
    # 
    # copy header keywords, recalculate the reference pixel coordinate using the left-most image (argmax_RA)
    # if the left-most image is not the bottom-most image, then also account for global offset in y from the 
    # bottom edge of the left-most image to the bottom edge of the bottom-most image. 
    output_hdu = fits.PrimaryHDU(data = output_data)
    output_header = copy.copy(output_hdu.header)
    if argmax_RA != argmin_Dec:
        corner_00_offset_in_Dec = (fits_corner_00_Dec_list[argmax_RA] - fits_corner_00_Dec_list[argmin_Dec]) # argmin_Dec defines the lower boudary of the output image, argmax_RA defines the left boudary. 
        corner_00_offset_in_y = corner_00_offset_in_Dec * 3600.0 / pixscale
        corner_00_offset_in_y = int(np.ceil(corner_00_offset_in_y))
        print2('corner_00_offset_in_Dec = %s'%(corner_00_offset_in_Dec))
        print2('corner_00_offset_in_y = %s'%(corner_00_offset_in_y))
    else:
        corner_00_offset_in_y = 0
    # 
    for i in range(1, naxis+1):
        for t in ['CTYPE', 'CUNIT', 'CRVAL', 'CRPIX', 'CDELT']:
            key = '%s%d'%(t, i)
            if key in fits_header_list[argmax_RA]:
                output_header[key] = fits_header_list[argmax_RA][key]
    # 
    output_header['CRPIX1'] = fits_header_list[argmax_RA]['CRPIX1']
    output_header['CRPIX2'] = fits_header_list[argmax_RA]['CRPIX2'] + corner_00_offset_in_y
    # 
    for key in ['EQUINOX', 'RADESYS', 'LONPOLE', 'LATPOLE', 'RESTFRQ', 'SPECSYS', 'ALTRVAL', 'ALTRPIX', 'VELREF', 
                'TELESCOP', 'INSTRUME', 'OBSERVER', 'DATE-OBS', 'TIMESYS', 'OBSRA', 'OBSDEC', 'OBSGEO-X', 'OBSGEO-Y', 'OBSGEO-Z', 
                'OBJECT', 'BMAJ', 'BMIN', 'BPA', 'BUNIT', 'BTYPE']:
        if key in fits_header_list[argmax_RA]:
            output_header[key] = fits_header_list[argmax_RA][key]
    # 
    output_header['DATE'] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + 'T' + time.strftime('%Z')
    output_header['ORIGIN'] = 'dzliu_linear_mosaic.dzliu_linear_mosaic()'
    # 
    #print(output_header)
    # 
    #for i in range(len()):
    #for t in ['RESTFRQ', 'SPECSYS', ]
    #for t in ['BMAJ', 'BMIN', 'BPA', ]
    # 
    # save to disk
    output_hdu = fits.PrimaryHDU(data = output_data, header = output_header)
    print2('writing cache blank fits file...')
    if os.path.isfile(output_fits_cube+'.cache.blank.fits'):
        os.remove(output_fits_cube+'.cache.blank.fits')
    output_hdu.writeto(output_fits_cube+'.cache.blank.fits')
    # 
    # Get WCS
    #output_wcs = WCS(output_header, naxis=naxis)
    # 
    # project each input data cube to the output data cube pixel coordinate
    #for i in range(ninput):
    #    project_fits_cube(input_fits_cube_list[i], 
    #                      output_fits_cube+'.cache.blank.fits', 
    #                      output_fits_cube+'.cache.projected.%d.fits'%(i), 
    #                      overwrite = True)
    # 
    # project and coadd
    for i in range(ninput):
        print2('-'*80)
        print2('Processing %d of %d input image'%(i+1, ninput))
        fits_image = None
        fits_pb = None
        with fits.open(input_fits_cube_list[i]) as hdulist:
            hdu = hdulist[0]
            fits_image = project_fits_cube_data(hdu.data, hdu.header, output_header)
        # 
        with fits.open(input_fits_pb_list[i]) as hdulist:
            hdu = hdulist[0]
            fits_pb = project_fits_cube_data(hdu.data, hdu.header, output_header)
        # 
        print2('fits_image.shape = %s'%(str(list(fits_image.shape))))
        print2('fits_pb.shape = %s'%(str(list(fits_pb.shape))))
        # 
        while len(fits_image.shape) > naxis:
            fits_image = fits_image[0]
        while len(fits_pb.shape) > naxis:
            fits_pb = fits_pb[0]
        # 
        print2('fits_image.shape = %s'%(str(list(fits_image.shape))))
        print2('fits_pb.shape = %s'%(str(list(fits_pb.shape))))
        # 
        mask = ~np.isnan(fits_image)
        output_cov[mask] += 1
        output_weights[i][mask] = (fits_pb[mask])**2
        output_slices[i][mask] = fits_image[mask]
        output_mean[mask] += fits_image[mask]
    # 
    # normalize weights
    mask_4D = (output_weights>0)
    mask_3D = (output_cov>0)
    output_weightsum_3D = np.sum(output_weights, axis=0)
    print('mask_4D.shape', mask_4D.shape)
    print('mask_3D.shape', mask_3D.shape, 'np.count_nonzero(mask_3D)', np.count_nonzero(mask_3D))
    print('output_weightsum_3D.shape', output_weightsum_3D.shape, 'np.count_nonzero(output_weightsum_3D)', np.count_nonzero(output_weightsum_3D))
    output_weights = output_weights / output_weightsum_3D
    # 
    # weighted mean
    output_data[mask_3D] = np.sum(output_slices * output_weights, axis=0)[mask_3D]
    # 
    # simple mean
    output_mean[mask_3D] = output_mean[mask_3D] / output_cov[mask_3D].astype(float)
    # 
    # check output directory 
    if os.sep in output_fits_cube:
        if not os.path.isdir(os.path.dirname(output_fits_cube)):
            os.makedirs(os.path.dirname(output_fits_cube))
    # 
    # Output the interpolated mask cube as a fits file
    print2('writing final coadded weighted-mean...')
    output_hdu = fits.PrimaryHDU(data = output_data, header = output_header)
    if os.path.isfile(output_fits_cube+'.coadd.fits'):
        shutil.move(output_fits_cube+'.coadd.fits', output_fits_cube+'.coadd.fits'+'.backup')
    output_hdu.writeto(output_fits_cube+'.coadd.fits')
    print2('Output to "%s"!'%(output_fits_cube+'.coadd.fits'))
    # 
    print2('writing final coadded simple-mean file...')
    output_hdu = fits.PrimaryHDU(data = output_mean, header = output_header)
    if os.path.isfile(output_fits_cube+'.mean.fits'):
        shutil.move(output_fits_cube+'.mean.fits', output_fits_cube+'.mean.fits'+'.backup')
    output_hdu.writeto(output_fits_cube+'.mean.fits')
    print2('Output to "%s"!'%(output_fits_cube+'.mean.fits'))
    # 
    print2('writing final coadded coverage file...')
    output_hdu = fits.PrimaryHDU(data = output_cov, header = output_header)
    if os.path.isfile(output_fits_cube+'.cov.fits'):
        shutil.move(output_fits_cube+'.cov.fits', output_fits_cube+'.cov.fits'+'.backup')
    output_hdu.writeto(output_fits_cube+'.cov.fits')
    print2('Output to "%s"!'%(output_fits_cube+'.cov.fits'))







############
#   main   #
############

dzliu_main_func_name = 'dzliu_linear_mosaic' # make sure this is the right main function in this script file

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




