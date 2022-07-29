# %%
import os
import numpy as np

from astropy import units as u
from astropy.nddata import CCDData
from astropy.io import fits
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.convolution import Gaussian2DKernel, convolve

from astroquery.astrometry_net import AstrometryNet
import astroalign as aa
import ccdproc
from ccdproc import Combiner

ast = AstrometryNet()


def avg_cube(directory, cubename, outname):
    """
    Create 2D fits img from sigma-clipped mean of 3D cube.
    """
    cube = CCDData.read(os.path.join(directory, cubename), unit=u.dimensionless_unscaled)
    # create ccdprocs combiner object
    combiner = Combiner(cube)
    # average combine with sigma clipping
    combiner.sigma_clipping(low_thresh=3, high_thresh=3)
    combined_average = combiner.average_combine()
    combined_average = np.asarray(combined_average)
    combined_average = combined_average.astype('float32')
    # write to 2D fits imgage with name=outname, header=original cube header
    fits.writeto(os.path.join(directory, outname), 
                combined_average, header=cube.header, overwrite=True)
    return


def calibrate_cube(directory, cubename, cali_dir, biasname, darkname, darkexptime, flatname):
    """
    Reduce 3D cube by master dark, bias, flat; output to subdir 'reduced_cubes'.
    """
    # load data cube and master calibration frames
    cube = CCDData.read(os.path.join(directory, cubename), unit=u.dimensionless_unscaled)
    master_bias = CCDData.read(os.path.join(cali_dir, biasname), unit=u.dimensionless_unscaled)
    master_dark = CCDData.read(os.path.join(cali_dir, darkname), unit=u.dimensionless_unscaled)
    master_flat = CCDData.read(os.path.join(cali_dir, flatname), unit=u.dimensionless_unscaled)

    # subtract master bias from master dark
    master_dark = ccdproc.subtract_bias(master_dark, master_bias)

    # subtract bias, dark, divide by flat
    reduced_cube = []
    nframes = int(cube.header['NAXIS3'])
    exptime = float(cube.header['EXPTIME'])
    for n in range(nframes):
        img = cube[n]
        bias_subtracted = ccdproc.subtract_bias(img, master_bias)
        dark_subtracted = ccdproc.subtract_dark(bias_subtracted, master_dark, 
                                                data_exposure=exptime*u.s, dark_exposure=darkexptime*u.s)
        reduced_img = ccdproc.flat_correct(dark_subtracted, master_flat)
        reduced_cube.append(reduced_img)

    # create 'reduced_cubes' subdirectory
    if not os.path.exists(os.path.join(directory, 'reduced_cubes')):
        print('Created /reduced_cubes subfolder')
        os.mkdir(os.path.join(directory, 'reduced_cubes'))
    
    # write to 3D cube in 'reduced_cubes', name=cubename_reduced, header=original cube header
    reduced_cube = np.asarray(reduced_cube)
    reduced_cube = reduced_cube.astype('float32')
    fits.writeto(os.path.join(directory,'reduced_cubes',cubename[:-5]+'_reduced.fits'), 
                reduced_cube, header=cube.header, overwrite=True)
    return


def solve_img(directory, imgname):
    """
    Submit to astrometry.net and solve image, return a WCS header
    """
    wcs_header = ast.solve_from_image(os.path.join(directory, imgname), force_image_upload=True)
    return wcs_header


def solve_and_align(directory, cubename):
    """
    Solve the first frame of a 3D cube, align other frames to it.
    Output an aligned WCS cube to subdir 'aligned_cubes'. 
    Note the wait time for astrometry.net varies greatly with time of day.
    """
    # read from 3D cube and obtain its 1st frame for solving and as reference frame
    hdulist = fits.open(os.path.join(directory, cubename))
    hdu = hdulist[0]
    ref_img = hdu.data[0]
    cube = np.float32(hdu.data)
    # IMPORTANT: np.float32 is to reset endian-ness for astroalign, do not modify
    # FITS from fits.writeto has a different endian-ness that is incompatable with astroalign
    # otherwise will raises ValueError: Big-endian buffer not supported on little-endian compiler
    # This issue is introduced in new versions of scikit-image that astroalign calls on

    # obtain essential keys from cube header
    nframes = int(hdu.header['NAXIS3'])
    assert nframes == len(cube)
    exptime = float(hdu.header['EXPTIME'])
    dateobs = hdu.header['DATE-OBS']

    # write 1st frame of cube to a temporary 2D fits and solve it
    fits.writeto(os.path.join(directory, 'solve_temporary.fits'), ref_img, header=hdu.header, overwrite=True)
    wcs_header = solve_img(directory, 'solve_temporary.fits')
    # os.remove(os.path.join(directory, 'solve_temporary.fits'))

    return wcs_header


            
# %%
