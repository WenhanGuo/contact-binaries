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


def calibrate_flat(cali_dir, flatname, darkname, biasname):
    """
    Subtract master_dark and master_bias from master_flat.
    """
    # load master calibration frames
    master_bias = CCDData.read(os.path.join(cali_dir, biasname), unit=u.dimensionless_unscaled)
    master_dark = CCDData.read(os.path.join(cali_dir, darkname), unit=u.dimensionless_unscaled)
    master_flat = CCDData.read(os.path.join(cali_dir, flatname), unit=u.dimensionless_unscaled)

    # subtract master_dark and master_bias from master_flat
    # 10s are dummy variables
    bias_subtracted = ccdproc.subtract_bias(master_flat, master_bias)
    dark_subtracted = ccdproc.subtract_dark(bias_subtracted, master_dark, 
                                                data_exposure=10*u.s, dark_exposure=10*u.s)

    fits.writeto(os.path.join(cali_dir, flatname[:-5]+'_calibrated.fits'), 
                    dark_subtracted, header=master_flat.header, overwrite=True)
    


def calibrate_cube(directory, cubename, cali_dir, biasname, darkname, darkexptime, 
        flatname=None, mode='multicolor', Cflat=None, gflat=None, rflat=None, iflat=None, out_dir=None):
    """
    Reduce 3D cube by master dark, bias, flat; output to out_dir.
    If mode == 'multicolor', must provide CLEAR, g, r, i master flats;
    If mode == 'monocolor', must provide flatname.
    """
    assert mode in ['monocolor', 'multicolor']

    if mode == 'multicolor':
        assert Cflat != None
        assert gflat != None
        assert rflat != None
        assert iflat != None
    if mode == 'monocolor':
        assert flatname != None

    # load data cube and master calibration frames
    cube = CCDData.read(os.path.join(directory, cubename), unit=u.dimensionless_unscaled)
    master_bias = CCDData.read(os.path.join(cali_dir, biasname), unit=u.dimensionless_unscaled)
    master_dark = CCDData.read(os.path.join(cali_dir, darkname), unit=u.dimensionless_unscaled)
    if mode == 'monocolor':
        master_flat = CCDData.read(os.path.join(cali_dir, flatname), unit=u.dimensionless_unscaled)
    elif mode == 'multicolor':
        master_Cflat = CCDData.read(os.path.join(cali_dir, Cflat), unit=u.dimensionless_unscaled)
        master_gflat = CCDData.read(os.path.join(cali_dir, gflat), unit=u.dimensionless_unscaled)
        master_rflat = CCDData.read(os.path.join(cali_dir, rflat), unit=u.dimensionless_unscaled)
        master_iflat = CCDData.read(os.path.join(cali_dir, iflat), unit=u.dimensionless_unscaled)

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
        # flat is normalized within the ccdproc.flat_correct function
        if mode == 'monocolor':
            reduced_img = ccdproc.flat_correct(dark_subtracted, master_flat)
        elif mode == 'multicolor':
            filt = cube.header['FILTER']
            assert filt in ['CLEAR', 'g', 'r', 'i']
            if filt == 'CLEAR':
                reduced_img = ccdproc.flat_correct(dark_subtracted, master_Cflat)
            elif filt == 'g':
                reduced_img = ccdproc.flat_correct(dark_subtracted, master_gflat)
            elif filt == 'r':
                reduced_img = ccdproc.flat_correct(dark_subtracted, master_rflat)
            elif filt == 'i':
                reduced_img = ccdproc.flat_correct(dark_subtracted, master_iflat)
        reduced_cube.append(reduced_img)

    # create 'reduced_cubes' subdirectory
    if not os.path.exists(out_dir):
        print('Created', out_dir)
        os.mkdir(out_dir)
    
    # write to 3D cube in 'reduced_cubes', name=cubename_reduced, header=original cube header
    reduced_cube = np.asarray(reduced_cube)
    reduced_cube = reduced_cube.astype('float32')
    fits.writeto(os.path.join(out_dir, cubename[:-5]+'_reduced.fits'), 
                    reduced_cube, header=cube.header, overwrite=True)

    return


def solve_img(directory, imgname):
    """
    Submit to astrometry.net and solve image, return a WCS header
    """
    wcs_header = ast.solve_from_image(os.path.join(directory, imgname), force_image_upload=True)
    return wcs_header


def align_frames(cube, ref_img, nframes):
    """
    Align every frame in cube to ref_img.
    """
    # convolve ref_img before astroalign source detection
    # convolved data only used for alignment; original data is used for writing to new cube
    sigma = 3.0 * gaussian_fwhm_to_sigma  # FWHM = 3.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    convolved_ref_img = convolve(ref_img, kernel, normalize_kernel=True)

    # convolve and align every frame to ref_img and obtain transformation parameters
    # reproject original frames using transformation parameters
    aligned_list = [ref_img]
    for n in range(1, nframes):
        img_to_align = cube[n]
        convolved_img_to_align = convolve(img_to_align, kernel, normalize_kernel=True)
        # try finding alignment transformation p. If failed, no transformation will be applied
        try: 
            p, (pos_img, pos_img_rot) = aa.find_transform(convolved_img_to_align, 
                                                        convolved_ref_img, detection_sigma=2.0)
        except:
            p = None
            aligned_img = img_to_align
            print("WARNING: failed to align frame", n, ", no transformation is applied.")
        
        if p != None:
            aligned_img, footprint = aa.apply_transform(p, img_to_align, ref_img)
            print("\nRotation: {:.2f} degrees".format(p.rotation * 180.0 / np.pi))
            print("Scale factor: {:.4f}".format(p.scale))
            print("Translation: (x, y) = ({:.2f}, {:.2f})".format(*p.translation))

        aligned_list.append(aligned_img)
        print('current array len =', len(aligned_list))

    # correct data type for aligned frames
    aligned_array = np.array(aligned_list)
    aligned_array = aligned_array.astype('float32')
    print('final array shape =', aligned_array.shape)

    return aligned_array




def solve_and_align(directory, cubename, out_dir):
    """
    Solve the first frame of a 3D cube, align other frames to it.
    Output an aligned WCS cube to out_dir. 
    Note the wait time for astrometry.net varies significantly with time of day.
    """
    # read from 3D cube and obtain its 1st frame for solving and as reference frame
    # IMPORTANT: np.float32 is to reset endian-ness for astroalign, do not modify
    # FITS from fits.writeto has a different endian-ness that is incompatable with astroalign
    # otherwise will raises ValueError: Big-endian buffer not supported on little-endian compiler
    # This issue is introduced in new versions of scikit-image that astroalign relies on
    hdulist = fits.open(os.path.join(directory, cubename))
    hdu = hdulist[0]
    cube = np.float32(hdu.data)
    ref_img = cube[0]

    # obtain essential cards from cube header
    nframes = int(hdu.header['NAXIS3'])
    assert nframes == len(cube)
    exptime = float(hdu.header['EXPTIME'])
    dateobs = hdu.header['DATE-OBS']

    # write 1st frame of cube to a temporary 2D fits and solve it
    fits.writeto(os.path.join(directory, 'solve_temporary.fits'), ref_img, header=hdu.header, overwrite=True)
    wcs_header = solve_img(directory, 'solve_temporary.fits')
    os.remove(os.path.join(directory, 'solve_temporary.fits'))

    # align frames to ref_img
    aligned_array = align_frames(cube=cube, ref_img=ref_img, nframes=nframes)

    # add essential cards from original header to WCS header; write to WCS fits
    wcs_header.set('EXPTIME', exptime, before='CTYPE1')
    wcs_header.set('DATE-OBS', dateobs, before='CTYPE1')

    # create out_dir directory and write WCS cube to it
    if not os.path.exists(out_dir):
        print('Created', out_dir)
        os.mkdir(out_dir)
    fits.writeto(os.path.join(out_dir, cubename[:-5]+'_aligned.fits'), 
                    aligned_array, header=wcs_header, overwrite=True)

    return


# %%
