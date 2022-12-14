# %%
import os
import numpy as np
import glob
from glob import glob1
import shutil

from astropy import units as u
from astropy.nddata import CCDData
from astropy.io import fits
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.convolution import Gaussian2DKernel, convolve

from astroquery.astrometry_net import AstrometryNet
# import astroalign as aa
import alipy

import sys
import six
sys.modules['astropy.extern.six'] = six
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


def calibrate_flat(cali_dir, flatname, darkbiasname):
    """
    Subtract dark and bias from flat, normalize.
    """
    # load master calibration frames
    darkbias = CCDData.read(os.path.join(cali_dir, darkbiasname), unit=u.dimensionless_unscaled)
    flat = CCDData.read(os.path.join(cali_dir, flatname), unit=u.dimensionless_unscaled)

    # subtract master_dark and master_bias from master_flat
    # 10s are dummy variables
    corrected = ccdproc.subtract_dark(flat, darkbias, 
                                                data_exposure=10*u.s, dark_exposure=10*u.s)
    # normalize flat
    normalized = corrected / np.median(corrected)

    fits.writeto(os.path.join(cali_dir, flatname[:-5]+'_norm.fits'), 
                    normalized, header=flat.header, overwrite=True)
    


def calibrate_cube(directory, cubename, cali_dir, darkbiasname, darkexptime, flatname=None, 
            mode='multicolor', Cflat=None, gflat=None, rflat=None, iflat=None, out_dir=None):
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
    darkbias = CCDData.read(os.path.join(cali_dir, darkbiasname), unit=u.dimensionless_unscaled)
    if mode == 'monocolor':
        master_flat = CCDData.read(os.path.join(cali_dir, flatname), unit=u.dimensionless_unscaled)
    elif mode == 'multicolor':
        master_Cflat = CCDData.read(os.path.join(cali_dir, Cflat), unit=u.dimensionless_unscaled)
        master_gflat = CCDData.read(os.path.join(cali_dir, gflat), unit=u.dimensionless_unscaled)
        master_rflat = CCDData.read(os.path.join(cali_dir, rflat), unit=u.dimensionless_unscaled)
        master_iflat = CCDData.read(os.path.join(cali_dir, iflat), unit=u.dimensionless_unscaled)

    # subtract bias, dark, divide by flat
    reduced_cube = []
    nframes = int(cube.header['NAXIS3'])
    exptime = float(cube.header['EXPTIME'])
    for n in range(nframes):
        img = cube[n]
        darkbias_subtracted = ccdproc.subtract_dark(img, darkbias, 
                                        data_exposure=exptime*u.s, dark_exposure=darkexptime*u.s)
        # flat is normalized within the ccdproc.flat_correct function
        if mode == 'monocolor':
            reduced_img = ccdproc.flat_correct(darkbias_subtracted, master_flat)
        elif mode == 'multicolor':
            filt = cube.header['FILTER']
            assert filt in ['CLEAR', 'g', 'r', 'i']
            if filt == 'CLEAR':
                reduced_img = ccdproc.flat_correct(darkbias_subtracted, master_Cflat)
            elif filt == 'g':
                reduced_img = ccdproc.flat_correct(darkbias_subtracted, master_gflat)
            elif filt == 'r':
                reduced_img = ccdproc.flat_correct(darkbias_subtracted, master_rflat)
            elif filt == 'i':
                reduced_img = ccdproc.flat_correct(darkbias_subtracted, master_iflat)
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


def slice_cube(directory, cubename, out_dir):
    """
    Time-increment slicing of FITS cube.
    If the cube is a Time Series cube, (has date-obs and exptime in header), 
    increment slices' headers so that their date-obs reflect actual frame exposure time.
    """
    # IMPORTANT: np.float32 is to reset endian-ness for astroalign, do not modify
    # FITS from fits.writeto has a different endian-ness that is incompatable with astroalign
    # otherwise will raises ValueError: Big-endian buffer not supported on little-endian compiler
    # This issue is introduced in new versions of scikit-image that astroalign relies on
    hdulist = fits.open(os.path.join(directory, cubename))
    hdu = hdulist[0]
    cube = np.float32(hdu.data)

    # obtain essential cards from cube header
    nframes = int(hdu.header['NAXIS3'])
    assert nframes == len(cube)
    exptime = float(hdu.header['EXPTIME'])
    dateobs = hdu.header['DATE-OBS']

    for n in range(nframes):
        img = cube[n]
        header = hdu.header
        actual_dateobs = np.datetime64(dateobs) + n * int(exptime)
        actual_dateobs = str(actual_dateobs)+'+00:00'
        header.set('DATE-OBS', actual_dateobs)
        
        # create out_dir directory and write WCS cube to it
        if not os.path.exists(out_dir):
            print('Created', out_dir)
            os.mkdir(out_dir)
        fits.writeto(os.path.join(out_dir, cubename[:-5]+'_sliced'+str(n+1)+'.fits'), 
                    img, header=hdu.header, overwrite=True)
    
    return


def solve_img(directory, imgname):
    """
    Submit to astrometry.net and solve image, return a WCS header
    """
    wcs_header = ast.solve_from_image(os.path.join(directory, imgname), force_image_upload=True)
    return wcs_header


def alipy_ident(directory, refname=None):
    """
    Adapted from do_alipy3. 
    Align 2D images in directory to ref.fits
    """
    images_to_align = sorted(glob.glob(os.path.join(directory, '*.fits')))
    if refname == None:
        ref_image = os.path.join(directory, images_to_align[0])
    else:
        ref_image = os.path.join(directory, refname)
    print('ref_image path =', ref_image)

    identifications = alipy.ident.run(ref_image, images_to_align, visu=False)
    # Put visu=True to get visualizations in form of png files (nice but much slower)
    # On multi-extension data, you will want to specify the hdu (see API doc).

    return ref_image, identifications



def alipy_align(ref_image, identifications, out_dir):

    success = 0
    # The output is a list of Identification objects, which contain the transforms :
    for id in identifications: # list of the same length as images_to_align.
        if id.ok == True: # i.e., if it worked
            success += 1
    print('Identification successful', success, '/', len(identifications), 'frames')
    print('--------------------------------------------------')

    outputshape = alipy.align.shape(ref_image)
    # This is simply a tuple (width, height)... you could specify any other shape.

    for id in identifications:
        if id.ok == True:
            # Align using scipy and the simple affine transorm :
            alipy.align.affineremap(id.ukn.filepath, id.trans, shape=outputshape, makepng=False)

    # By default, the aligned images are written into a directory "alipy_out".
    # Move alipy_out to desired output dir.
    if os.path.exists(out_dir):
        os.rmdir(out_dir)
    shutil.copytree('./alipy_out', out_dir)
    shutil.rmtree('./alipy_out')

    return



def solve_and_migrate_header(directory, refname=None):
    """
    Solve an image.
    Apply the WCS transformation cards to all images in folder.
    Preserve each slice's other header info.
    Needs all slices already aligned by alipy_align()
    """
    imgnames = sorted(glob.glob(os.path.join(directory, '*.fits')))

    if refname == None:
        refname = sorted(glob1(directory, '*fits'))[0]
    print('ref image name =', refname)

    wcs_header = solve_img(directory=directory, imgname=refname)
    del wcs_header['COMMENT']

    for n in range(len(imgnames)):
        with fits.open(imgnames[n]) as hdulist:
            hdu = hdulist[0]
            hdr = hdu.header
            # obtain essential headers from each slice
            exptime = float(hdu.header['EXPTIME'])
            dateobs = hdu.header['DATE-OBS']
            filt = hdu.header['FILTER']
            sky = hdu.header['SKY']
            # append essential headers to WCS header, overwrite file
            new_hdr = wcs_header
            new_hdr.set('EXPTIME', exptime, before='CTYPE1')
            new_hdr.set('DATE-OBS', dateobs, before='CTYPE1')
            new_hdr.set('FILTER', filt, before='CTYPE1')
            new_hdr.set('SKY', sky, before='CTYPE1')
            fits.writeto(imgnames[n], hdu.data, header=new_hdr, overwrite=True)
        print('Header rewrite complete', n+1, '/', len(imgnames), 'frames')

    return


# ------------------------------------ LEGACY FUNCTIONS ------------------------------------
'''
def legacy_align_frames(cube, ref_img, nframes, do_convolve=True):
    """
    This is legacy align_frames function using astroalign. 
    Align every frame in cube to ref_img.
    """
    if do_convolve == True:
        # convolve ref_img before astroalign source detection
        # convolved data only used for alignment; original data is used for writing to new cube
        sigma = 3.0 * gaussian_fwhm_to_sigma  # FWHM = 3.
        kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
        convolved_ref_img = convolve(ref_img, kernel, normalize_kernel=True)

        # convolve and align every frame to ref_img and obtain transformation parameters
        # reproject original frames using transformation parameters
        aligned_list = []
        for n in range(nframes):
            img_to_align = cube[n]
            convolved_img_to_align = convolve(img_to_align, kernel, normalize_kernel=True)
            # try finding alignment transformation p. If failed, no transformation will be applied
            try: 
                p, (pos_img, pos_img_rot) = aa.find_transform(convolved_img_to_align, 
                                                            convolved_ref_img, detection_sigma=1.0)
            except:
                p = None
                aligned_img = img_to_align
                print("WARNING: failed to align frame", n, ", no transformation is applied.")
    elif do_convolve == False:
        aligned_list = []
        for n in range(nframes):
            img_to_align = cube[n]
            try: 
                p, (pos_img, pos_img_rot) = aa.find_transform(img_to_align, ref_img, detection_sigma=2.0)
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




def legacy_solve_and_align(directory, cubename, out_dir, do_convolve=True):
    """
    This is legacy solve_and_align function using astroalign.
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
    filt = hdu.header['FILTER']

    # write 1st frame of cube to a temporary 2D fits and solve it
    fits.writeto(os.path.join(directory, 'solve_temporary.fits'), ref_img, header=hdu.header, overwrite=True)
    wcs_header = solve_img(directory, 'solve_temporary.fits')
    os.remove(os.path.join(directory, 'solve_temporary.fits'))

    # align frames to ref_img
    aligned_array = legacy_align_frames(cube=cube, ref_img=ref_img, nframes=nframes, do_convolve=do_convolve)

    # add essential cards from original header to WCS header; write to WCS fits
    wcs_header.set('EXPTIME', exptime, before='CTYPE1')
    wcs_header.set('DATE-OBS', dateobs, before='CTYPE1')
    wcs_header.set('FILTER', filt, before='CTYPE1')

    # create out_dir directory and write WCS cube to it
    if not os.path.exists(out_dir):
        print('Created', out_dir)
        os.mkdir(out_dir)
    fits.writeto(os.path.join(out_dir, cubename[:-5]+'_aligned.fits'), 
                    aligned_array, header=wcs_header, overwrite=True)

    return


# %%
def legacy_batch_solve_and_align(directory, out_dir, do_convolve=True):
    """
    This is legacy batch_solve_and_align using astroalign.
    Solve and align every cube in folder to ref_img.
    """
    cubelist = sorted(glob1(directory, '*.fits'))

    # compact code to set ref_img = 1st frame of 1st cube; equivalent code broken down below
    ref_img = np.float32(fits.open(os.path.join(directory, cubelist[0]))[0].data)[0]
    
    # write ref_img to a temporary 2D fits and solve it
    fits.writeto(os.path.join(directory, 'solve_temporary.fits'), ref_img, overwrite=True)
    wcs_header = solve_img(directory, 'solve_temporary.fits')
    os.remove(os.path.join(directory, 'solve_temporary.fits'))


    for cubename in cubelist:
        # IMPORTANT: np.float32 is to reset endian-ness for astroalign, do not modify
        # FITS from fits.writeto has a different endian-ness that is incompatable with astroalign
        # otherwise will raises ValueError: Big-endian buffer not supported on little-endian compiler
        # This issue is introduced in new versions of scikit-image that astroalign relies on
        hdulist = fits.open(os.path.join(directory, cubename))
        hdu = hdulist[0]
        cube = np.float32(hdu.data)

        # obtain essential cards from cube header
        nframes = int(hdu.header['NAXIS3'])
        assert nframes == len(cube)
        exptime = float(hdu.header['EXPTIME'])
        dateobs = hdu.header['DATE-OBS']
        filt = hdu.header['FILTER']

        # align frames to ref_img
        aligned_array = legacy_align_frames(cube=cube, ref_img=ref_img, nframes=nframes, do_convolve=do_convolve)

        # add essential cards from original header to WCS header; write to WCS fits
        wcs_header.set('EXPTIME', exptime, before='CTYPE1')
        wcs_header.set('DATE-OBS', dateobs, before='CTYPE1')
        wcs_header.set('FILTER', filt, before='CTYPE1')

        # create out_dir directory and write WCS cube to it
        if not os.path.exists(out_dir):
            print('Created', out_dir)
            os.mkdir(out_dir)
        fits.writeto(os.path.join(out_dir, cubename[:-5]+'_aligned.fits'), 
                        aligned_array, header=wcs_header, overwrite=True)

    return
'''
