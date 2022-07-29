# %%
import os
from glob import glob1
import numpy as np
from astropy import units as u
from astropy.nddata import CCDData
import ccdproc
from ccdproc import Combiner
from astropy.io import fits


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
    fits.writeto(os.path.join(directory, outname), \
                    combined_average, header=cube.header, overwrite=True)
    return


def calibrate_cube(directory, cubename, cali_dir, darkname, biasname, flatname):
    """
    Reduce 3D cube by master dark, bias, flat; output to subdir 'reduced_cubes'.
    """
    # load data cube and master calibration frames
    cube = CCDData.read(os.path.join(directory, cubename), unit=u.dimensionless_unscaled)
    master_dark = CCDData.read(os.path.join(cali_dir, darkname), unit=u.dimensionless_unscaled)
    master_bias = CCDData.read(os.path.join(cali_dir, biasname), unit=u.dimensionless_unscaled)
    master_flat = CCDData.read(os.path.join(cali_dir, flatname), unit=u.dimensionless_unscaled)

    # subtract bias, dark, divide by flat
    bias_subtracted = ccdproc.subtract_bias(cube, master_bias)
    dark_subtracted = ccdproc.subtract_dark(bias_subtracted, master_dark)
    reduced_cube = ccdproc.flat_correct(dark_subtracted, master_flat)

    # create 'reduced_cubes' subdirectory
    if not os.path.exists(os.path.join(directory, 'reduced_cubes')):
        print('Created /reduced_cubes subfolder')
        os.mkdir(os.path.join(directory, 'reduced_cubes'))
    
    # write to 3D cube in 'reduced_cubes', name=cubename_reduced, header=original cube header
    reduced_cube = np.asarray(reduced_cube)
    reduced_cube = reduced_cube.astype('float32')
    fits.writeto(os.path.join(directory,'reduced_cubes',cubename[:-5]+'_reduced.fits'), \
				reduced_cube, header=cube.header, overwrite=True)
    return


    




            
# %%
