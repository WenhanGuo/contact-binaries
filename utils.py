# %%
import os
from glob import glob1
import numpy as np
from astropy import units as u
from astropy.nddata import CCDData
from ccdproc import Combiner
from astropy.io import fits


def avg_cube(directory, cubename, outname):
    """
    Create 2D fits img from sigma-clipped mean of a 3D cube.
    """
    cube = CCDData.read(os.path.join(directory, cubename), unit=u.dimensionless_unscaled)

    # create ccdprocs combiner object
    combiner = Combiner(cube)

    # average combine with sigma clipping
    combiner.sigma_clipping(low_thresh=3, high_thresh=3)
    combined_average = combiner.average_combine()
    combined_average = np.asarray(combined_average)
    combined_average = combined_average.astype('float32')

    # write to fits with name=outname, header=original cube header
    fits.writeto(os.path.join(directory, outname), \
                    combined_average, header=cube.header, overwrite=True)
    return







            
# %%
