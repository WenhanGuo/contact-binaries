# %%
import os
from glob import glob1
import numpy as np
from astropy import units as u
from astropy.nddata import CCDData
from ccdproc import Combiner
from astropy.io import fits


# %%
Directory = '/Users/danny/Mirror/ASTRO/JPL_NEO/Calibration'
CubeName = 'bias_20220726_205921.fits'

cube = CCDData.read(os.path.join(Directory, CubeName), unit=u.dimensionless_unscaled)
# nframes = cube.header['NAXIS3']


# %%
# Create ccdprocs combiner object
combiner = Combiner(cube)

# Average combine with sigma clipping
combiner.sigma_clipping(low_thresh=3, high_thresh=3)
combined_average = combiner.average_combine()

fits.writeto(os.path.join(Directory, 'master_bias_20220726.fits'), \
            np.asarray(combined_average), overwrite=True)

            
# %%
