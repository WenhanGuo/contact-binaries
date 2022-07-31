# %%
import os
from glob import glob1
from astropy.io import fits
from astropy.wcs import WCS
from phot_functions import *
from astropy.timeseries import TimeSeries
from astropy.table import vstack
import matplotlib.pyplot as plt


directory = '/Users/danny/Mirror/ASTRO/JPL_NEO/Contact_Binary/data/CSS_034852/20220719/aligned_cubes'
cubelist = sorted(glob1(directory, '*_aligned.fits'))

out_dir = '/Users/danny/Mirror/ASTRO/JPL_NEO/Contact_Binary/data/CSS_034852/20220719'


# %%
hdulist = fits.open(os.path.join(directory, cubelist[0]))
hdu = hdulist[0]
img = hdu.data[0]
wcs = WCS(hdu.header, naxis=2)

sky_aperture, sky_annulus = detect_stars(img, wcs, radius=15.0, r_in=30.0, r_out=45.0, 
                                        xbounds=[50,1450], ybounds=[50,1450], 
                                        target_coord='16h38m19.65s +03d48m52.0s', 
                                        visu=False, visu_dir='/Users/danny/Desktop', high_res=True)


# %%
for cubename in cubelist:
    target_ts = simple_photometry(directory, cubename, 
                target_coord='16h38m19.65s +03d48m52.0s', 
                radius=15.0, annulus=True, r_in=30.0, r_out=45.0)

    tablename = os.path.join(out_dir, 'target_flux.ecsv')
    if not os.path.exists(tablename):
        # write ts table to ecsv
        target_ts.write(tablename, overwrite=False)
    else:
        existing_target_ts = TimeSeries.read(tablename, time_column='time')
        new_target_ts = vstack([existing_target_ts, target_ts])
        new_target_ts.write(tablename, overwrite=True)


# %%
for cubename in cubelist:
    ref_ts = multithread_photometry(directory, cubename, 
                sky_aperture=sky_aperture, 
                annulus=True, sky_annulus=sky_annulus)

    tablename = os.path.join(out_dir, 'ref_flux.ecsv')
    if not os.path.exists(tablename):
        # write ts table to ecsv
        ref_ts.write(tablename, overwrite=False)
    else:
        existing_ref_ts = TimeSeries.read(tablename, time_column='time')
        new_ref_ts = vstack([existing_ref_ts, ref_ts])
        new_ref_ts.write(tablename, overwrite=True)


# %%
table, RFM = quicklook_lightcurve(directory=out_dir, visu_dir='/Users/danny/Desktop', 
                                    layout=[4,6], normalize=True)


# %%
