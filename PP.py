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

RA_Dec = '16h38m19.65s +03d48m52.0s'
target_sky_aperture, target_sky_annulus = sky_aperture_from_RADec(wcs=wcs, RA_Dec=RA_Dec, 
                                    radius=15.0, annulus=True, r_in=30.0, r_out=45.0)

ref_sky_aperture, ref_sky_annulus = detect_stars(img, wcs, detection_threshold=5, 
                                    radius=15.0, annulus=True, r_in=30.0, r_out=45.0, 
                                    xbounds=[100,1250], ybounds=[150,1400], 
                                    mask_coord=RA_Dec, 
                                    visu=True, visu_dir='/Users/danny/Desktop', high_res=True)


# %%
target_tablename = os.path.join(out_dir, 'target_flux.ecsv')
ref_tablename = os.path.join(out_dir, 'ref_flux.ecsv')

assert os.path.exists(target_tablename) == False
assert os.path.exists(ref_tablename) == False

for cubename in cubelist:
    print('\n---------------------------------------------')
    print('Performing photometry on', cubelist.index(cubename)+1, 'th cube')
    print('---------------------------------------------')

    target_ts, ref_ts = dual_thread_photometry(directory=directory, cubename=cubename,
        target_sky_aperture=target_sky_aperture, ref_sky_aperture=ref_sky_aperture, 
        annulus=True, target_sky_annulus=target_sky_annulus, ref_sky_annulus=ref_sky_annulus)

    if not os.path.exists(target_tablename):
        # write target ts table to ecsv
        target_ts.write(target_tablename, overwrite=False)
    else:
        existing_target_ts = TimeSeries.read(target_tablename, time_column='time')
        new_target_ts = vstack([existing_target_ts, target_ts])
        new_target_ts.write(target_tablename, overwrite=True)
    
    if not os.path.exists(ref_tablename):
        # write ref ts table to ecsv
        ref_ts.write(ref_tablename, overwrite=False)
    else:
        existing_ref_ts = TimeSeries.read(ref_tablename, time_column='time')
        new_ref_ts = vstack([existing_ref_ts, ref_ts])
        new_ref_ts.write(ref_tablename, overwrite=True)



# %%
table, RFM = LCs_visualizer(directory=out_dir, visu_dir='/Users/danny/Desktop', 
                                    mode='full', layout=[4,4])


# %%
reflist = [12, 10, 8, 11, 7]

diff_lc_name = os.path.join(out_dir, 'diff_lc.ecsv')

diff_lc = differential_photometry(directory=out_dir, reflist=reflist)
diff_lc.write(diff_lc_name, overwrite=True)


# %%