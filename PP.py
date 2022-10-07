# %%
import os
from glob import glob1
from astropy.io import fits
from astropy.wcs import WCS
from phot_functions import *
from astropy.timeseries import TimeSeries
from astropy.table import vstack
import matplotlib.pyplot as plt
 I want to update

directory = '/Users/danny/Mirror/ASTRO/Contact_Binary/data/CSS_cgri_test/aligned'
# directory = '/Users/danny/Mirror/ASTRO/ASTR101/lab2/data/aligned'
cubelist = sorted(glob1(directory, '*.fits'))

out_dir = '/Users/danny/Mirror/ASTRO/Contact_Binary/data/CSS_cgri_test'
# out_dir = '/Users/danny/Mirror/ASTRO/ASTR101/lab2/data'


# %%
hdulist = fits.open(os.path.join(directory, cubelist[0]))
hdu = hdulist[0]
# img = hdu.data[0]
img = hdu.data
wcs = WCS(hdu.header, naxis=2)

# RA_Dec = '16h38m19.65s +03d48m52.0s'   # for CSS_034852
# RA_Dec = '22h37m47.92s +01d32m02.25s'   # for astro lab2
RA_Dec = '21h36m54.12s +02d28m53.8s'
target_sky_aperture, target_sky_annulus = sky_aperture_from_RADec(wcs=wcs, RA_Dec=RA_Dec, 
                                    radius=15.0, annulus=True, r_in=30.0, r_out=45.0)

ref_sky_aperture, ref_sky_annulus = detect_stars(img, wcs, skymap=True, detection_threshold=8, 
                                    radius=15.0, annulus=True, r_in=30.0, r_out=45.0, 
                                    xbounds=[100,1000], ybounds=[150,1200], 
                                    mask_coord=RA_Dec, 
                                    visu=True, visu_dir=out_dir, high_res=True)


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
                                                target_sky_aperture=target_sky_aperture, 
                                                ref_sky_aperture=ref_sky_aperture, 
                                                annulus=True, 
                                                target_sky_annulus=target_sky_annulus, 
                                                ref_sky_annulus=ref_sky_annulus)
    update_table(ts=target_ts, mother_tablename=target_tablename)
    update_table(ts=ref_ts, mother_tablename=ref_tablename)

normalize_target_table(directory=out_dir, target_table='target_flux.ecsv')
normalize_ref_table(directory=out_dir, ref_table='ref_flux.ecsv')


# %%
table, RFM = LCs_visualizer(directory=out_dir, visu_dir='/Users/danny/Desktop', 
                                    mode='full', layout=[2,2])


# %%
reflist = [1, 2]

diff_lc_name = os.path.join(out_dir, 'diff_lc.ecsv')

diff_lc = differential_photometry(directory=out_dir, reflist=reflist)
diff_lc.write(diff_lc_name, overwrite=True)


# %%
# obj_dir = '/Users/danny/Mirror/ASTRO/JPL_NEO/Contact_Binary/data/CSS_034852'
obj_dir = '/Users/danny/Mirror/ASTRO/ASTR101/lab2/data'

fold_lc(obj_dir=obj_dir, table='diff_lc.ecsv')





# %%