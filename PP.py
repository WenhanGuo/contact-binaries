# %%
# This is the Photometry Pipeline for obtaining light curves from calibrated data
# Input: WCS-aligned frames from PPP
# Output: light curve data and visualizations
# !!! Use photometry environment (python 3.8) for this pipeline !!!
import os
from glob import glob1
from astropy.io import fits
from astropy.wcs import WCS
from phot_functions import *
from astropy.timeseries import TimeSeries
from astropy.table import vstack
import matplotlib.pyplot as plt
# %matplotlib widget

# Input alignd data directory here
# directory = '/Users/danny/My_Root/ASTRO/Contact_Binary/data/CSS_cgri_test/aligned'
directory = '/Users/danny/My_Root/ASTRO/ASTR101/final_project/data/20221030/aligned'
imglist = sorted(glob1(directory, '*.fits'))

# Input desired visualization path (typically just data dir)
# out_dir = '/Users/danny/My_Root/ASTRO/Contact_Binary/data/CSS_cgri_test'
out_dir = '/Users/danny/My_Root/ASTRO/ASTR101/final_project/data/20221030'



# %%
# Source Detection; define apertures and annuli
hdulist = fits.open(os.path.join(directory, imglist[0]))
hdu = hdulist[0]
img = hdu.data
wcs = WCS(hdu.header, naxis=2)

# Input target object's RA and Dec
# RA_Dec = '16h38m19.65s +03d48m52.0s'   # for CSS_034852
# RA_Dec = '22h37m47.92s +01d32m02.25s'   # for astro lab2 n3s2
# RA_Dec = '21h36m54.12s +02d28m53.5s'   # for CSS cgri test
# RA_Dec = '19h16m03.66s +30d15m40.6s'   # for astro lab2 n3s1
RA_Dec = '18h13m11.13s +42d51m50.4s'   # for astro final project

# Define target object's aperture and annuli
target_sky_aperture, target_sky_annulus = sky_aperture_from_RADec(wcs=wcs, RA_Dec=RA_Dec, 
                                    radius=15.0, annulus=True, r_in=30.0, r_out=45.0)

# Source detection; define reference stars' apertures and annuli
# Run this cell multiple times to try out different detection sigma, xybounds, aperture & annulus sizes
ref_sky_aperture, ref_sky_annulus = detect_stars(img, wcs, skymap=True, detection_threshold=8, 
                                    radius=15.0, annulus=True, r_in=30.0, r_out=45.0, 
                                    xbounds=[100,1200], ybounds=[150,1450], 
                                    mask_coord=RA_Dec, 
                                    visu=True, visu_dir=out_dir, high_res=True)


# %%
# Perform annulus photometry on all aligned images in directory
# Make sure photometry tables do NOT exist; we will create these with this cell
target_tablename = os.path.join(out_dir, 'target_flux.ecsv')
ref_tablename = os.path.join(out_dir, 'ref_flux.ecsv')
assert os.path.exists(target_tablename) == False
assert os.path.exists(ref_tablename) == False

# Photometry for every image in imglist; concatenate results to target_flux.ecsv and ref_flux.ecsv
for imgname in imglist:
    print('Performing photometry on', imglist.index(imgname)+1, '/', len(imglist), 'images')
    target_ts, ref_ts = dual_thread_photometry(directory=directory, cubename=imgname,
                                                target_sky_aperture=target_sky_aperture, 
                                                ref_sky_aperture=ref_sky_aperture, 
                                                annulus=True, 
                                                target_sky_annulus=target_sky_annulus, 
                                                ref_sky_annulus=ref_sky_annulus)
    update_table(ts=target_ts, mother_tablename=target_tablename)
    update_table(ts=ref_ts, mother_tablename=ref_tablename)

# Normalize both flux tables by star
normalize_target_table(directory=out_dir, target_table='target_flux.ecsv')
normalize_ref_table(directory=out_dir, ref_table='ref_flux.ecsv')


# %%
# everything below is ongoing -- I have to switch from mono-filter to cgri
"""
table, RFM = LCs_visualizer(directory=out_dir, visu_dir=out_dir, 
                                    mode='full', layout=[3,5])


# %%
reflist = [1, 2]

diff_lc_name = os.path.join(out_dir, 'diff_lc.ecsv')

diff_lc = differential_photometry(directory=out_dir, reflist=reflist)
diff_lc.write(diff_lc_name, overwrite=True)


# %%
# obj_dir = '/Users/danny/Mirror/ASTRO/JPL_NEO/Contact_Binary/data/CSS_034852'
obj_dir = '/Users/danny/Mirror/ASTRO/ASTR101/lab2/data'

fold_lc(obj_dir=obj_dir, table='diff_lc.ecsv')
"""
