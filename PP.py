# %%
import os
from glob import glob1
from astropy.io import fits
from astropy.wcs import WCS
from phot_functions import estimate_background, detect_stars


directory = '/Users/danny/Mirror/ASTRO/JPL_NEO/Contact_Binary/data/CSS_034852/20220719/aligned_cubes'
cubelist = sorted(glob1(directory, '*_aligned.fits'))

cubename = cubelist[0]


# %%

hdulist = fits.open(os.path.join(directory, cubename))
hdu = hdulist[0]

img = hdu.data[0]

# skymap = estimate_background(img, mode='2D', visu=True, visu_dir='/Users/danny/Desktop', high_res=False)
wcs = WCS(hdu.header, naxis=2)
sky_apertures, sky_annulus = detect_stars(img, wcs, radius=15, r_in=25, r_out=35, 
                                        xbounds=[100,1400], ybounds=[100,1400], 
                                        target_coord='16h38m19.65s +03d48m52.0s', 
                                        visu=True, visu_dir='/Users/danny/Desktop', high_res=True)

# %%
