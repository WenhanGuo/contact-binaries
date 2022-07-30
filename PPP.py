# %%
import os
from essentials import calibrate_cube, solve_and_align
from glob import glob1
import astropy.units as u


directory = '/Users/danny/Mirror/ASTRO/JPL_NEO/Contact_Binary/data/CSS_034852/20220719'
cubelist = sorted(glob1(directory, '*.fits'))


for cubename in cubelist:
    calibrate_cube(directory=directory, cubename=cubename, 
                    cali_dir='/Users/danny/Mirror/ASTRO/JPL_NEO/Calibration', 
                    biasname='master_bias_20220726.fits', 
                    darkname='master_dark_20220726.fits', 
                    darkexptime=10.0, 
                    flatname='master_flat11_20220713.fits', 
                    out_dir=os.path.join(directory, 'reduced_cubes'))


# %%

reduced_cubelist = sorted(glob1(os.path.join(directory, 'reduced_cubes'), '*_reduced.fits'))

for cubename in reduced_cubelist: 
    solve_and_align(directory=os.path.join(directory, 'reduced_cubes'), 
                    cubename=cubename, 
                    out_dir=os.path.join(directory, 'aligned_cubes'))


# %%
