# %%
# !!! Use astroconda environment for this pipeline !!!
import os
from essentials import *
from glob import glob1

# Input data directory here
directory = '/Users/danny/Mirror/ASTRO/ASTR101/lab2/n3s1data'
# directory = '/Users/danny/Mirror/ASTRO/Contact_Binary/data/CSS_cgri_test'
cubelist = sorted(glob1(directory, '*.fits'))

# Input super calibs directory here
cali_dir='/Users/danny/Mirror/ASTRO/Calibration/CALIBS_20220921'


# %%
# Calibrate cubes with dark, bias, and colored flats
for cubename in cubelist:
    calibrate_cube(directory=directory, cubename=cubename, 
                    cali_dir=cali_dir, 
                    darkbiasname='Master_DarkBias.fits', 
                    darkexptime=10.0, 
                    mode='multicolor',
                    Cflat='Master_Flat_C_norm.fits',
                    gflat='Master_Flat_g_norm.fits',
                    rflat='Master_Flat_r_norm.fits',
                    iflat='Master_Flat_i_norm.fits',
                    out_dir=os.path.join(directory, 'reduced_cubes'))


# %%
# Slice cubes and increment timestamps
redu_cubelist = sorted(glob1(directory+'/reduced_cubes', '*.fits'))
for cubename in redu_cubelist:
    slice_cube(directory=os.path.join(directory, 'reduced_cubes'), 
                cubename=cubename,
                out_dir=os.path.join(directory, 'sliced'))


# %%
# Solve frame and align all images; migrate WCS headers
alipy_align(directory=os.path.join(directory, 'sliced'),
            out_dir=os.path.join(directory, 'aligned'))

solve_and_migrate_header(directory=os.path.join(directory, 'aligned'))


# %%
# ------------------------------------ LEGACY FUNCTIONS ------------------------------------

"""
reduced_cubelist = sorted(glob1(os.path.join(directory, 'reduced_cubes'), '*_reduced.fits'))

for cubename in reduced_cubelist: 
    solve_and_align(directory=os.path.join(directory, 'reduced_cubes'), 
                    cubename=cubename, 
                    out_dir=os.path.join(directory, 'aligned_cubes'))

batch_solve_and_align(directory=os.path.join(directory, 'reduced_cubes'), 
                        out_dir=os.path.join(directory, 'aligned_cubes'),
                        do_convolve=False)
"""
