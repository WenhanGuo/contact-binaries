# %%
import os
from essentials import *
from glob import glob1


directory = '/Users/danny/Mirror/ASTRO/ASTR101/lab2/data'
# directory = '/Users/danny/Mirror/ASTRO/Contact_Binary/data/CSS_034852/20220719'
cali_dir='/Users/danny/Mirror/ASTRO/Calibration/CALIBS_20220921'
cubelist = sorted(glob1(directory, '*.fits'))


# %%
# For creating sigma-clipped master frames from calibration cubes

avg_cube(directory=cali_dir, 
        cubename='Darks_20220921_001521.fits',
        outname='Master_DarkBias.fits')
avg_cube(directory=cali_dir, 
        cubename='Biases_20220921_004356.fits',
        outname='Master_Bias.fits')
avg_cube(directory=cali_dir, 
        cubename='twilights_7C_20220920_190121.fits',
        outname='Master_Flat_C.fits')
avg_cube(directory=cali_dir, 
        cubename='twilights_4g_20220920_185345.fits',
        outname='Master_Flat_g.fits')
avg_cube(directory=cali_dir, 
        cubename='twilights_2r_20220920_184356.fits',
        outname='Master_Flat_r.fits')
avg_cube(directory=cali_dir, 
        cubename='twilights_1i_20220920_183559.fits',
        outname='Master_Flat_i.fits')


# %%
# calibrate master dark and master flats

calibrate_flat(cali_dir=cali_dir, 
                flatname='Master_Flat_C.fits',
                darkbiasname='Master_DarkBias.fits')
calibrate_flat(cali_dir=cali_dir, 
                flatname='Master_Flat_g.fits',
                darkbiasname='Master_DarkBias.fits')
calibrate_flat(cali_dir=cali_dir, 
                flatname='Master_Flat_r.fits',
                darkbiasname='Master_DarkBias.fits')
calibrate_flat(cali_dir=cali_dir, 
                flatname='Master_Flat_i.fits',
                darkbiasname='Master_DarkBias.fits')


# %%
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
for cubename in cubelist:
    slice_cube(directory=directory, 
                cubename=cubename,
                out_dir=os.path.join(directory, 'sliced'))

# %%

alipy_align(directory=os.path.join(directory, 'sliced'),
            out_dir=os.path.join(directory, 'aligned'))


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
                        do_convolve=False)"""

