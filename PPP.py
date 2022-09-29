# %%
import os
from essentials import avg_cube, calibrate_cube, solve_and_align
from glob import glob1


directory = '/Users/danny/Mirror/ASTRO/Contact_Binary/data/CSS_034852/cgri_test'
cali_dir='/Users/danny/Mirror/ASTRO/Calibration/CALIBS_20220921'
cubelist = sorted(glob1(directory, '*.fits'))


# %%
# For creating sigma-clipped master frames
'''
avg_cube(directory=cali_dir, 
        cubename='Darks_20220921_001521.fits',
        outname='Master_Dark_20220921.fits')
avg_cube(directory=cali_dir, 
        cubename='Biases_20220921_004356.fits',
        outname='Master_Bias_20220921.fits')
avg_cube(directory=cali_dir, 
        cubename='twilights_7C_20220920_190121.fits',
        outname='Master_Flat_C_20220921.fits')
avg_cube(directory=cali_dir, 
        cubename='twilights_4g_20220920_185345.fits',
        outname='Master_Flat_g_20220921.fits')
avg_cube(directory=cali_dir, 
        cubename='twilights_2r_20220920_184356.fits',
        outname='Master_Flat_r_20220921.fits')
avg_cube(directory=cali_dir, 
        cubename='twilights_1i_20220920_183559.fits',
        outname='Master_Flat_i_20220921.fits')
'''

# %%
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
