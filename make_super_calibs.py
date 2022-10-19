# %%
# This is the pipeline to create super calibration frames from raw TMO cubes.
# !!! Use astroconda environment for this pipeline !!!
from essentials import *

# Input calibration directory here
cali_dir='/Users/danny/Mirror/ASTRO/Calibration/CALIBS_20220921'


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