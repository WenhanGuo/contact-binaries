# %%
# This is the Pre-Process Pipeline for calibrating, solving, and aligning cubes
# Input: raw data cubes in folder; super dark, bias, flats (multiple if filtered)
# Output: aligned frames with WCS headers

# !!! Use astroconda environment (python 3.7) for this pipeline !!!
# alipy3-header version needs to be installed in the environment
import os
from essentials import *
from glob import glob1

# Input data directory here
directory = '/Volumes/TMO_Data_4TB/cb_data/C01+13/20221106_newflat'
cubelist = sorted(glob1(directory, '*.fits'))

# Input super calibs directory here
cali_dir='/Users/TMObserver/Documents/CALIBS_20221106'


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

# for cubename in cubelist:
#     calibrate_cube(directory=directory, cubename=cubename, 
#                     cali_dir=cali_dir, 
#                     darkbiasname='Master_DarkBias.fits', 
#                     darkexptime=10.0, 
#                     mode='monocolor', 
#                     flatname='Master_Flat_C_norm.fits',
#                     out_dir=os.path.join(directory, 'reduced_cubes'))


# %%
# Slice cubes and increment timestamps
redu_cubelist = sorted(glob1(directory+'/reduced_cubes', '*.fits'))
for cubename in redu_cubelist:
    slice_cube(directory=os.path.join(directory, 'reduced_cubes'), 
                cubename=cubename,
                out_dir=os.path.join(directory, 'sliced'))


# %%
# Alipy source detection on all images
ref_image, identifications = alipy_ident(directory=os.path.join(directory, 'sliced'))

# %%
# Align all images to ref img, preserve headers
alipy_align(ref_image, identifications, out_dir=os.path.join(directory, 'aligned'))

# %%
# Solve frame and migrate WCS headers
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

# %%
