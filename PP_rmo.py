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
from astropy.timeseries import TimeSeries, aggregate_downsample
from astropy.table import vstack
import astropy.units as u
from pycaret.anomaly import *

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import seaborn as sns
sns.set()

directory = '/Volumes/TMO_Data_4TB/cb_data/C01+13/20221121/aligned'
cubelist = sorted(glob1(directory, '*.fits'))
out_dir = '/Volumes/TMO_Data_4TB/cb_data/C01+13/20221121'

obj_dir = '/Volumes/TMO_Data_4TB/cb_data/C01+13'
save_dir = '/Users/TMObserver/Documents/cb_scripts'

# %%
hdulist = fits.open(os.path.join(directory, cubelist[0]))
hdu = hdulist[0]
img = hdu.data
wcs = WCS(hdu.header, naxis=2)

RA_Dec = '01h18m48.5395s +13d21m07.62s'
target_sky_aperture, target_sky_annulus = sky_aperture_from_RADec(wcs=wcs, RA_Dec=RA_Dec, 
                                    radius=10.0, annulus=True, r_in=30.0, r_out=45.0)

ref_sky_aperture, ref_sky_annulus = detect_stars(img, wcs, skymap=True, detection_threshold=5, 
                                    radius=10.0, annulus=True, r_in=30.0, r_out=45.0, 
                                    xbounds=[100,1200], ybounds=[150,1200], 
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
table, RFM = LCs_visualizer(directory=out_dir, visu_dir=out_dir, 
                                    mode='full', layout=[3,5])

# %%
# reflist = [1, 2]
# diff_lc_name = os.path.join(out_dir, 'diff_lc.ecsv')

# diff_lc = differential_photometry(directory=out_dir, reflist=reflist)
# diff_lc.write(diff_lc_name, overwrite=True)

# %%
# obj_dir = '/Users/danny/Mirror/ASTRO/JPL_NEO/Contact_Binary/data/CSS_034852'
# fold_lc(obj_dir=obj_dir, table='diff_lc.ecsv')

# %%
# Differential photometry on target
diff_phot_target('/Volumes/TMO_Data_4TB/cb_data/C01+13/20221106', 5)
diff_phot_target('/Volumes/TMO_Data_4TB/cb_data/C01+13/20221120', 4)
diff_phot_target('/Volumes/TMO_Data_4TB/cb_data/C01+13/20221121', 5)

# %%
date_folders = glob1(obj_dir, '[0-9]*')
joined_ts = 1.0   # init value placeholder
table = 'diff_lc.ecsv'

for date in date_folders:
    ts = TimeSeries.read(os.path.join(obj_dir, date, table), time_column='time')
    ts['obs_no'] = date
    if type(joined_ts) == float:
        joined_ts = ts
    else:
        joined_ts = vstack([joined_ts, ts])

dates = date_folders
colors = ['crimson', 'tab:blue', 'k']
joined_ts['color'] = 'placeholder'
for row in joined_ts:
    for i in range(len(dates)):
        if row['obs_no'] == dates[i]:
            row['color'] = colors[i]

joined_ts.write(obj_dir+'/ts_3nights.csv', overwrite=True)
folded_ts = joined_ts.fold(period=0.3439788 * u.day)

plt.scatter(folded_ts.time.jd, folded_ts['diff_mag'], c=folded_ts['color'], 
                        s=10, marker='x', alpha=1, linewidth=0.5)
plt.xlabel('phase')
plt.ylabel('dmag')
plt.ylim(0.9,0.3)
plt.savefig(save_dir+'/lc_r10.pdf', bbox_inches='tight', overwrite=True)

# %%
# Outlier rejection and lc binning
url = 'https://raw.githubusercontent.com/WenhanGuo/contact-binaries/master/ts_3nights.csv'
url_out = 'https://raw.githubusercontent.com/WenhanGuo/contact-binaries/master/outliers.csv'

ts = reject_outliers(save_dir, url, 0.3)
bin_lc(save_dir, ts)

# %%
# Differential photometry on a reference star of choice
diff_phot_ref('/Volumes/TMO_Data_4TB/cb_data/C01+13/20221106', 5, 2)
diff_phot_ref('/Volumes/TMO_Data_4TB/cb_data/C01+13/20221120', 4, 1)
diff_phot_ref('/Volumes/TMO_Data_4TB/cb_data/C01+13/20221121', 5, 1)

# %%
date_folders = glob1(obj_dir, '[0-9]*')
joined_ts = 1.0   # init value placeholder
table = 'diff_lc_ref.ecsv'

for date in date_folders:
    ts = TimeSeries.read(os.path.join(obj_dir, date, table), time_column='time')
    ts['obs_no'] = date
    ts = ts['time', 'target_mag', 'ref_mag', 'diff_mag', 'obs_no']
    if type(joined_ts) == float:
        joined_ts = ts
    else:
        joined_ts = vstack([joined_ts, ts])

dates = date_folders
colors = ['k', 'crimson', 'tab:blue']
joined_ts['color'] = 'placeholder'
for row in joined_ts:
    for i in range(len(dates)):
        if row['obs_no'] == dates[i]:
            row['color'] = colors[i]

joined_ts.write(obj_dir+'/ts_3nights_ref.csv', overwrite=True)
folded_ts = joined_ts.fold(period=0.3439788 * u.day)

plt.scatter(folded_ts.time.jd, folded_ts['diff_mag'], c=folded_ts['color'], 
                        s=10, marker='x', alpha=1, linewidth=0.5)
plt.xlabel('phase')
plt.ylabel('diff mag')
plt.ylim(-1.85,-2.1)
plt.savefig(obj_dir+'/vari_ref.pdf')