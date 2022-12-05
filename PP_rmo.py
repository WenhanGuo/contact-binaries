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

directory = '/Volumes/TMO_Data_4TB/cb_data/C01+13/20221121/aligned'

cubelist = sorted(glob1(directory, '*.fits'))

out_dir = '/Volumes/TMO_Data_4TB/cb_data/C01+13/20221121'



# %%
hdulist = fits.open(os.path.join(directory, cubelist[0]))
hdu = hdulist[0]
# img = hdu.data[0]
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
reflist = [1, 2]

diff_lc_name = os.path.join(out_dir, 'diff_lc.ecsv')

diff_lc = differential_photometry(directory=out_dir, reflist=reflist)
diff_lc.write(diff_lc_name, overwrite=True)


# %%
obj_dir = '/Users/danny/Mirror/ASTRO/JPL_NEO/Contact_Binary/data/CSS_034852'

fold_lc(obj_dir=obj_dir, table='diff_lc.ecsv')


# %%
import seaborn as sns
sns.set()

# to be beautified
def diff_phot(out_dir, refnum):
    ts = TimeSeries.read(os.path.join(out_dir, 'norm_target_flux.ecsv'), time_column='time')
    ref = TimeSeries.read(os.path.join(out_dir, 'norm_ref_flux.ecsv'), time_column='time')

    ts['target_mag'] = -2.5 * np.log10(ts['target_flux']) + 25
    ts['target_flux_err'] = np.sqrt(ts['target_flux'])
    ts['target_mag_err'] = 2.512 * ts['target_flux_err'] / (ts['target_flux'] * np.log(10))

    ref['ref_mag'] = -2.5 * np.log10(ref['ref_flux']) + 25
    ref['ref_flux_err'] = np.sqrt(ref['ref_flux'])
    ref['ref_mag_err'] = 2.512 * ref['ref_flux_err'] / (ref['ref_flux'] * np.log(10))

    g_ts = ts[ts['filter'] == 'g']
    g_ref = ref[ref['filter'] == 'g']

    fig, axes = plt.subplots(1, 2, figsize=(20,12))
    axes[0].scatter(g_ts.time.datetime64, g_ts['target_mag'], label='target g')
    axes[0].scatter(g_ref.time.datetime64, g_ref['ref_mag'][:,refnum-1], label='ref1')

    g_diff = g_ts['target_mag'] - (g_ref['ref_mag'][:,refnum-1])
    # g_diff = g_diff - np.mean(g_diff) .  # normalization is what causes problems in folding
    g_ts['diff_mag'] = g_diff
    g_ts.write(out_dir+'/diff_lc.ecsv', overwrite=True)

    axes[1].scatter(g_ts.time.datetime64, g_diff)
    axes[1].set_ylim(2.9, 2.3)

    axes[0].legend()
    axes[0].legend()
    axes[0].legend()
    axes[0].legend()
    plt.tight_layout()
    plt.savefig(out_dir+'/test.pdf')

    return

diff_phot('/Volumes/TMO_Data_4TB/cb_data/C01+13/20221106', 5)
diff_phot('/Volumes/TMO_Data_4TB/cb_data/C01+13/20221120', 4)
diff_phot('/Volumes/TMO_Data_4TB/cb_data/C01+13/20221121', 5)

# %%
# del ts[1448]
# del ts[1449]
ts['diff_mag'] = g_diff
ts['diff_mag'] = ts['diff_mag'] - np.mean(ts['diff_mag'])
ts_folded = ts.fold(period=0.3439788 * u.day)
del ts_folded['filter']

from astropy.timeseries import aggregate_downsample
ts_binned = aggregate_downsample(ts_folded, time_bin_size=60 * u.s)
fig, ax = plt.subplots(1,1, figsize=(12,10))
plt.scatter(ts_folded.time.jd, ts_folded['diff_mag'])
plt.plot(ts_binned.time_bin_start.jd, ts_binned['diff_mag'], 'r-', drawstyle='steps-post')
plt.ylim(0.3,-0.3)
plt.savefig(out_dir+'/60bin.pdf')

# %%
obj_dir = '/Volumes/TMO_Data_4TB/cb_data/C01+13'
date_folders = glob1(obj_dir, '[0-9]*')
# date_folders = ['20221106', '20221121']
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
colors = ['tab:red', 'tab:green', 'tab:blue', 'k']
joined_ts['color'] = 'placeholder'
for row in joined_ts:
    for i in range(len(dates)):
        if row['obs_no'] == dates[i]:
            row['color'] = colors[i]

joined_ts.write(obj_dir+'/joined_ts_3nights.csv', overwrite=True)
folded_ts = joined_ts.fold(period=0.3439788 * u.day)

plt.scatter(folded_ts.time.jd, folded_ts['diff_mag'], c=folded_ts['color'], 
                        s=10, marker='x', alpha=1, linewidth=0.5)

plt.ylim(0.9,0.35)
plt.savefig(obj_dir+'/test.pdf')


# %%
def diff_phot_ref(out_dir, refnum1, refnum2):
    ref = TimeSeries.read(os.path.join(out_dir, 'norm_ref_flux.ecsv'), time_column='time')

    ref['target_mag'] = -2.5 * np.log10(ref['ref_flux'][:,refnum1-1]) + 25
    ref['ref_mag'] = -2.5 * np.log10(ref['ref_flux'][:,refnum2-1]) + 25

    g_ref = ref[ref['filter'] == 'g']

    fig, axes = plt.subplots(1, 2, figsize=(20,12))
    axes[0].scatter(g_ref.time.datetime64, g_ref['target_mag'], label='target')
    axes[0].scatter(g_ref.time.datetime64, g_ref['ref_mag'], label='ref1')

    g_diff = (g_ref['target_mag']) - (g_ref['ref_mag'])
    # g_diff = g_diff - np.mean(g_diff) .  # normalization is what causes problems in folding
    g_ref['diff_mag'] = g_diff
    g_ref.write(out_dir+'/diff_lc_ref.ecsv', overwrite=True)

    axes[1].scatter(g_ref.time.datetime64, g_diff)
    # axes[1].set_ylim(2.9, 2.3)

    axes[0].legend()
    axes[0].legend()
    axes[0].legend()
    axes[0].legend()
    plt.tight_layout()
    plt.savefig(out_dir+'/test_ref.pdf')

    return

diff_phot_ref('/Volumes/TMO_Data_4TB/cb_data/C01+13/20221106', 3, 5)
diff_phot_ref('/Volumes/TMO_Data_4TB/cb_data/C01+13/20221120', 2, 4)
diff_phot_ref('/Volumes/TMO_Data_4TB/cb_data/C01+13/20221121', 2, 5)

# %%
obj_dir = '/Volumes/TMO_Data_4TB/cb_data/C01+13'
date_folders = glob1(obj_dir, '[0-9]*')
# date_folders = ['20221106', '20221121']
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
colors = ['tab:red', 'tab:green', 'tab:blue', 'k']
joined_ts['color'] = 'placeholder'
for row in joined_ts:
    for i in range(len(dates)):
        if row['obs_no'] == dates[i]:
            row['color'] = colors[i]

joined_ts.write(obj_dir+'/joined_ts_3nights_ref.csv', overwrite=True)
folded_ts = joined_ts.fold(period=0.3439788 * u.day)

plt.scatter(folded_ts.time.jd, folded_ts['diff_mag'], c=folded_ts['color'], 
                        s=10, marker='x', alpha=1, linewidth=0.5)

plt.ylim(-1.85,-2.1)
plt.savefig(obj_dir+'/test_ref.pdf')