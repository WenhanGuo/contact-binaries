# %%
import os
import numpy as np
from glob import glob1

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.stats import SigmaClip, sigma_clipped_stats, gaussian_fwhm_to_sigma
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.timeseries import TimeSeries, aggregate_downsample
from astropy.table import vstack, Table
import astropy.units as u

from skimage.draw import disk
from sklearn.preprocessing import normalize as sknorm
from photutils.background import Background2D, SExtractorBackground
from photutils.segmentation import detect_sources, deblend_sources, SourceCatalog
from photutils.aperture import CircularAperture, CircularAnnulus, ApertureStats
from scipy.signal import savgol_filter
from pycaret.anomaly import *

import pandas as pd
import numpy as np

from visu_functions import *



def estimate_background(img, mode='2D', visu=False, visu_dir=None, high_res=False):
    """
    Background estimation for an image using photutils. 

    mode=='2D': 2D background subtraction
    mode=='1D': uniform sigma-clip median background estimation
    visu==True: visualize estimated skymap
    """
    if mode != '1D' and mode != '2D':
        print('mode keyword only accept \'1D\' or \'2D\'.')
    elif mode == '2D':
        # 2D background estimation using box size = 100,100
        sigma_clip = SigmaClip(sigma=3.)
        bkg_estimator = SExtractorBackground()
        bkg = Background2D(img, box_size=(100, 100), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    elif mode == '1D':
        # uniform sigma-clip median background estimation
        _, median, _ = sigma_clipped_stats(img, sigma=3.0)
        return median
    
    # visualize estimated skymap
    if visu == True:
        visualize_background(img, skymap=bkg.background, visu_dir=visu_dir, high_res=high_res)

    return bkg.background



def detect_stars(img, wcs, skymap=True, detection_threshold=5, radius=None, annulus=False, r_in=None, r_out=None, 
                                                        xbounds=None, ybounds=None, mask_coord=None, 
                                                        visu=False, visu_dir=None, high_res=False):
    """
    Detect sources within a bounding box, while masking target RA-Dec.
    Return a list of centroid RA-Decs of detected sources.

    img: 2D numpy data array
    wcs: WCS header of cube
    detection_threshold: number of std above background to qualify as a star
    radius: aperture photometry radius
    r_in: annulus photometry inner radius
    r_out: annulus photometry outer radius
    xbounds: horizontal bounding box [xmin, xmax] pixel positions
    ybounds: vertical bounding box [ymin, ymax] pixel positions
    mask_coord: target RA-Dec position to create circular aperture mask of r=30; exclude from detection
    visu==True: visualize source detection and skymap estimation
    visu_dir: directory to save visualization results
    high_res: if want to save high resolution results; takes longer time
    """
    if not radius:
        radius = 15

    # subtract 2D background
    # ----------------------
    if skymap == True:
        skymap = estimate_background(img, mode='2D', visu=False)
        img -= skymap

    # create bounding box mask and source mask
    # ----------------------------------------
    assert type(xbounds) == type(ybounds)
    # create bounding box mask
    if xbounds and ybounds:
        # create mask=True for all pixel, then mask=False for region within bounding box
        img_mask = ~np.zeros(img.shape, dtype=bool)
        img_mask[ybounds[0]:ybounds[1]+1, xbounds[0]:xbounds[1]+1] = False
    else:
        # create null img_mask (img_mask=False for every pixel)
        img_mask = np.zeros(img.shape, dtype=bool)
    # create masked image indicating bounding box
    masked_img = np.ma.array(img, mask=img_mask)

    # create source mask on target coord
    if mask_coord:
        # circular aperture mask of r=30 pix around target RA-Dec
        # create mask=False for all pixel, then mask=True for target aperture
        # https://scikit-image.org/docs/stable/api/skimage.draw.html?highlight=disk#skimage.draw.disk
        pixel_coord = wcs.world_to_pixel(SkyCoord(mask_coord))
        source_mask = np.zeros(img.shape, dtype=bool)
        rr, cc = disk(center=(pixel_coord[1], pixel_coord[0]), radius=30)
        source_mask[rr, cc] = True
    else:
        # create null source_mask (source_mask=False for every pixel)
        source_mask = np.zeros(img.shape, dtype=bool)
        pixel_coord = None

    # combine bounding box mask (img_mask) and source mask to form master mask
    master_mask = img_mask + source_mask

    # source detection with photutils image segmentation
    # --------------------------------------------------
    # define detection threshold
    _, _, std = sigma_clipped_stats(img, sigma=3.0)
    threshold = detection_threshold * std

    # convolve image prior to segmentation
    sigma = 3.0 * gaussian_fwhm_to_sigma  # FWHM = 3.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    convolved_img = convolve(img, kernel)

    # npixels: how many connected pixels, each above threshold, should an area have to qualify as a source
    npixels = 10
    segment_map = detect_sources(convolved_img, threshold, npixels=npixels, mask=master_mask)

    # generate greyzone_mask: 10 pixels zone within bounding box
    greyzone_mask = ~np.zeros(img.shape, dtype=bool)
    greyzone_mask[ybounds[0]+10:ybounds[1]-9, xbounds[0]+10:xbounds[1]-9] = False
    # remove detected sources in greyzone_mask: edge effect on mask boundary causes centroid inaccuracies
    segment_map.remove_masked_labels(greyzone_mask, partial_overlap=True, relabel=True)

    # deblend sources
    segm_deblend = deblend_sources(convolved_img, segment_map, npixels=npixels, 
                                    nlevels=32, contrast=0.001)
    # photutils SourceCatalog table for source detection
    cat = SourceCatalog(img, segm_deblend, convolved_data=convolved_img)
    columns = ['label','xcentroid','ycentroid','fwhm','gini','eccentricity','orientation','kron_flux']
    table = cat.to_table(columns=columns)

    # format table
    table['xcentroid'].info.format = '.1f'
    table['ycentroid'].info.format = '.1f'
    table['fwhm'].info.format = '.2f'
    table['gini'].info.format = '.4f'
    table['eccentricity'].info.format = '.2f'
    table['orientation'].info.format = '.1f'
    table['kron_flux'].info.format = '.0f'
    print(table)

    # visualization
    # ---------------------
    if visu == True:
        visualize_detection(img=img, masked_img=masked_img, skymap=skymap, table=table, 
                            radius=radius, r_in=r_in, r_out=r_out, 
                            xbounds=xbounds, ybounds=ybounds, pixel_coord=pixel_coord, 
                            visu_dir=visu_dir, high_res=high_res)

    # create sky_aperture object and sky_annulus object
    # -------------------------------------------------
    positions = []
    for n in range(len(table)):
        positions.append((table[n]['xcentroid'], table[n]['ycentroid']))
        
    if annulus == False:
        sky_aperture = sky_aperture_from_pix(wcs=wcs, positions=positions, radius=radius, 
                            annulus=False, r_in=r_in, r_out=r_out)
        return sky_aperture

    elif annulus == True:
        sky_aperture, sky_annulus = sky_aperture_from_pix(wcs=wcs, positions=positions, radius=radius, 
                            annulus=True, r_in=r_in, r_out=r_out)
        return sky_aperture, sky_annulus

                        

def sky_aperture_from_pix(wcs, positions, radius, annulus=False, r_in=None, r_out=None):
    aperture = CircularAperture(positions, r=radius)
    sky_aperture = aperture.to_sky(wcs)
    if annulus == True:
        assert type(r_in) == float and type(r_out) == float
        annulus_aperture = CircularAnnulus(positions, r_in=r_in, r_out=r_out)
        sky_annulus = annulus_aperture.to_sky(wcs)
        return sky_aperture, sky_annulus
    elif annulus == False:
        return sky_aperture



def sky_aperture_from_RADec(wcs, RA_Dec, radius, annulus=True, r_in=None, r_out=None):
    # convert RA-Dec into xy pixel position for target star
    pixel_coord = wcs.world_to_pixel(SkyCoord(RA_Dec))
    position = [(pixel_coord[0], pixel_coord[1])]
    # define aperture object
    aperture = CircularAperture(position, r=radius)
    sky_aperture = aperture.to_sky(wcs)
    if annulus == True:
        assert type(r_in) == float and type(r_out) == float
        # define annulus aperture object
        annulus_aperture = CircularAnnulus(position, r_in=r_in, r_out=r_out)
        sky_annulus = annulus_aperture.to_sky(wcs)
        return sky_aperture, sky_annulus
    elif annulus == False:
        return sky_aperture



def simple_photometry(directory, cubename, target_coord, radius, annulus=False, r_in=None, r_out=None):
    """
    Simple time-resolved circular aperture/annulus photometry on one RA-Dec. 
    """
    # open cube and obtain wcs info
    hdulist = fits.open(os.path.join(directory, cubename))
    hdu = hdulist[0]
    nframes = int(hdu.header['NAXIS3'])
    exptime = float(hdu.header['EXPTIME']) * u.s
    dateobs = hdu.header['DATE-OBS'][:-6]
    wcs = WCS(hdu.header, naxis=2)

    # convert RA-Dec into xy pixel position for target star
    pixel_coord = wcs.world_to_pixel(SkyCoord(target_coord))
    position = [(pixel_coord[0], pixel_coord[1])]
    print('Target star pixel position =', position)

    # define aperture object
    aperture = CircularAperture(position, r=radius)
    if annulus == True:
        # define annulus aperture object
        annulus_aperture = CircularAnnulus(position, r_in=r_in, r_out=r_out)
        sigclip = SigmaClip(sigma=3.0, maxiters=10)

    # aperture photometry for frames in cube
    fluxlist = []
    for n in range(nframes):
        # subtract 2D background
        img = hdu.data[n]
        skymap = estimate_background(img, mode='2D', visu=False)
        img -= skymap

        # obtain statistics for each aperture
        aper_stats = ApertureStats(img, aperture, sigma_clip=None)  

        if annulus == True:
            assert type(r_in) == float and type(r_out) == float
            # sigma-clipped median within a circular annulus
            bkg_stats = ApertureStats(img, annulus_aperture, sigma_clip=sigclip)
            total_bkg = bkg_stats.median * aper_stats.sum_aper_area.value
            apersum_bkgsub = aper_stats.sum - total_bkg
            fluxlist.append(apersum_bkgsub[0])
            # print('median, total_bkg =', bkg_stats.median, total_bkg)
        else:
            fluxlist.append(aper_stats.sum[0])
    
    # output astropy TimeSeries object
    ts = TimeSeries(time_start=dateobs, time_delta=exptime, n_samples=nframes)
    ts['target_flux'] = fluxlist

    return ts



def parallel_photometry(directory, cubename, sky_aperture, annulus=False, sky_annulus=None):
    """
    Time-resolved circular aperture/annulus photometry on multiple RA-Dec coordinates.
    Must provide a sky_aperture object containing RA-Decs. 
    """
    # open cube and obtain wcs info
    hdulist = fits.open(os.path.join(directory, cubename))
    hdu = hdulist[0]
    nframes = int(hdu.header['NAXIS3'])
    exptime = float(hdu.header['EXPTIME']) * u.s
    dateobs = hdu.header['DATE-OBS'][:-6]
    wcs = WCS(hdu.header, naxis=2)

    # define aperture object by converting from sky to pixel
    aperture = sky_aperture.to_pixel(wcs)
    if annulus == True:
        assert sky_annulus != None
        # convert sky annulus to pixel
        annulus_aperture = sky_annulus.to_pixel(wcs)
        sigclip = SigmaClip(sigma=3.0, maxiters=10)
    print('aperture =', aperture)
    print('annulus_aperture =', annulus_aperture)

    # aperture photometry for frames in cube
    fluxlist = []
    for n in range(nframes):
        # subtract 2D background
        img = hdu.data[n]
        skymap = estimate_background(img, mode='2D', visu=False)
        img -= skymap

        # obtain statistics for each aperture
        aper_stats = ApertureStats(img, aperture, sigma_clip=None)  

        if annulus == True:
            # sigma-clipped median within a circular annulus
            bkg_stats = ApertureStats(img, annulus_aperture, sigma_clip=sigclip)
            total_bkg = bkg_stats.median * aper_stats.sum_aper_area.value
            apersum_bkgsub = aper_stats.sum - total_bkg
            fluxlist.append(apersum_bkgsub)
        else:
            fluxlist.append(aper_stats.sum)
    
    # output astropy TimeSeries object
    ts = TimeSeries(time_start=dateobs, time_delta=exptime, n_samples=nframes)
    ts['ref_flux'] = fluxlist
     
    return ts



def dual_thread_photometry(directory, cubename, target_sky_aperture, ref_sky_aperture,
                            annulus=False, target_sky_annulus=None, ref_sky_annulus=None):
    """
    Combining simple_photometry and parallel_photometry. 
    Time-resolved circular aperture/annulus photometry on target star and reference stars,
    Avoid opening FITS cubes twice.
    Must provide target RA-Dec and ref star sky_aperture object.
    """
    # open cube and obtain wcs info
    hdulist = fits.open(os.path.join(directory, cubename))
    hdu = hdulist[0]
    if 'NAXIS3' in hdu.header:
        # the image is a 3D cube
        nframes = int(hdu.header['NAXIS3'])
    else:
        # the image is a 2D slice
        nframes = 1
    exptime = float(hdu.header['EXPTIME']) * u.s
    dateobs = hdu.header['DATE-OBS'][:-6]
    filt = hdu.header['FILTER']
    wcs = WCS(hdu.header, naxis=2)

    # define target aperture object by converting from sky to pixel
    target_aperture = target_sky_aperture.to_pixel(wcs)
    ref_aperture = ref_sky_aperture.to_pixel(wcs)
    if annulus == True:
        assert target_sky_annulus != None and ref_sky_annulus != None
        # convert sky annulus to pixel
        target_annulus_aperture = target_sky_annulus.to_pixel(wcs)
        ref_annulus_aperture = ref_sky_annulus.to_pixel(wcs)
        sigclip = SigmaClip(sigma=3.0, maxiters=10)
    print('target aperture =', target_aperture)
    print('target annulus =', target_annulus_aperture)
    print('ref star aperture len =', len(ref_aperture))
    print('ref star annulus len =', len(ref_annulus_aperture))

    # aperture photometry for frames in cube
    target_fluxlist = []
    ref_fluxlist = []
    for n in range(nframes):
        # subtract 2D background
        if nframes >= 2:
            img = hdu.data[n]
        if nframes == 1:
            img = hdu.data
        skymap = estimate_background(img, mode='2D', visu=False)
        img -= skymap

        # obtain statistics for each aperture
        target_aper_stats = ApertureStats(img, target_aperture, sigma_clip=None)
        ref_aper_stats = ApertureStats(img, ref_aperture, sigma_clip=None)

        if annulus == True:
            # sigma-clipped median within a circular annulus
            target_bkg_stats = ApertureStats(img, target_annulus_aperture, sigma_clip=sigclip)
            target_total_bkg = target_bkg_stats.median * target_aper_stats.sum_aper_area.value
            target_apersum_bkgsub = target_aper_stats.sum - target_total_bkg
            target_fluxlist.append(target_apersum_bkgsub[0])

            ref_bkg_stats = ApertureStats(img, ref_annulus_aperture, sigma_clip=sigclip)
            ref_total_bkg = ref_bkg_stats.median * ref_aper_stats.sum_aper_area.value
            ref_apersum_bkgsub = ref_aper_stats.sum - ref_total_bkg
            ref_fluxlist.append(ref_apersum_bkgsub)
        else:
            target_fluxlist.append(target_aper_stats.sum[0])
            ref_fluxlist.append(ref_aper_stats.sum)

    # output astropy TimeSeries object
    target_ts = TimeSeries(time_start=dateobs, time_delta=exptime, n_samples=nframes)
    target_ts['target_flux'] = target_fluxlist
    target_ts['filter'] = filt

    ref_ts = TimeSeries(time_start=dateobs, time_delta=exptime, n_samples=nframes)
    ref_ts['ref_flux'] = ref_fluxlist
    ref_ts['filter'] = filt
    
    return target_ts, ref_ts



def update_table(ts, mother_tablename):
    """
    Utility function for updating (concatenating) TimeSeries tables.
    """
    if not os.path.exists(mother_tablename):
        # write ts table to ecsv for the first time
        ts.write(mother_tablename, overwrite=False)
    else:
        # concatenate ts table onto existing mother table
        existing_ts = TimeSeries.read(mother_tablename, time_column='time')
        new_ts = vstack([existing_ts, ts])
        new_ts.write(mother_tablename, overwrite=True)

    return



def normalize_target_table(directory, target_table='target_flux.ecsv'):
    """
    Normalize target flux in each filter separately.
    """
    table = Table.read(os.path.join(directory, target_table))
    table_by_filt = table.group_by('filter')
    # number of filters = number of group key names
    nfilt = len(table_by_filt.groups.keys)

    norm_table = None
    for n in range(nfilt):
        # get a subtable containing 1 specific filter
        subtable = table_by_filt.groups[n]
        # normalize target_flux column
        tf = subtable['target_flux']
        subtable['norm_target_flux'] = tf / np.mean(tf)
        # concatenate normalized subtables to form new norm_table
        if norm_table is None:
            norm_table = subtable
        else:
            norm_table = vstack([norm_table, subtable])
    norm_table.write(os.path.join(directory, 'norm_target_flux.ecsv'), overwrite=True)

    return



def normalize_ref_table(directory, ref_table='ref_flux.ecsv'):
    """
    Normalize reference flux in each filter separately.
    The column 'norm_ref_flux' is a 2d numpy array, same shape as RFM.
    x axis of norm_ref_flux: different reference stars
    y axis of norm_ref_flux: flux for 1 ref star, normalized separately by filters
    For an explanation of RFM please see LCs_visualizer()
    """
    table = Table.read(os.path.join(directory, ref_table))
    table_by_filt = table.group_by('filter')
    # number of filters = number of group key names
    nfilt = len(table_by_filt.groups.keys)

    norm_table = None
    for n in range(nfilt):
        # get a subtable containing 1 specific filter
        subtable = table_by_filt.groups[n]
        # normalize ref_flux column (2d)
        RFM = subtable['ref_flux']
        subtable['norm_ref_flux'] = RFM
        for nref in range(RFM.shape[1]):
            subtable['norm_ref_flux'][:,nref] = RFM[:,nref] / np.mean(RFM[:,nref])
        # concatenate normalized subtables to form new norm_table
        if norm_table is None:
            norm_table = subtable
        else:
            norm_table = vstack([norm_table, subtable])
    norm_table.write(os.path.join(directory, 'norm_ref_flux.ecsv'), overwrite=True)

    return



def LCs_visualizer(directory, visu_dir, mode='simple', ref_flux='ref_flux.ecsv', 
                            ylim=[0.85, 1.15], layout=[4, 6]):
    """
    Visualize lightcurves for reference stars. 
    """
    assert mode in ['simple', 'full', 'raw']

    table = TimeSeries.read(os.path.join(directory, ref_flux), time_column='time')
    # each row of table['ref_flux'] column contains a list of reference flux
    # table['ref_flux'] is thus a 2d numpy array; we can treat it as an image in all cases
    # from now on, table['ref_flux'] is refered to as 'reference flux map' (RFM)
    # x axis of RFM = number of reference stars
    # y axis of RFM = time-dependent flux
    RFM = table['ref_flux']

    # 1D RFM visualization
    if mode in ['simple', 'full']:
        visualize_1D_RFM(table=table, mode='norm', visu_dir=visu_dir, ylim=ylim, layout=layout)
    elif mode == 'raw':
        visualize_1D_RFM(table=table, mode='raw', visu_dir=visu_dir, ylim=ylim, layout=layout)

    # 2D RFM visualization
    visualize_2D_RFM(table=table, visu_dir=visu_dir, ylim=ylim)

    # 3D RFM visualization
    if mode in ['simple', 'raw']:
        visualize_3D_RFM(table=table, mode='simple', visu_dir=visu_dir, ylim=ylim)
    elif mode == 'full':
        visualize_3D_RFM(table=table, mode='full scatter', visu_dir=visu_dir, ylim=ylim)

    return table, RFM



def differential_photometry(directory, reflist, target_flux='target_flux.ecsv', ref_flux='ref_flux.ecsv'):
    """
    Differential photometry on target star by selected reference stars. 
    Output diff_lc.
    """
    # read target_flux, ref_flux; initialize diff_lc table
    target_table = TimeSeries.read(os.path.join(directory, target_flux), time_column='time')
    ref_table = TimeSeries.read(os.path.join(directory, ref_flux), time_column='time')
    diff_lc = target_table

    # get specific columns from ref_table['ref_flux'] according to reflist
    raw_ref_flux = ref_table['ref_flux'][:, np.array(reflist)-1]
    diff_lc['raw_ref_flux'] = raw_ref_flux
    # normalize raw ref flux by column (by each ref star)
    norm_ref_flux = sknorm(raw_ref_flux, axis=0)
    norm_ref_flux = norm_ref_flux / np.mean(norm_ref_flux)
    diff_lc['norm_ref_flux'] = norm_ref_flux

    # smooth norm ref flux by column
    smooth_ref_flux = np.empty((0, len(diff_lc)))
    for n in range(len(reflist)):
        norm_nref_flux = diff_lc['norm_ref_flux'][:, n]
        smooth_nref_flux = savgol_filter(norm_nref_flux, window_length=51, polyorder=1)
        smooth_ref_flux = np.append(smooth_ref_flux, np.array([smooth_nref_flux]), axis=0)
    diff_lc['smooth_ref_flux'] = smooth_ref_flux.transpose()

    # mean ref flux column
    diff_lc['mean_ref_flux'] = np.mean(diff_lc['smooth_ref_flux'], axis=1)

    # differential photometry
    diff_lc['diff_flux'] = diff_lc['target_flux'] / diff_lc['mean_ref_flux']
    diff_lc['norm_diff_flux'] = diff_lc['diff_flux'] / np.mean(diff_lc['diff_flux'])

    return diff_lc



def fold_lc(obj_dir, table='diff_lc.ecsv'):
    """
    Visualize lightcurve from differential photometry. 
    """
    date_folders = glob1(obj_dir, '[0-9]*')
    joined_ts = None
    for date in date_folders:
        ts = TimeSeries.read(os.path.join(obj_dir, date, table), time_column='time')
        # plt.scatter(ts.time.datetime64, ts['diff_flux'], s=2)
        if joined_ts == None:
            joined_ts = ts
        else:
            joined_ts = vstack([joined_ts, ts])
    folded_ts = joined_ts.fold(period=0.2053320 * u.day)
    plt.scatter(folded_ts.time.jd, folded_ts['norm_diff_flux'], 
                            s=10, marker='x', alpha=1, linewidth=0.5)
    return joined_ts



def reject_outliers(directory, table, fraction):
    """
    Outlier detection and rejection.
    table: differential photometry csv file with 'time' column of format='yyyy-mm-ddTHH:MM:SS.000'.
    fraction: fraction of data set that wants to be identified as outliers.
    """
    df = pd.read_csv(table)
    df.set_index(pd.DatetimeIndex(df['time']), inplace=True)
    df.drop(['time'], axis=1)

    df['diff_mag'] = np.float64(df['diff_mag'])
    df['avg_mag'] = df['diff_mag'].rolling(6).mean()
    
    df.drop('avg_mag', axis=1, inplace=True)
    
    s = setup(df, session_id=123)
    iforest = create_model('iforest', fraction=fraction)
    iforest_results = assign_model(iforest)
    
    outliers = iforest_results[iforest_results['Anomaly'] == 1]
    out = outliers.to_csv(directory+'outliers.csv')

    df_out = pd.read_csv(os.path.join(directory, 'outliers.csv'))
    df_out.set_index(pd.DatetimeIndex(df_out['time']), inplace=True)
    df_out = df_out.drop(columns=['time', 'Anomaly', 'Anomaly_Score'])
    df_out = df_out.loc[:, ~df_out.columns.str.contains('^Unnamed')]
    
    df = pd.concat([df, df_out])
    df = df[~df.index.duplicated(keep=False)]
    df.set_index(pd.DatetimeIndex(df['time']), inplace=True)
    df = df.drop(['time'], axis=1)
    ts = TimeSeries.from_pandas(df)
    
    return ts



def bin_lc(directory, ts):
    """
    Fold and bin light curve.
    """
    ts_folded = ts.fold(period=0.3439788 * u.day)
    del ts_folded['filter', 'color']

    ts_binned = aggregate_downsample(ts_folded, time_bin_size=60*u.s)

    plt.scatter(ts_folded.time.jd, ts_folded['diff_mag'], c='k', s=5)
    plt.scatter(ts_binned.time_bin_start.jd, ts_binned['diff_mag'], c='crimson', s=3, marker='x')
    plt.xlabel('phase')
    plt.ylabel('normalized dmag')
    plt.ylim(0.9, 0.3)
    plt.savefig(directory+'/lc_binned.pdf', bbox_inches='tight', overwrite=True)

    return



def diff_phot_target(directory, refnum):
    """
    Differential photometry on target.
    """
    ts = TimeSeries.read(os.path.join(directory, 'norm_target_flux.ecsv'), time_column='time')
    ref = TimeSeries.read(os.path.join(directory, 'norm_ref_flux.ecsv'), time_column='time')

    ts['target_mag'] = -2.5 * np.log10(ts['target_flux']) + 25
    ts['target_flux_err'] = np.sqrt(ts['target_flux'])
    ts['target_mag_err'] = 2.512 * ts['target_flux_err'] / (ts['target_flux'] * np.log(10))

    ref['ref_mag'] = -2.5 * np.log10(ref['ref_flux']) + 25
    ref['ref_flux_err'] = np.sqrt(ref['ref_flux'])
    ref['ref_mag_err'] = 2.512 * ref['ref_flux_err'] / (ref['ref_flux'] * np.log(10))

    g_ts = ts[ts['filter'] == 'g']
    g_ref = ref[ref['filter'] == 'g']

    fig, axes = plt.subplots(2, 1, figsize=(20,12))
    axes[0].scatter(g_ts.time.datetime64, g_ts['target_mag'], c='k', s=5, label='target mag')
    axes[0].scatter(g_ref.time.datetime64, g_ref['ref_mag'][:,refnum-1], c='crimson', s=5, label='ref mag')

    g_ts['diff_mag'] = g_ts['target_mag'] - (g_ref['ref_mag'][:,refnum-1])
    g_ts.write(directory+'/diff_lc.ecsv', overwrite=True)

    axes[1].scatter(g_ts.time.datetime64, g_ts['diff_mag'], c='k', s=5, label='dmag')

    axes[0].legend()
    axes[1].legend()
    plt.savefig(directory+'/diff_phot.pdf', bbox_inches='tight')

    return



def diff_phot_ref(directory, refnum1, refnum2):
    ref = TimeSeries.read(os.path.join(directory, 'norm_ref_flux.ecsv'), time_column='time')

    ref['ref1_mag'] = -2.5 * np.log10(ref['ref_flux'][:,refnum1-1]) + 25
    ref['ref2_mag'] = -2.5 * np.log10(ref['ref_flux'][:,refnum2-1]) + 25

    g_ref = ref[ref['filter'] == 'g']

    fig, axes = plt.subplots(2, 1, figsize=(20,12))
    axes[0].scatter(g_ref.time.datetime64, g_ref['ref1_mag'], c='k', s=5, label='target ref')
    axes[0].scatter(g_ref.time.datetime64, g_ref['ref2_mag'], c='crimson', s=5, label='ref')

    g_ref['diff_mag'] = (g_ref['ref1_mag']) - (g_ref['ref2_mag'])
    g_ref.write(directory+'/diff_lc_ref.ecsv', overwrite=True)

    axes[1].scatter(g_ref.time.datetime64, g_ref['diff_mag'], c='k', s=5, label='dmag')
    
    axes[0].legend()
    axes[1].legend()
    plt.savefig(directory+'/diff_phot_ref.pdf', bbox_inches='tight')

    return

# %%
