# %%
import os
import numpy as np

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.stats import SigmaClip, sigma_clipped_stats, gaussian_fwhm_to_sigma
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.timeseries import TimeSeries
from astropy.table import QTable, vstack
import astropy.units as u

from skimage.draw import disk
from sklearn.preprocessing import normalize as sknorm
from photutils.background import Background2D, SExtractorBackground
from photutils.segmentation import detect_sources, deblend_sources, SourceCatalog
from photutils.aperture import CircularAperture, CircularAnnulus, ApertureStats
from scipy.signal import savgol_filter

import matplotlib.pylab as plt
from matplotlib.patches import Rectangle, Circle
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import Normalize
import matplotlib.dates as mdates
from astropy.visualization import ZScaleInterval, MinMaxInterval, PowerStretch, LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize



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
        assert type(visu_dir) == str
        fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(20,8))
        
        axes[0].set_title('Original Image')
        im1 = axes[0].imshow(img, cmap='RdYlBu_r', origin='lower', vmin=0, vmax=500)
        fig.colorbar(im1, ax=axes[0], location='bottom', pad=0.1)

        axes[1].set_title('Skymap Estimation')
        im2 = axes[1].imshow(bkg.background, cmap='RdYlBu_r', origin='lower')
        fig.colorbar(im2, ax=axes[1], location='bottom', pad=0.1)

        axes[2].set_title('Background Subtracted Image')
        im3 = axes[2].imshow(img - bkg.background, cmap='RdYlBu_r', origin='lower', vmin=0, vmax=500)
        fig.colorbar(im3, ax=axes[2], location='bottom', pad=0.1)

        if high_res == True:
            plt.savefig(os.path.join(visu_dir, 'skymap_visualization.pdf'), dpi=1200)
        else: 
            plt.savefig(os.path.join(visu_dir, 'skymap_visualization.pdf'))
        plt.show()

    return bkg.background



def detect_stars(img, wcs, radius=None, r_in=None, r_out=None, xbounds=None, ybounds=None, target_coord=None, 
                    visu=False, visu_dir=None, high_res=False):
    """
    Detect sources within a bounding box, while masking target RA-Dec.
    Return a list of centroid RA-Decs of detected sources.

    img: 2D numpy data array
    xbounds: bounding box [xmin, xmax] pixel positions
    ybounds: bounding box [ymin, ymax] pixel positions
    target_coord: target RA-Dec position to create circular aperture mask of r=15
    visu==True: visualize source detection and skymap estimation
    """
    if not radius:
        radius = 15

    # subtract 2D background
    skymap = estimate_background(img, mode='2D', visu=False)
    img -= skymap

    if xbounds and ybounds:
        # create image cutout according to pixel bounding box
        # note that Cutout2D takes position=(x,y) but size=(ny,nx)
        # https://docs.astropy.org/en/stable/nddata/utils.html#saving-a-2d-cutout-to-a-fits-file-with-an-updated-wcs
        cutout_center = ((xbounds[0]+xbounds[1])/2, (ybounds[0]+ybounds[1])/2)
        cutout_size = u.Quantity((ybounds[1]-ybounds[0], xbounds[1]-xbounds[0]), u.pixel)
        cutout = Cutout2D(img, position=cutout_center, size=cutout_size, wcs=wcs)
        img_crop = cutout.data
        # update wcs
        wcs = cutout.wcs
    else:
        img_crop = img

    if target_coord:
        # create circular aperture mask of r=30 pix around target RA-Dec
        # create mask=False for all pixel, then mask=True for target aperture
        # https://scikit-image.org/docs/stable/api/skimage.draw.html?highlight=disk#skimage.draw.disk
        pixel_coord = wcs.world_to_pixel(SkyCoord(target_coord))
        mask = np.zeros(img_crop.shape, dtype=bool)
        rr, cc = disk(center=(pixel_coord[1], pixel_coord[0]), radius=30)
        mask[rr, cc] = True
    else:
        # create null mask (mask=False for every pixel)
        mask = np.zeros(img_crop.shape, dtype=bool)

    # define detection threshold and convolve data
    _, _, std = sigma_clipped_stats(img_crop, sigma=3.0)
    threshold = 3 * std

    sigma = 3.0 * gaussian_fwhm_to_sigma  # FWHM = 3.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    convolved_img = convolve(img_crop, kernel)

    # npixels: how many connected pixels, each above threshold, should an area have to qualify as a source
    npixels = 10
    segment_map = detect_sources(convolved_img, threshold, npixels=npixels, mask=mask)
    # remove sources that partially overlap with 10 pix from border
    segment_map.remove_border_labels(border_width=10)
    segm_deblend = deblend_sources(convolved_img, segment_map, npixels=npixels, 
                                    nlevels=32, contrast=0.001)
    cat = SourceCatalog(img_crop, segm_deblend, convolved_data=convolved_img)
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


    if visu == True:
        assert type(visu_dir) == str
        # create grid for subplots
        fig = plt.figure()
        ax1 = plt.subplot2grid(shape=(3, 4), loc=(0, 0), colspan=3, rowspan=3)
        ax2 = plt.subplot2grid(shape=(3, 4), loc=(0, 3), colspan=1, rowspan=1)
        ax3 = plt.subplot2grid(shape=(3, 4), loc=(1, 3), colspan=1, rowspan=1)
        ax4 = plt.subplot2grid(shape=(3, 4), loc=(2, 3), colspan=1, rowspan=1)
        axes = [ax1, ax2, ax3, ax4]
        for ax in axes[1:]:
            ax.tick_params(labelbottom=False, labelleft=False)

        plt.suptitle('Background Subtraction and Source Detection', fontsize=12)

        # drawing main source detection plot
        # ----------------------------------
        fig.set_figheight(6)
        fig.set_figwidth(6)
        ax1.set_title('Source Detection (aperture mask, border rejection)', fontsize=10)
        # image normalization: interval and stretch
        norm = ImageNormalize(img, interval=ZScaleInterval(), stretch=PowerStretch(1.5))
        # show im0 (background image) in grey scale
        im0 = ax1.imshow(img, cmap='Greys', origin='lower', norm=norm)
        # define bounding box for im1 (image in front) and show im1 using colored cmap
        extent = [xbounds[0]+2, xbounds[1]+2, ybounds[0]+2, ybounds[1]+2]
        im1 = ax1.imshow(img_crop, cmap='RdYlBu_r', origin='lower', extent=extent, norm=norm)
        # scatter xy-centroid positions of detected sources
        ax1.scatter(table['xcentroid']+xbounds[0], table['ycentroid']+ybounds[0], 
                    s=10, marker='+', color='lime', linewidth=0.5)
        # draw bounding box
        ax1.add_patch(Rectangle((xbounds[0], ybounds[0]), xbounds[1]-xbounds[0], ybounds[1]-ybounds[0],
                                linestyle = 'dashed', edgecolor = 'black', 
                                fill=False, lw=0.3))
        # draw circular apertures for detected sources
        for n in range(len(table)):
            # define center pixel positions (cenx, ceny) and draw circles around each center
            cenx = table[n]['xcentroid']+xbounds[0]
            ceny = table[n]['ycentroid']+ybounds[0]
            circle = Circle((cenx, ceny), radius=radius, fill=False, color='white', lw=0.3)
            ax1.add_patch(circle)
            # annotate reference star number next to detection circles
            ax1.annotate(str(n+1), (cenx+25, ceny-30), color='black', fontsize=5, ha='center', va='center')
            if r_in and r_out:
                # draw circular annulus for detected sources
                ax1.add_patch(Circle((cenx, ceny), radius=r_in, fill=False, color='red', lw=0.3, alpha=0.8))
                ax1.add_patch(Circle((cenx, ceny), radius=r_out, fill=False, color='red', lw=0.3, alpha=0.8))
        if target_coord:
            # draw circular aperture for target RA-Dec
            ax1.add_patch(Circle((pixel_coord[0]+xbounds[0], pixel_coord[1]+ybounds[0]), 
                                    radius=radius, fill=False, color='fuchsia', lw=0.3))
            ax1.scatter(pixel_coord[0]+xbounds[0], pixel_coord[1]+ybounds[0], 
                                    s=10, marker=(5,2), color='yellow', lw=0.5)
            if r_in and r_out:
                # draw circular annulus for target RA-Dec
                ax1.add_patch(Circle((pixel_coord[0]+xbounds[0], pixel_coord[1]+ybounds[0]), 
                                    radius=r_in, fill=False, color='midnightblue', lw=0.3, alpha=0.8))
                ax1.add_patch(Circle((pixel_coord[0]+xbounds[0], pixel_coord[1]+ybounds[0]), 
                                    radius=r_out, fill=False, color='midnightblue', lw=0.3, alpha=0.8))
        # set tick parameters and color bar
        ax1.tick_params(axis='x', labelsize=8)
        ax1.tick_params(axis='y', labelsize=8)
        cbar1 = fig.colorbar(im1, ax=ax1, location='left', aspect=30)
        cbar1.ax.tick_params(labelsize=8)

        # plot original image with interval=[0,500]
        # -----------------------------------------
        fig.set_figheight(6)
        fig.set_figwidth(9)
        ax2.set_title('Original Image', fontsize=8)
        im2 = ax2.imshow(img + skymap, cmap='RdYlBu_r', origin='lower', vmin=0, vmax=500)
        # draw colorbar
        cbar2 = fig.colorbar(im2, ax=ax2, location='right')
        cbar2.ax.tick_params(labelsize=6)

        # plot background subtracted image with interval=[0,500]
        # ------------------------------------------------------
        ax3.set_title('Background Subtracted Image', fontsize=8)
        im3 = ax3.imshow(img, cmap='RdYlBu_r', origin='lower', vmin=0, vmax=500)
        # draw colorbar
        cbar3 = fig.colorbar(im3, ax=ax3, location='right')
        cbar3.ax.tick_params(labelsize=6)

        # plot skymap estimation from original image
        # ------------------------------------------
        ax4.set_title('Skymap Estimation', fontsize=8)
        im4 = ax4.imshow(skymap, cmap='RdYlBu_r', origin='lower')
        # draw colorbar
        cbar4 = fig.colorbar(im4, ax=ax4, location='right')
        cbar4.ax.tick_params(labelsize=6)

        # tight layout and output pdf
        plt.tight_layout()
        if high_res == True:
            plt.savefig(os.path.join(visu_dir, 'source_detection.pdf'), dpi=1200)
        else:
            plt.savefig(os.path.join(visu_dir, 'source_detection.pdf'))
        plt.show()


    # create sky_aperture object and sky_annulus object
    positions = []
    for n in range(len(table)):
        positions.append((table[n]['xcentroid'], table[n]['ycentroid']))
    aperture = CircularAperture(positions, r=radius)
    sky_aperture = aperture.to_sky(wcs)
    if r_in == None or r_out == None:
        return sky_aperture
    else:
        annulus_aperture = CircularAnnulus(positions, r_in=r_in, r_out=r_out)
        sky_annulus = annulus_aperture.to_sky(wcs)

    return sky_aperture, sky_annulus



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



def multithread_photometry(directory, cubename, sky_aperture, annulus=False, sky_annulus=None):
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



def LCs_visualizer(directory, visu_dir, mode='simple', ref_flux_table='ref_flux.ecsv', 
                            ylim=[0.85, 1.15], layout=[4, 6]):
    """
    Visualize lightcurves for reference stars. 
    """
    assert mode in ['simple', 'full', 'raw']

    table = TimeSeries.read(os.path.join(directory, ref_flux_table), time_column='time')
    # each row of table['ref_flux'] column contains a list of reference flux
    # table['ref_flux'] is thus a 2d numpy array; we can treat it as an image in all cases
    # from now on, table['ref_flux'] is refered to as 'reference flux map' (RFM)
    # x axis of RFM = number of reference stars
    # y axis of RFM = time-dependent flux
    RFM = table['ref_flux']

    # number of reference stars = width of RFM
    nrefs = RFM.shape[1]

    nrows, ncols = layout[0], layout[1]
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(28, 20))
    fig.suptitle('Reference Stars Lightcurves', fontsize=30, y=0.92)
    for nrow in range(nrows):
        for ncol in range(ncols):
            # n = reference star number
            n = nrow * ncols + ncol
            if mode in ['simple', 'full']:
                # normalize flux for each reference star = normalize RFM by column median
                norm_nflux = RFM[:, n] / np.median(RFM[:, n])
                # set the same ylim for all normalized flux plots
                axes[nrow][ncol].set_ylim(ylim)
                # scatter normalized flux write original median flux to legend
                axes[nrow][ncol].scatter(table.time.datetime64, norm_nflux, 
                                s=10, marker='x', color='k', linewidth=0.5, 
                                label='ref'+str(n+1)+', median≈'+str(int(round(np.median(RFM[:, n]),-2))))
                # smooth normalized flux by a Savitzky-Golay rolling mean of deg=1
                smooth_nflux = savgol_filter(norm_nflux, window_length=51, polyorder=1)
                # plot smoothed normalized flux curve in red
                axes[nrow][ncol].plot(table.time.datetime64, smooth_nflux, 
                                markersize=0, marker='.', color='tab:red', linewidth=2)
            elif mode == 'raw':
                # scatter raw flux for each reference star
                axes[nrow][ncol].scatter(table.time.datetime64, RFM[:, n], 
                                s=10, marker='x', color='k', linewidth=0.5, label='ref'+str(n+1))
                # smooth raw flux by a Savitzky-Golay rolling mean of deg=1
                smooth_nflux = savgol_filter(RFM[:, n], window_length=51, polyorder=1)
                # plot smoothed raw flux curve in red
                axes[nrow][ncol].plot(table.time.datetime64, smooth_nflux, 
                                markersize=0, marker='.', color='tab:red', linewidth=2)
            # format plot and axis labels
            axes[nrow][ncol].grid(True)
            axes[nrow][ncol].legend(loc='upper left', fontsize=14)
            axes[nrow][ncol].xaxis.set_major_locator(MaxNLocator(7))
            axes[nrow][ncol].xaxis.set_major_formatter(
                                    mdates.ConciseDateFormatter(axes[nrow][ncol].xaxis.get_major_locator()))
    plt.savefig(os.path.join(visu_dir, 'reference_stars_lc.pdf'))
    plt.close()


    # if we transpose RFM, its row = each reference star; column = time-dependent flux
    # the transposed RFM is refered to as 'spectrum'
    # the wavelength axis of a normal spectrum is now a time axis for the RFM spectrum
    spectrum = RFM.transpose()

    # sort spectrum by row median (median flux for each reference star) in descending order; get sort_index
    sort_index = np.median(spectrum, axis=1).argsort()[::-1]
    # normalize spectrum by row
    normalized_spectrum = sknorm(spectrum, axis=1)
    sorted_spec = normalized_spectrum[sort_index]
    sorted_spec = sorted_spec * (1 / np.median(sorted_spec))

    if mode == 'full':
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
        fig.suptitle('2D Reference Flux Map')
        # plot the flux-sorted spectrum as an image
        # set spectrum color range = ylim for normalized reference star lc
        im1 = ax.imshow(sorted_spec, cmap='viridis', aspect='auto', 
                        interpolation='None', vmin=ylim[0], vmax=ylim[1])
        fig.colorbar(im1, ax=ax, location='right', aspect=30, pad=0.03, 
                            label='normalized flux')

        # configure yticks and yticklabels to represent correct reference star number
        ax.set_yticks(np.arange(start=0, stop=nrefs))
        labels = [str(n+1) for n in sort_index]
        labels = ['ref'+m for m in labels]
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_ylabel('reference star number  - - - median flux + + +')
        ax.set_xlabel('Time', fontsize=10, labelpad=10)
        ax.tick_params(labelbottom=False)

        # tight layout and output png
        plt.tight_layout()
        plt.savefig(os.path.join(visu_dir, '2d_RFM.png'), dpi=1200)
        plt.close()


    fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(16,12))
    fig.suptitle('3D Reference Flux Map', y=0.92, fontsize=16)

    X = np.arange(start=1, stop=nrefs+1, step=1)
    Y = mdates.date2num(table.time.datetime64)

    for n in range(nrefs):
        # n = reference star number
        norm_nflux = RFM[:, n] / np.median(RFM[:, n])
        smooth_nflux = savgol_filter(norm_nflux, window_length=51, polyorder=1)
        if n == 0:
            ax.plot(np.linspace(start=n,stop=n,num=RFM.shape[0]), Y, smooth_nflux, 
                c='red', linewidth=0.8, alpha=0.8, label='Savitzky-Golay smoothing')
        else:
            ax.plot(np.linspace(start=n,stop=n,num=RFM.shape[0]), Y, smooth_nflux, 
                c='red', linewidth=0.8, alpha=0.8)
    ax.w_yaxis.set_major_locator(MaxNLocator(7))
    ax.w_yaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.w_yaxis.get_major_locator()))

    # configure yticks and yticklabels to represent correct reference star number
    ax.set_xticks(np.arange(start=0, stop=nrefs))
    labels = [str(n+1) for n in sort_index]
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel('DateTime', fontsize=10, labelpad=10)
    ax.set_xlabel('Reference Star Number  + + + median flux - - -', labelpad=10)
    ax.set_zlabel('Normalized Flux', labelpad=10)

    if mode == 'full':
        X, Y = np.meshgrid(X, Y)
        Z = sorted_spec.transpose()
        ax.set_zlim(ylim[0], ylim[1])
        norm = Normalize(vmin=ylim[0], vmax=ylim[1])
        surf = ax.scatter(X-1, Y, Z, c=Z, s=1, norm=norm, cmap='viridis', linewidths=0, alpha=1, label='raw flux')
        fig.colorbar(surf, ax=ax, location='left', shrink=0.5, aspect=30, label='normalized flux')
    
    ax.legend()
    plt.savefig(os.path.join(visu_dir, '3d_RFM.pdf'), dpi=1200)
    plt.show()

    return table, RFM



def differential_photometry(directory, target_flux='target_flux.ecsv', ref_flux='ref_flux.ecsv', reflist=None):
    """
    Differential photometry from selected reference stars. 
    Output diff_lc.
    """
    return



def visualize_lightcurve(directory, diff_lc='diff_lc.ecsv'):
    """
    Visualize lightcurve from differential photometry. 
    """
    return













