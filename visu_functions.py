import os

import numpy as np
from scipy.signal import savgol_filter
from sklearn.preprocessing import normalize as sknorm

import mpld3
import matplotlib.pylab as plt
from matplotlib.patches import Rectangle, Circle
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import Normalize
import matplotlib.dates as mdates
from astropy.visualization import ZScaleInterval, MinMaxInterval, PowerStretch, LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize

from IPython.display import set_matplotlib_formats
set_matplotlib_formats('retina')


def visualize_background(img, skymap, visu_dir, high_res=False):
    """
    Visualize background subtraction (original image, subtracted image, skymap estimation)

    img: original image (2D numpy data array). if using hdu: hdu.data
    skymap: 2D background estimation, must be the same shape as img
    visu_dir: directory to save visualization results
    high_res: if want to save high resolution results; takes longer time
    """
    assert type(visu_dir) == str

    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(20,8))
    
    axes[0].set_title('Original Image')
    im1 = axes[0].imshow(img, cmap='RdYlBu_r', origin='lower', vmin=0, vmax=500)
    fig.colorbar(im1, ax=axes[0], location='bottom', pad=0.1)

    axes[1].set_title('Skymap Estimation')
    im2 = axes[1].imshow(skymap, cmap='RdYlBu_r', origin='lower')
    fig.colorbar(im2, ax=axes[1], location='bottom', pad=0.1)

    axes[2].set_title('Background Subtracted Image')
    im3 = axes[2].imshow(img - skymap, cmap='RdYlBu_r', origin='lower', vmin=0, vmax=500)
    fig.colorbar(im3, ax=axes[2], location='bottom', pad=0.1)

    if high_res == True:
        plt.savefig(os.path.join(visu_dir, 'skymap_visualization.pdf'), dpi=1200)
    else: 
        plt.savefig(os.path.join(visu_dir, 'skymap_visualization.pdf'))
    plt.show()

    return



def visualize_detection(img, masked_img, skymap, table, 
                        radius, r_in, r_out, xbounds, ybounds, pixel_coord, 
                        visu_dir, high_res=False):
    """
    Visualize source detection and background subtraction for one input image.

    img: original image (2D numpy data array). if using hdu: hdu.data
    masked_img: masked image indicating bounding box
    skymap: 2D background estimation, must be the same shape as img
    table: photutils source detection table, must have columns of 'xcentroid' and 'ycentroid'
    radius: visualize apertures on source centroids; aperture photometry radius in pixels
    r_in: annulus photometry inner radius
    r_out: annulus photometry outer radius
    xbounds: horizontal bounding box [xmin, xmax] pixel positions
    ybounds: vertical bounding box [ymin, ymax] pixel positions
    pixel_coord: if not None: pixel coordinates to draw target star independent from source detection
    visu_dir: directory to save visualization results
    high_res: if want to save high resolution results; takes longer time
    """
    assert type(visu_dir) == str
    if xbounds == None and ybounds == None:
        xbounds = [0, img.shape[1]]
        ybounds = [0, img.shape[0]]

    # create grid for subplots
    # ------------------------
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
    # extent = [xbounds[0]-0.5, xbounds[1]-0.5, ybounds[0]-1, ybounds[1]-1]
    extent=[-0.5, img.shape[1]-0.5, -0.5, img.shape[0]-0.5]
    im1 = ax1.imshow(masked_img, cmap='RdYlBu_r', origin='lower', 
                    extent=extent, norm=norm, interpolation='none')
    # scatter xy-centroid positions of detected sources
    ax1.scatter(table['xcentroid'], table['ycentroid'], 
                s=10, marker='+', color='lime', linewidth=0.5)
    # draw bounding box
    ax1.add_patch(Rectangle((xbounds[0], ybounds[0]), xbounds[1]-xbounds[0], ybounds[1]-ybounds[0],
                            linestyle = 'dashed', edgecolor = 'black', 
                            fill=False, lw=0.3))
    # draw circular apertures for detected sources
    for n in range(len(table)):
        # define center pixel positions (cenx, ceny) and draw circles around each center
        cenx = table[n]['xcentroid']
        ceny = table[n]['ycentroid']
        circle = Circle((cenx, ceny), radius=radius, fill=False, color='white', lw=0.3)
        ax1.add_patch(circle)
        # annotate reference star number next to detection circles
        ax1.annotate(str(n+1), (cenx+25, ceny-30), color='black', fontsize=5, ha='center', va='center')
        if r_in and r_out:
            # draw circular annulus for detected sources
            ax1.add_patch(Circle((cenx, ceny), radius=r_in, fill=False, color='red', lw=0.3, alpha=0.8))
            ax1.add_patch(Circle((cenx, ceny), radius=r_out, fill=False, color='red', lw=0.3, alpha=0.8))
    if pixel_coord is not None:
        # draw circular aperture for target RA-Dec
        ax1.add_patch(Circle((pixel_coord[0], pixel_coord[1]), 
                        radius=radius, fill=False, color='fuchsia', lw=0.3))
        ax1.scatter(pixel_coord[0], pixel_coord[1], s=10, marker=(5,2), color='yellow', lw=0.5)
        if r_in and r_out:
            # draw circular annulus for target RA-Dec
            ax1.add_patch(Circle((pixel_coord[0], pixel_coord[1]), 
                                radius=r_in, fill=False, color='midnightblue', lw=0.3, alpha=0.8))
            ax1.add_patch(Circle((pixel_coord[0], pixel_coord[1]), 
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
        plt.savefig(os.path.join(visu_dir, 'source_detection.png'), dpi=1200)
        
    else:
        plt.savefig(os.path.join(visu_dir, 'source_detection.pdf'))
    plt.show()

    return



def visualize_1D_RFM(table, mode, visu_dir, ylim=[0.85, 1.15], layout=[4, 6]):
    """
    Visualize reference stars' lightcurves in traditional time-flux subplots.
    Please see phot_functions.py for an explanation of RFM
    """
    assert mode in ['norm', 'raw']

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
            if n < nrefs:
                if mode == 'norm':
                    # normalize flux for each reference star = normalize RFM by column median
                    norm_nflux = RFM[:, n] / np.median(RFM[:, n])
                    # set the same ylim for all normalized flux plots
                    axes[nrow][ncol].set_ylim(ylim)
                    # scatter normalized flux write original median flux to legend
                    axes[nrow][ncol].scatter(table.time.datetime64, norm_nflux, 
                                s=10, marker='x', color='k', linewidth=0.5, 
                                label='ref'+str(n+1)+', median???'+str(int(round(np.median(RFM[:, n]),-2))))
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
    # save in pdf and in html (interactive)
    plt.savefig(os.path.join(visu_dir, 'reference_stars_lc.pdf'))
    mpld3.save_html(fig, os.path.join(visu_dir, 'reference_stars_lc.html'))
    plt.close()

    return



def visualize_2D_RFM(table, visu_dir, ylim=[0.85, 1.15]):
    """
    Visualize reference stars' lightcurves in a 2D image: x-axis=time, y-axis=nref, color=flux
    Please see phot_functions.py for an explanation of RFM
    """
    RFM = table['ref_flux']
    nrefs = RFM.shape[1]

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

    return



def visualize_3D_RFM(table, mode, visu_dir, ylim):
    """
    Visualize reference stars' lightcurves in 3D. 
    Please see phot_functions.py for an explanation of RFM
    """
    assert mode in ['simple', 'full scatter']
    
    # initialization; please see visualize_2D_RFM for explanation
    RFM = table['ref_flux']
    nrefs = RFM.shape[1]
    spectrum = RFM.transpose()
    sort_index = np.median(spectrum, axis=1).argsort()[::-1]
    normalized_spectrum = sknorm(spectrum, axis=1)
    sorted_spec = normalized_spectrum[sort_index]
    sorted_spec = sorted_spec * (1 / np.median(sorted_spec))


    fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(16,12))
    fig.suptitle('3D Reference Flux Map', y=0.92, fontsize=16)

    # X-axis: reference star number. the np.arange here is a locator; does not correlate to ref number
    # Y-axis: mpl num conversion of np.datetime64 timestamps
    X = np.arange(start=1, stop=nrefs+1, step=1)
    Y = mdates.date2num(table.time.datetime64)

    # loop through reference stars according to sort order, brightest first
    for n in sort_index:
        # n = reference star number
        # normalize each star's flux by its median flux
        norm_nflux = RFM[:, n] / np.median(RFM[:, n])
        # smooth normalized flux by a Savitzky-Golay rolling mean of deg=1
        smooth_nflux = savgol_filter(norm_nflux, window_length=51, polyorder=1)
        # reinitialize an X-axis locator for each star to plot on specific x
        Xlocator = sort_index.tolist().index(n)
        if n == sort_index[0]:
            # plot smoothed curve for brightest star with label
            ax.plot(np.linspace(start=Xlocator,stop=Xlocator,num=RFM.shape[0]), Y, smooth_nflux, 
                c='red', linewidth=0.8, alpha=0.8, label='Savitzky-Golay smoothing')
        else:
            # plot smoothed curve for nth brightest star one by one
            ax.plot(np.linspace(start=Xlocator,stop=Xlocator,num=RFM.shape[0]), Y, smooth_nflux, 
                c='red', linewidth=0.8, alpha=0.8)
    # set Y-axis format to concise DateTime; limit max=7 time labels for Y-axis
    ax.w_yaxis.set_major_locator(MaxNLocator(7))
    ax.w_yaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.w_yaxis.get_major_locator()))

    # configure xticks and xticklabels to represent correct reference star number
    ax.set_xticks(np.arange(start=0, stop=nrefs))
    labels = [str(n+1) for n in sort_index]
    ax.set_xticklabels(labels, fontsize=8)

    # set XYZ axes labels
    ax.set_ylabel('DateTime', fontsize=10, labelpad=10)
    ax.set_xlabel('Reference Star Number  + + + median flux - - -', labelpad=10)
    ax.set_zlabel('Normalized Flux', labelpad=10)

    # scatter raw flux data points if mode == full scatter
    if mode == 'full scatter':
        # meshgrid for X, Y to get img shape; correlation: one (x,y) pair --> one z
        X, Y = np.meshgrid(X, Y)
        Z = sorted_spec.transpose()

        # set zlim (flux/display color) to the same as 1D ylim and 2D color vmin-vmax
        ax.set_zlim(ylim[0], ylim[1])
        # normalize colorbar
        norm = Normalize(vmin=ylim[0], vmax=ylim[1])
        # scatter raw flux, color --> z, disable depth of field (alpha=1) for performance speed
        surf = ax.scatter(X-1, Y, Z, c=Z, s=1, norm=norm, cmap='viridis', 
                            linewidths=0, alpha=1, label='raw flux')
        fig.colorbar(surf, ax=ax, location='left', shrink=0.5, aspect=30, label='normalized flux')
    
    # legend and save pdf
    ax.legend()
    plt.savefig(os.path.join(visu_dir, '3d_RFM.pdf'), dpi=1200)
    plt.show()
    
    return



