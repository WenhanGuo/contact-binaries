import os

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