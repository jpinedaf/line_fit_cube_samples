import pyspeckit
from pyspeckit.spectrum.models import nh2d
import numpy as np
import astropy.units as u
from astropy.io import fits
from skimage.morphology import remove_small_objects,closing,disk,opening

data_dir = 'data/'
fit_dir = 'fit/'
# primary beam corrected file and in Kelvin units
file_in_K = data_dir + 'oNH2D.fits'

snr_min = 5
file_thick = fit_dir + 'oNH2D_fit_thick_par_snr{0}.fits'.format(snr_min)
file_thin = fit_dir + 'oNH2D_fit_thin_par_snr{0}.fits'.format(snr_min)
rms_file = data_dir + 'oNH2D_rms.fits'
SNR_file = data_dir + 'oNH2D_SNR.fits'
Tpeak_file = data_dir + 'oNH2D_Tpeak.fits'
mask_file = data_dir + 'oNH2D_mask.fits'

# rest-freq used in Harju et al. (2020) for ortho-NH2D
freq_line = 85.92627*u.GHz


cube = pyspeckit.Cube(file_in_K)
cube.xarr.refX = freq_line
cube.xarr.velocity_convention = 'radio'
cube.xarr.convert_to_unit('km/s')

#
# source dependent parameters
#
# pixel to start fit and to make sample plot of the line fit 
xmax = 269; ymax = 240
# velocity range with signal (of the main hf-component) to 
# calculate moment maps
vmin = 3.4; vmax = 5.0
# default velocity used as guess in fit
vmean = 4.3
# velocity range used to plot velocity map
vmin_plot = 4.0; vmax_plot = 4.4
# velocity range used to calculate rms
vrms_0 = 0.; vrms_1 = 1.5

multicore = 40
plot_dv = False

# Prepare extra files
Prepare_Files = False
# Optically thin
Optically_Thin = False
Show_Optically_Thin = False
# Optically thick
Optically_Thick = False
Show_Optically_Thick = False

if Prepare_Files:
    rms_map = cube.slice(vrms_0, vrms_1, unit='km/s').cube.std(axis=0)
    Tpeak =  cube.slice(vmin, vmax, unit='km/s').cube.max(axis=0)
    peaksnr =  Tpeak / rms_map
    # create mask using a signal-to-noise cut
    planemask = (peaksnr > snr_min)
    # remove some small and isolated pixels, this is data dependent
    planemask = remove_small_objects(planemask, min_size=80)
    planemask = opening(planemask, disk(1))

    hd_cube=cube.header.copy()
    key_remove=['NAXIS3','CRPIX3','CDELT3','CUNIT3','CTYPE3','CRVAL3','SPECSYS']
    for key_i in key_remove:
        hd_cube.remove(key_i)
    hd_cube['WCSAXES'] = 2
    fits.writeto(rms_file, rms_map, hd_cube, overwrite=True)
    fits.writeto(Tpeak_file, Tpeak, hd_cube, overwrite=True)
    fits.writeto(SNR_file, peaksnr, hd_cube, overwrite=True)
    fits.writeto(mask_file, planemask.astype(int), hd_cube, overwrite=True)
else:
    planemask = fits.getdata(mask_file)
    Tpeak = fits.getdata(Tpeak_file)
    rms_map = fits.getdata(rms_file)
    peaksnr =  Tpeak/rms_map

import matplotlib.pyplot as plt
plt.ion()


if Optically_Thin:
    cube.Registry.add_fitter('nh2d_vtau', pyspeckit.spectrum.models.nh2d.nh2d_vtau_fitter, 4)

    print('start optically thin fit')
    cube.fiteach(fittype='nh2d_vtau',  guesses=[15.0, 0.1, vmean, 0.12],
                 verbose_level=1, signal_cut=snr_min,
                 limitedmin=[True, True, True, True],
                 limitedmax=[True, False, True, True],
                 minpars=[2.8, 0, vmin, 0.05],
                 maxpars=[250.0, 0, vmax, 1.0],
                 fixed=[False, True, False, False], 
                 use_neighbor_as_guess=True, 
                 start_from_point=(xmax, ymax),
                 errmap=rms_map, 
                 maskmap=planemask,
                 multicore=multicore)
    cube.write_fit(file_thin, overwrite=True)


if Optically_Thick:
    cube.Registry.add_fitter('nh2d_vtau', pyspeckit.spectrum.models.nh2d.nh2d_vtau_fitter, 4)

    print('start optically thick fit')
    cube.fiteach(fittype='nh2d_vtau',  guesses=[7.0, 2.0, vmean, 0.12],
                 verbose_level=2, signal_cut=snr_min,
                 limitedmin=[True, True, True, True],
                 limitedmax=[True, True, True, True],
                 minpars=[2.8, 0, vmin, 0.05],
                 maxpars=[20., 50, vmax, 1.0],
                 fixed=[False, False, False, False], 
                 use_neighbor_as_guess=True, 
                 start_from_point=(xmax, ymax),
                 errmap=rms_map, 
                 maskmap=planemask,
                 multicore=multicore)
    cube.write_fit(file_thick, overwrite=True)


if Show_Optically_Thin:
    #
    cube.Registry.add_fitter('nh2d_vtau', pyspeckit.spectrum.models.nh2d.nh2d_vtau_fitter, 4)
    cube.load_model_fit(file_thin, npars=4, npeaks=1, _temp_fit_loc=(xmax, ymax))
    cube.mapplot()
    cube.plot_spectrum(xmax, ymax, plot_fit=True)
    if plot_dv:
        cube.mapplot.plane = cube.parcube[3, :, :]
        cube.mapplot(estimator=None, vmin=0.05, vmax=0.25)
    else:
        cube.mapplot.plane = cube.parcube[2, :, :]
        cube.mapplot(estimator=None, vmin=vmin_plot, vmax=vmax_plot,
            cmap='RdYlBu_r')
    plt.draw()
    plt.show()
    plt.savefig('NH2D_pyspeckit_poor_thin_fit.png')


if Show_Optically_Thick:
    #
    cube.Registry.add_fitter('nh2d_vtau', pyspeckit.spectrum.models.nh2d.nh2d_vtau_fitter, 4)
    cube.load_model_fit( file_thick, npars=4, npeaks=1, _temp_fit_loc=(xmax, ymax))
    cube.mapplot()
    cube.plot_spectrum(xmax, ymax, plot_fit=True)
    if plot_dv:
        cube.mapplot.plane = cube.parcube[3, :, :]
        cube.mapplot(estimator=None, vmin=0.05, vmax=0.15)
    else:
        cube.mapplot.plane = cube.parcube[2, :, :]
        cube.mapplot(estimator=None, vmin=vmin_plot, vmax=vmax_plot,
            cmap='RdYlBu_r')
    plt.draw()
    plt.show()
    plt.savefig('NH2D_pyspeckit_poor_thick_fit.png')
