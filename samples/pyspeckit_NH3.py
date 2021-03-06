import pyspeckit
import astropy.io.fits as fits
import numpy as np

from spectral_cube import SpectralCube
from radio_beam import Beam
import astropy.units as u
from skimage.morphology import remove_small_objects,closing,disk,opening

from pyspeckit.spectrum.models import ammonia


fit_dir = 'fit/'
data_dir = 'data/'

OneOneIntegrated = data_dir + 'NH3_11_TdV.fits'
OneOneMom1 = data_dir + 'NH3_11_mom1.fits'
OneOneMom2 = data_dir + 'NH3_11_mom2.fits'
OneOneFile = data_dir + 'NH3_11.fits'
OneOnePeak = data_dir + 'NH3_11_Tpeak.fits'
RMSFile_11 = data_dir + 'NH3_11_rms.fits'
TwoTwoFile = data_dir + 'NH3_22.fits'
RMSFile_22 = data_dir + 'NH3_22_rms.fits'
# Define frequencies used to determine the velocity
freq11 = 23.694506*u.GHz
freq22 = 23.722633335*u.GHz


def mom_map(do_plot=False):
    """
    AAA
    """
    cube11sc = SpectralCube.read(OneOneFile)
    cube11_v = cube11sc.with_spectral_unit(u.km/u.s, velocity_convention='radio')
    chan = np.arange(cube11sc.header['NAXIS3'])

    slab = cube11_v.spectral_slab(3.44*u.km/u.s, 5.0*u.km/u.s)
    w11 = slab.moment(order=0, axis=0)
    moment1 = slab.moment1(axis=0)
    moment2 = slab.moment2(axis=0)
    Tpeak = slab.max(axis=0)
    w11.write(OneOneIntegrated, format='fits', overwrite=True)
    moment1.write(OneOneMom1, format='fits', overwrite=True)
    moment2.write(OneOneMom2, format='fits', overwrite=True)
    Tpeak.write(OneOnePeak, format='fits', overwrite=True)

    a=[ 5,  72, 251, 373, 490, 666]
    b=[45, 220, 335, 450, 630, 683]
    index_rms = np.hstack(np.arange(start, stop+1) for start, stop in zip(a, b))

    mask_rms = np.zeros(cube11_v.shape, dtype=bool)
    mask_rms[index_rms] = True
    mask_rms = mask_rms & np.isfinite((cube11_v.unmasked_data[:,:,:]).value)
    cube11_rms = cube11_v.with_mask(mask_rms)
    rms11 = cube11_rms.std(axis=0)
    rms11.write(RMSFile_11, overwrite=True)
    if do_plot:
        import matplotlib.pyplot as plt
        # compare full_cube and the one only with noise (no signal)
        ii = 153; jj = 144
        plt.plot(chan, cube11sc[:, jj, ii].value, 'red', chan, cube11_rms[:, jj, ii], 'blue')
        plt.show()
        plt.imshow(Tpeak.value/rms11.value, origin='lower', vmin=0.3, vmax=10)
        plt.show()
        plt.close()
    cube22sc = SpectralCube.read(TwoTwoFile)
    spectral_axis_22 = cube22sc.with_spectral_unit(u.km/u.s, velocity_convention='radio').spectral_axis
    good_channels_22 = (spectral_axis_22 < 3.5*u.km/u.s) | (spectral_axis_22 > 5.0*u.km/u.s)
    masked_cube22 = cube22sc.with_mask(good_channels_22[:, np.newaxis, np.newaxis])
    rms22 = masked_cube22.std(axis=0)
    rms22.write(RMSFile_22, overwrite=True)


def hmm1_cubefit(vmin=3.4, vmax=5.0, tk_ave=10., do_plot=False, snr_min=5.0,
        multicore=1, do_thin=False):
    """
    Fit NH3(1,1) and (2,2) cubes for H-MM1.
    It fits all pixels with SNR larger than requested. 
    Initial guess is based on moment maps and neighboring pixels. 
    The fitting can be done in parallel mode using several cores, 
    however, this is dangerous for large regions, where using a 
    good initial guess is important. 
    It stores the result in a FITS cube. 

    TODO:
    -convert FITS cube into several FITS files
    -Improve initial guess
    
    Parameters
    ----------
    vmin : numpy.float
        Minimum centroid velocity to plot, in km/s.
    vmax : numpy.float
        Maximum centroid velocity to plot, in km/s.
    tk_ave : numpy.float
        Mean kinetic temperature of the region, in K.
    do_plot : bool
        If True, then a map of the region to map is shown.
    snr_min : numpy.float
        Minimum signal to noise ratio of the spectrum to be fitted.
    multicore : int
        Numbers of cores to use for parallel processing. 
    """

    cube11sc = SpectralCube.read(OneOneFile)
    cube22sc = SpectralCube.read(TwoTwoFile)
    cube11_v = cube11sc.with_spectral_unit(u.km/u.s,velocity_convention='radio',
                                           rest_value=freq11)
    cube22_v = cube22sc.with_spectral_unit(u.km/u.s,velocity_convention='radio',
                                           rest_value=freq22)
    from pyspeckit.spectrum.units import SpectroscopicAxis
    spec11 = SpectroscopicAxis(cube11_v.spectral_axis, refX=freq11,
                                  velocity_convention='radio')
    spec22 = SpectroscopicAxis(cube22_v.spectral_axis, refX=freq22,
                                  velocity_convention='radio')

    errmap11 = fits.getdata(RMSFile_11)
    errmap22 = fits.getdata(RMSFile_22)
    errmap_K = errmap11 #[errmap11, errmap22]
    Tpeak11 = fits.getdata(OneOnePeak)

    moment1 = fits.getdata(OneOneMom1)
    moment2 = (fits.getdata(OneOneMom2))**0.5

    snr = cube11sc.filled_data[:].value/errmap11
    peaksnr = Tpeak11 / errmap11

    planemask = (peaksnr > snr_min) # *(errmap11 < 0.15)
    planemask = remove_small_objects(planemask, min_size=40)
    planemask = opening(planemask, disk(1))
    #planemask = (peaksnr>20) * (errmap11 < 0.2)

    mask = (snr>3)*planemask

    maskcube = cube11sc.with_mask(mask.astype(bool))
    maskcube = maskcube.with_spectral_unit(u.km/u.s, velocity_convention='radio')
    slab = maskcube.spectral_slab(vmax*u.km/u.s, vmin*u.km/u.s)
    w11 = slab.moment(order=0, axis=0).value
    peakloc = np.nanargmax(w11)
    ymax, xmax = np.unravel_index(peakloc, w11.shape)

    moment2[np.isnan(moment2)] = 0.2
    moment2[moment2 < 0.2] = 0.2

    ## Load FITS files
    cube11 = pyspeckit.Cube(OneOneFile, maskmap=planemask)
    cube22 = pyspeckit.Cube(TwoTwoFile, maskmap=planemask)
    # Stack files
    cubes = pyspeckit.CubeStack([cube11, cube22], maskmap=planemask)
    cubes.unit = "K"
    # Define initial guess
    guesses = np.zeros((6,) + cubes.cube.shape[1:])
    moment1[moment1 < vmin] = vmin+0.2
    moment1[moment1 > vmax] = vmax-0.2
    guesses[0, :, :] = tk_ave                # Kinetic temperature 
    guesses[1, :, :] = 7                     # Excitation  Temp
    guesses[2, :, :] = 14.5                  # log(column)
    guesses[3, :, :] = moment2  # Line width / 5 (the NH3 moment overestimates linewidth)               
    guesses[4, :, :] = moment1  # Line centroid              
    guesses[5, :, :] = 0.5                   # F(ortho) - ortho NH3 fraction (fixed)
    if do_plot:
        import matplotlib.pyplot as plt
        plt.imshow(w11*planemask, origin='lower')
        plt.show()
    print('start fit')
    cubes.specfit.Registry.add_fitter('cold_ammonia',
            ammonia.cold_ammonia_model(), 6)
    if do_thin:
        file_out = "{0}H-MM1_cold_parameter_maps_snr{1}_thin_v1.fits".format(fit_dir, snr_min)
    else:
        file_out = "{0}H-MM1_cold_parameter_maps_snr{1}_thick_v1.fits".format(fit_dir, snr_min)
    cubes.fiteach(fittype='cold_ammonia', guesses=guesses,
                  integral=False, verbose_level=3,
                  fixed=[do_thin, False, False, False, False, True],
                  signal_cut=2,
                  limitedmax=[True, False, False, False, True, True],
                  maxpars=[20, 15, 20, 0.4, vmax, 1],
                  limitedmin=[True, True, True, True, True, True],
                  minpars=[5, 2.8, 12.0, 0.05, vmin, 0],
                  start_from_point=(xmax, ymax),
                  use_neighbor_as_guess=True,
                  position_order=1/peaksnr,
                  errmap=errmap_K, multicore=multicore)
    # Store fits into FITS cube
    fitcubefile = fits.PrimaryHDU(
                    data=np.concatenate([cubes.parcube, cubes.errcube]),
                    header=cubes.header)
    fitcubefile.header.set('PLANE1', 'TKIN')
    fitcubefile.header.set('PLANE2', 'TEX')
    fitcubefile.header.set('PLANE3', 'COLUMN')
    fitcubefile.header.set('PLANE4', 'SIGMA')
    fitcubefile.header.set('PLANE5', 'VELOCITY')
    fitcubefile.header.set('PLANE6', 'FORTHO')
    fitcubefile.header.set('PLANE7', 'eTKIN')
    fitcubefile.header.set('PLANE8', 'eTEX')
    fitcubefile.header.set('PLANE9', 'eCOLUMN')
    fitcubefile.header.set('PLANE10', 'eSIGMA')
    fitcubefile.header.set('PLANE11', 'eVELOCITY')
    fitcubefile.header.set('PLANE12', 'eFORTHO')
    fitcubefile.header.set('CDELT3', 1)
    fitcubefile.header.set('CTYPE3', 'FITPAR')
    fitcubefile.header.set('CRVAL3', 0)
    fitcubefile.header.set('CRPIX3', 1)
    fitcubefile.writeto(file_out, overwrite=True)
