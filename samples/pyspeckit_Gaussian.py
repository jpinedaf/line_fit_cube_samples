import numpy as np
import pyspeckit
from astropy.io import fits
import astropy.units as u

data_dir = 'data/'
fit_dir = 'fit/'

file_in = data_dir + 'HC3N_9-8.fits'
rms_file = data_dir + 'HC3N_9-8_rms.fits'
file_par = fit_dir + 'HC3N_9-8_fitted.fits'

freq_line = 81.87992490 * u.GHz
snr_min = 8.0
xmax = 28; ymax = 32
# range of parameters
vmin = 1.0; vmax = 4.1
tpeak_min = 0.0; tpeak_max = 80e-3
dv_min = 0.02; dv_max = 0.9

# Load cube and setup the spectral axis
cube = pyspeckit.Cube(file_in)
cube.xarr.refX = freq_line
cube.xarr.velocity_convention = 'radio'
cube.xarr.convert_to_unit('km/s')

err_map = fits.getdata(rms_file)
# Initial guess value
guesses = np.array([1.5, 5.5, 0.3])

cube.fiteach(guesses=guesses,
    start_from_point=(xmax, ymax),
    fittype='gaussian', blank_value=np.nan,
    use_neighbor_as_guess=True,
    limitedmax=[True, True, True],
    limitedmin=[True, True, False],
    maxpars=[tpeak_max, vmax, dv_max],
    minpars=[tpeak_min, vmin, dv_min],
    multicore=4, errmap=err_map, signal_cut=snr_min)

cube.write_fit(file_par, overwrite=True)
