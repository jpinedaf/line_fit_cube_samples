from spectral_cube import SpectralCube
import numpy as np
import astropy.units as u

# file without the primary beam corrected
# then primary beam response, and primary beam corrected one
file_in = 'cube_no_PBcorr.fits'
file_in_PB = 'cube_PB_response.fits'
file_in_K = 'cube_ready.fits'

cube = SpectralCube.read(file_in)
PB = fits.getdata(file_in_PB)
kcube = cube.to(u.K) / np.squeeze(PB)
kcube.write(file_in_K, overwrite=True)
