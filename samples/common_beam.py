import radio_beam
import astropy.units as u
from astropy.io import fits
from spectral_cube import SpectralCube

# takes 2 files, and find the common beam
# and convolve both to that beam
# then write out the new files

# files to be convolved
fits_files = ['cube1.fits', 'cube2.fits']
# suffix for the new files
suffix = '_convolved'

# load the parameters of the beams
bmaj_list, bmin_list, bpa_list = [], [], []
for file_i in fits_files:
    hd = fits.getheader(file_i)
    bmaj_list.append(hd['BMAJ'])
    bmin_list.append(hd['BMIN'])
    bpa_list.append(hd['BPA'])

# create the Beam objects
my_beams = radio_beam.Beams(
    major=bmaj_list * u.deg, minor=bmin_list * u.deg, pa=bpa_list * u.deg)
common_beam = my_beams.common_beam()

# write out the smoothed cubes
for file_in in fits_files:
    cube = SpectralCube.read(file_in)
    new_cube = cube.convolve_to(common_beam)
    new_cube.write(file_in.replace(
        '.fits', f'{suffix}.fits'), overwrite=True)
