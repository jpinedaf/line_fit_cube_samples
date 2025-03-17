from astropy.io import fits
from spectral_cube import SpectralCube

# files to be regridded
fits_files = ['cube1.fits', 'cube2.fits']
# suffix for the new files
suffix = '_regridded'
# template file
template_file = 'template.fits'

# load the template cube
# check if the template file has the necessary header keys
# there might be more keys that need to be checked
template_hd = fits.getheader(template_file)
list_key = ['NAXIS1', 'NAXIS2', 'CRVAL1', 'CRVAL2', 'CDELT1', 'CDELT2',
            'CRPIX1', 'CRPIX2', 'CTYPE1', 'CTYPE2', 'CUNIT1', 'CUNIT2',]
for key in list_key:
    if key not in template_hd:
        raise KeyError(f'Header key {key} not found in the template file.')
# load all cubes and regrid them
for f in fits_files:
    cube = SpectralCube.read(f)
    target_hd = cube.header
    for key in list_key:
        target_hd[key] = template_hd[key]
    regridded_cube = cube.reproject(target_hd)
    regridded_cube.write(
        f.replace('.fits', f'{suffix}.fits'), overwrite=True)
