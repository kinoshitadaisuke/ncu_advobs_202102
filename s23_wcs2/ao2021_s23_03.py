#!/usr/pkg/bin/python3.9

#
# Time-stamp: <2021/06/10 22:18:23 (CST) daisuke>
#

# importing argparse module
import argparse

# importing sys module
import sys

# importing pathlib module
import pathlib

# importing numpy module
import numpy

# importing astropy module
import astropy.coordinates
import astropy.modeling
import astropy.units
import astropy.wcs

# choice of PSF model
list_psf_model = ['2dg', '2dm']

# constructing parser object
desc = "measuring position on images"
parser = argparse.ArgumentParser (description=desc)

# adding arguments
parser.add_argument ('-m', '--psf-model', choices=list_psf_model, \
                     default='2dg', \
                     help='PSF model [2dg=Gaussian, 2dm=Moffat] (default: 2dg)')
parser.add_argument ('-p', '--position', default='', \
                     help='RA and Dec in hh:mm:ss.ss,+/-dd:mm:ss.s format')
parser.add_argument ('-w', '--half-width', type=int, default=5, \
                     help='half-width of centroid box (default: 30)')
parser.add_argument ('FITS_files', nargs='+', help='FITS files')

# command-line argument analysis
args = parser.parse_args ()

# input parameters
psf_model  = args.psf_model
position   = args.position
half_width = args.half_width
files_fits = args.FITS_files

# check of input position
if not ( (',' in position) and (':' in position) ):
    print ("RA and Dec must be given in hh:mm:ss.ss,+/-dd:mm:ss.s format.", \
           file=sys.stderr)
    sys.exit ()

# function to read a FITS file
def read_fits (file_fits):
    # counter
    i = 0
    # reading FITS file
    with astropy.io.fits.open (file_fits) as hdu:
        # reading header and image
        header = hdu[i].header
        image  = hdu[i].data
        # if no image in PrimaryHDU, then read next HDU
        while (header['NAXIS'] == 0):
            i += 1
            header = hdu[i].header
            image  = hdu[i].data
    # returning header and image
    return (header, image)

# units
u_ha  = astropy.units.hourangle
u_deg = astropy.units.deg
    
# RA and Dec
ra_str, dec_str = position.split (',')

# making SkyCoord object
coord = astropy.coordinates.SkyCoord (ra_str, dec_str, \
                                      unit=(u_ha, u_deg), frame='icrs')

# RA and Dec in deg
ra_deg  = coord.ra.deg
dec_deg = coord.dec.deg

# processing each FITS file
for file_fits in files_fits:
    # making pathlib object
    path_fits = pathlib.Path (file_fits)
    # check of FITS file
    if not ( (path_fits.exists () ) and (path_fits.is_file () ) \
             and (file_fits[-5:] == '.fits') ):
        print ("Input file must be a FITS file.", file=sys.stderr)
        continue
    # reading FITS file
    (header, image) = read_fits (file_fits)
    # WCS
    wcs = astropy.wcs.WCS (header)
    # calculation of (x, y) from RA and Dec
    (x, y) = wcs.world_to_pixel (coord)
    
    # region of calculation
    box_xmin = int (x) - half_width
    box_xmax = int (x) + half_width
    box_ymin = int (y) - half_width
    box_ymax = int (y) + half_width
    box = image[box_ymin:box_ymax, box_xmin:box_xmax]

    # rough background subtraction
    box -= numpy.median (box)

    # PSF fitting
    box_y, box_x = numpy.indices (box.shape)
    if (psf_model == '2dg'):
        psf_init = astropy.modeling.models.Gaussian2D (x_mean=half_width, \
                                                       y_mean=half_width)
    elif (psf_model == '2dm'):
        psf_init = astropy.modeling.models.Moffat2D (x_0=half_width, \
                                                     y_0=half_width)
    fit = astropy.modeling.fitting.LevMarLSQFitter ()
    psf_fitted = fit (psf_init, box_x, box_y, box, maxiter=1000)

    # results of fitting
    if (psf_model == '2dg'):
        x_centre = psf_fitted.x_mean.value
        y_centre = psf_fitted.y_mean.value
        fwhm_x   = psf_fitted.x_fwhm
        fwhm_y   = psf_fitted.y_fwhm
        fwhm     = (fwhm_x + fwhm_y) / 2.0
    elif (psf_model == '2dm'):
        x_centre = psf_fitted.x_0.value
        y_centre = psf_fitted.y_0.value
        fwhm     = psf_fitted.fwhm

    # (x, y) of centre of star
    x_c = x_centre + box_xmin
    y_c = y_centre + box_ymin

    # measured RA and Dec
    coord_measured = wcs.pixel_to_world (x_c, y_c)

    # printing result
    print ("%-26s %12.8f %+13.8f %-32s # (x,y)=(%10.4f,%10.4f)" \
           % (path_fits.name, coord_measured.ra.deg, coord_measured.dec.deg, \
              coord_measured.to_string ('hmsdms'), x_c, y_c) )
