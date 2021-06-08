#!/usr/pkg/bin/python3.9

#
# Time-stamp: <2021/06/08 21:39:35 (CST) daisuke>
#

# importing argparse module
import argparse

# importing sys module
import sys

# importing datetime module
import datetime

# importing numpy
import numpy

# importing astropy module
import astropy.coordinates
import astropy.io.fits
import astropy.units

# date/time
now = datetime.datetime.now ()

# constructing parser object
desc   = 'adding FITS keywords EPOCH and SECPIX to a FITS file'
parser = argparse.ArgumentParser (description=desc)

# adding arguments
parser.add_argument ('-i', '--file-input', default='', \
                     help='input FITS file name')
parser.add_argument ('-o', '--file-output', default='', \
                     help='output FITS file name')
parser.add_argument ('-e1', '--keyword-epoch', default='EPOCH', \
                     help='FITS keyword for epoch (default: EPOCH)')
parser.add_argument ('-e2', '--keyword-equinox', default='EQUINOX', \
                     help='FITS keyword for equinox (default: EQUINOX)')
parser.add_argument ('-f', '--keyword-focallength', default='FOCALLEN', \
                     help='FITS keyword for focal length (default: FOCALLEN)')
parser.add_argument ('-s1', '--keyword-pixelscale1', default='SECPIX1', \
                     help='FITS keyword for pixel scale (default: SECPIX1)')
parser.add_argument ('-s2', '--keyword-pixelscale2', default='SECPIX2', \
                     help='FITS keyword for pixel scale (default: SECPIX2)')
parser.add_argument ('-x', '--keyword-xpixsz', default='XPIXSZ', \
                     help='FITS keyword for pixel size in X (default: XPIXSZ)')
parser.add_argument ('-y', '--keyword-ypixsz', default='YPIXSZ', \
                     help='FITS keyword for pixel size in Y (default: YPIXSZ)')

# command-line argument analysis
args = parser.parse_args ()

# file names
file_input  = args.file_input
file_output = args.file_output

# FITS keywords
keyword_epoch       = args.keyword_epoch
keyword_equinox     = args.keyword_equinox
keyword_focallength = args.keyword_focallength
keyword_secpix1     = args.keyword_pixelscale1
keyword_secpix2     = args.keyword_pixelscale2
keyword_xpixsz      = args.keyword_xpixsz
keyword_ypixsz      = args.keyword_ypixsz

# check of input file name
if not (file_input[-5:] == '.fits'):
    print ("Input file must be a FITS file.")
    sys.exit ()

# check of catalogue file name
if not (file_output[-5:] == '.fits'):
    print ("Output file must be a FITS file.")
    sys.exit ()

# units
u_ha = astropy.units.hourangle
u_deg = astropy.units.deg
 
# function to read a FITS file
def read_fits (file_fits):
    # reading FITS file
    with astropy.io.fits.open (file_fits) as hdu:
        # reading header and image
        header = hdu[0].header
        image  = hdu[0].data
        # if no image in PrimaryHDU, then read next HDU
        if (header['NAXIS'] == 0):
            header = hdu[1].header
            image  = hdu[1].data
    # returning header and image
    return (header, image)

# reading header and image of a FITS file
(header, image) = read_fits (file_input)

# focal length in mm
focallength_mm = header[keyword_focallength]
focallength_m  = focallength_mm * 10**-3

# pixel sizes
pixelsize_x_micron = header[keyword_xpixsz]
pixelsize_y_micron = header[keyword_ypixsz]
pixelsize_x_m      = pixelsize_x_micron * 10**-6
pixelsize_y_m      = pixelsize_y_micron * 10**-6

# pixel scale
pixelscale_x = pixelsize_x_m / focallength_m / numpy.pi * 180.0 * 3600.0 * -1.0
pixelscale_y = pixelsize_y_m / focallength_m / numpy.pi * 180.0 * 3600.0

# CDELT1 and CDELT2
cdelt1 = pixelsize_x_m / focallength_m / numpy.pi * 180.0 * -1.0
cdelt2 = pixelsize_y_m / focallength_m / numpy.pi * 180.0 * -1.0

# modifying format of RA and DEC value
ra_str  = header['RA']
dec_str = header['DEC']
ra_str  = ra_str.replace (' ', ':')
dec_str = dec_str.replace (' ', ':')

# RA and Dec
coord = astropy.coordinates.SkyCoord (ra_str, dec_str, unit=(u_ha, u_deg), \
                                      frame='icrs')
ra_deg  = coord.ra.deg
dec_deg = coord.dec.deg

# modifying format of HA and ST
ha_str = header['HA']
st_str = header['ST']
ha_str = ha_str.replace (' ', ':')
st_str = st_str.replace (' ', ':')
    
# updating FITS header
header['comment'] = "Updated on %s" % (now)
header['comment'] = "Added some keywords necessary for WCSTools"
header['comment'] = "New file: %s" % (file_output)
header[keyword_epoch]   = 2000.0
header[keyword_equinox] = 2000.0
header[keyword_secpix1] = pixelscale_x
header[keyword_secpix2] = pixelscale_y
header['RA']            = ra_str
header['DEC']           = dec_str
header['HA']            = ha_str
header['ST']            = st_str
header['CRVAL1']        = ra_deg
header['CRVAL2']        = dec_deg
header['CRPIX1']        = header['NAXIS1'] / 2.0
header['CRPIX2']        = header['NAXIS2'] / 2.0
header['CROTA1']        = 0.0
header['CROTA2']        = 0.0
header['CDELT1']        = cdelt1
header['CDELT2']        = cdelt2

# writing a FITS file
astropy.io.fits.writeto (file_output, image, header=header)
