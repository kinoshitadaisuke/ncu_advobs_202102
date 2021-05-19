#!/usr/pkg/bin/python3.9

#
# Time-stamp: <2021/05/19 09:12:53 (CST) daisuke>
#

# importing argparse module
import argparse

# importing sys module
import sys

# importing numpy module
import numpy

# constructing parser object
desc   = "calculating sky background brightness"
parser = argparse.ArgumentParser (description=desc)

# adding arguments
parser.add_argument ('-b', '--skybg', type=float, default=-999.99, \
                     help='sky background level in ADU')
parser.add_argument ('-e', '--skybg-err', type=float, default=-999.99, \
                     help='error of sky background level in ADU')
parser.add_argument ('-s', '--star', type=float, default=-999.99, \
                     help='net flux of standard star in ADU')
parser.add_argument ('-m', '--mag', type=float, default=-999.99, \
                     help='apparent magnitude of standard star')
parser.add_argument ('-p', '--pixelscale', type=float, default=-999.99, \
                     help='pixel scale in arcsec/pix')

# parsing arguments
args = parser.parse_args ()

# input parameters
skybg_flux = args.skybg
skybg_err  = args.skybg_err
star_flux  = args.star
star_mag   = args.mag
pixelscale = args.pixelscale

# check of input parameters
if (skybg_flux < 0.0):
    print ("Something is wrong with sky background level in ADU.")
    sys.exit ()
if (skybg_err < 0.0):
    print ("Something is wrong with error of sky background level in ADU.")
    sys.exit ()
if (star_flux < 0.0):
    print ("Something is wrong with net flux of star in ADU.")
    sys.exit ()
if (star_mag < 0.0):
    print ("Something is wrong with magnitude of star.")
    sys.exit ()
if (pixelscale < 0.0):
    print ("Something is wrong with pixel scale in arcsec/pix.")
    sys.exit ()

# calculation of sky background brightness in mag/arcsec^2
skybg_flux_per_sqarcsec = skybg_flux / pixelscale**2
skybg_mag = star_mag - 2.5 * numpy.log10 (skybg_flux_per_sqarcsec / star_flux)
skybg_err = 2.5 / numpy.log (10) * skybg_err / skybg_flux

# printing result
print ("sky background brightness = %6.3f +/- %6.3f mag/arcsec^2" \
       % (skybg_mag, skybg_err) )
