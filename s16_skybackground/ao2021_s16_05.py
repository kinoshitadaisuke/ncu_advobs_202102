#!/usr/pkg/bin/python3.9

#
# Time-stamp: <2021/05/19 09:12:19 (CST) daisuke>
#

# importing argparse module
import argparse

# importing sys module
import sys

# importing datetime module
import datetime

# importing numpy module
import numpy

# importing astropy module
import astropy.coordinates

# importing photutils module
import photutils.centroids
import photutils.aperture

# constructing parser object
desc = 'aperture photometry of a star at given RA and Dec'
parser = argparse.ArgumentParser (description=desc)

# adding argument
parser.add_argument ('-i', '--input', default='', \
                     help='input FITS file name')
parser.add_argument ('-o', '--output', default='', \
                     help='output file name')
parser.add_argument ('-r', '--ra', type=float, default=-999.999, \
                     help='RA in degree')
parser.add_argument ('-d', '--dec', type=float, default=-999.999, \
                     help='Dec in degree')
parser.add_argument ('-a', '--aperture', type=float, default=1.5, \
                     help='aperture radius in FWHM (default: 1.5)')
parser.add_argument ('-w', '--halfwidth', type=int, default=10, \
                     help='half-width for centroid measurement (default: 10)')
parser.add_argument ('-s1', '--skyannulus1', type=float, default=4.0, \
                     help='inner sky annulus radius in FWHM (default: 4)')
parser.add_argument ('-s2', '--skyannulus2', type=float, default=7.0, \
                     help='outer sky annulus radius in FWHM (default: 7)')
parser.add_argument ('-t', '--threshold', type=float, default=3.0, \
                     help='threshold for sigma-clipping in sigma (default: 3)')
parser.add_argument ('-n', '--maxiters', type=int, default=10, \
                     help='maximum number of iterations (default: 10)')
parser.add_argument ('-e', '--keyword-exptime', default='EXPTIME', \
                     help='FITS keyword for exposure time (default: EXPTIME)')
parser.add_argument ('-f', '--keyword-filter', default='FILTER', \
                     help='FITS keyword for filter name (default: FILTER)')
parser.add_argument ('-m', '--keyword-airmass', default='AIRMASS', \
                     help='FITS keyword for airmass (default: AIRMASS)')

# command-line argument analysis
args = parser.parse_args ()

# input parameters
file_fits             = args.input
file_output           = args.output
target_ra_deg         = args.ra
target_dec_deg        = args.dec
aperture_radius_fwhm  = args.aperture
halfwidth             = args.halfwidth
skyannulus_inner_fwhm = args.skyannulus1
skyannulus_outer_fwhm = args.skyannulus2
threshold             = args.threshold
maxiters              = args.maxiters
keyword_exptime       = args.keyword_exptime
keyword_filter        = args.keyword_filter
keyword_airmass       = args.keyword_airmass

# check of RA and Dec
if ( (target_ra_deg < 0.0) or (target_ra_deg > 360.0) \
     or (target_dec_deg < -90.0) or (target_dec_deg > 90.0) ):
    print ("Something is wrong with RA or Dec!")
    print ("Check RA and Dec you specify.")
    print ("RA  = %f deg" % target_ra_deg)
    print ("Dec = %f deg" % target_dec_deg)
    sys.exit ()

# check of FITS file name
if not (file_fits[-5:] == '.fits'):
    print ("The file \"%s\" is not a FITS file!" % file_fits)
    print ("Check the file name.")
    sys.exit ()

# check of output file
if (file_output == ''):
    print ("Output file name must be given.")
    sys.exit ()

# date/time
now = datetime.datetime.now ().isoformat ()
    
# opening FITS file
with astropy.io.fits.open (file_fits) as hdu_list:
    # reading header information
    header = hdu_list[0].header
    # WCS information
    wcs = astropy.wcs.WCS (header)
    # reading image data
    data   = hdu_list[0].data

# extraction of information from FITS header
exptime      = header[keyword_exptime]
filter_name  = header[keyword_filter]
airmass      = header[keyword_airmass]
    
# sky coordinate
coord_sky = astropy.coordinates.SkyCoord (target_ra_deg, target_dec_deg, \
                                          unit='deg')

# conversion from sky coordinate into pixel coordinate
coord_pix = wcs.world_to_pixel (coord_sky)

# initial guess of target position
init_x = coord_pix[0]
init_y = coord_pix[1]

# region of sub-frame for centroid measurement
subframe_xmin = int (init_x - halfwidth)
subframe_xmax = int (init_x + halfwidth + 1)
subframe_ymin = int (init_y - halfwidth)
subframe_ymax = int (init_y + halfwidth + 1)

# extracting sub-frame for centroid measurement
subframe = data[subframe_ymin:subframe_ymax, subframe_xmin:subframe_xmax]

# sky subtraction
subframe_skysub = subframe - numpy.median (subframe)

# centroid measurement
(com_x, com_y) = photutils.centroids.centroid_com (subframe_skysub)

# PSF fitting
subframe_y, subframe_x = numpy.indices (subframe_skysub.shape)
psf_init = astropy.modeling.models.Gaussian2D (x_mean=com_x, y_mean=com_y)
fit = astropy.modeling.fitting.LevMarLSQFitter ()
psf_fitted = fit (psf_init, subframe_x, subframe_y, subframe_skysub, \
                  maxiter=1000)

# fitted PSF parameters
centre_x = psf_fitted.x_mean.value + subframe_xmin
centre_y = psf_fitted.y_mean.value + subframe_ymin
theta    = psf_fitted.theta.value
fwhm_x   = psf_fitted.x_fwhm
fwhm_y   = psf_fitted.y_fwhm
fwhm     = (fwhm_x + fwhm_y) / 2.0

# position of centre of star in pixel coordinate
position_pix = (centre_x, centre_y)

# aperture radius in pixel
aperture_radius_pix  = fwhm * aperture_radius_fwhm
skyannulus_inner_pix = fwhm * skyannulus_inner_fwhm
skyannulus_outer_pix = fwhm * skyannulus_outer_fwhm

# making aperture
apphot_aperture \
    = photutils.aperture.CircularAperture (position_pix, r=aperture_radius_pix)
apphot_annulus \
    = photutils.aperture.CircularAnnulus (position_pix, \
                                          r_in=skyannulus_inner_pix, \
                                          r_out=skyannulus_outer_pix)

# making masked data for sky annulus
skyannulus_data       = apphot_annulus.to_mask (method='center').multiply (data)
skyannulus_mask       = skyannulus_data <= 0.0
skyannulus_maskeddata = numpy.ma.array (skyannulus_data, mask=skyannulus_mask)

# sky background estimate using sigma-clipping algorithm
skybg_mean, skybg_median, skybg_stddev \
    = astropy.stats.sigma_clipped_stats (skyannulus_maskeddata, \
                                         sigma=threshold, maxiters=maxiters, \
                                         cenfunc='median')
skybg_per_pix = skybg_median

# aperture photometry
noise        = numpy.sqrt (data)
phot_star    = photutils.aperture.aperture_photometry (data, apphot_aperture, \
                                                       error=noise)
net_flux     = phot_star['aperture_sum'] - skybg_per_pix * apphot_aperture.area
net_flux_err = phot_star['aperture_sum_err']

# instrumental magnitude
instmag     = -2.5 * numpy.log10 (net_flux / exptime)
instmag_err = 2.5 / numpy.log (10) * net_flux_err / net_flux

# writing results into a file
with open (file_output, 'w') as fh:
    fh.write ("#\n")
    fh.write ("# Result of Aperture Photometry\n")
    fh.write ("#\n")
    fh.write ("#  Date/Time of Analysis\n")
    fh.write ("#   Date/Time = %s\n" % now)
    fh.write ("#\n")
    fh.write ("#  Input Parameters\n")
    fh.write ("#   FITS file                    = %s\n" % file_fits)
    fh.write ("#   RA                           = %f deg\n" % target_ra_deg)
    fh.write ("#   Dec                          = %f deg\n" % target_dec_deg)
    fh.write ("#   aperture radius              = %f in FWHM\n" \
              % aperture_radius_fwhm)
    fh.write ("#   inner sky annulus            = %f in FWHM\n" \
              % skyannulus_inner_fwhm)
    fh.write ("#   outer sky annulus            = %f in FWHM\n" \
              % skyannulus_outer_fwhm)
    fh.write ("#   half-width for centroid      = %f pixel\n" % halfwidth)
    fh.write ("#   threshold for sigma-clipping = %f in sigma\n" % threshold)
    fh.write ("#   number of max iterations     = %d\n" % maxiters)
    fh.write ("#   keyword for airmass          = %s\n" % keyword_airmass)
    fh.write ("#   keyword for exposure time    = %s\n" % keyword_exptime)
    fh.write ("#   keyword for filter name      = %s\n" % keyword_filter)
    fh.write ("#\n")
    fh.write ("#  Calculated and Measured Quantities\n")
    fh.write ("#   (init_x, init_y)     = (%f, %f)\n" % (init_x, init_y) )
    fh.write ("#   (com_x, com_y)       = (%f, %f)\n" \
              % (com_x + subframe_xmin, com_y + subframe_ymin) )
    fh.write ("#   (centre_x, centre_y) = (%f, %f)\n" % (centre_x, centre_y) )
    fh.write ("#   FWHM of stellar PSF  = %f pixel\n" % fwhm)
    fh.write ("#   sky background level = %f ADU per pixel\n" % skybg_per_pix)
    fh.write ("#   net flux             = %f ADU\n" % net_flux)
    fh.write ("#   net flux err         = %f ADU\n" % net_flux_err)
    fh.write ("#   instrumental mag     = %f\n" % instmag)
    fh.write ("#   instrumental mag err = %f\n" % instmag_err)
    fh.write ("#\n")
    fh.write ("#  Results\n")
    fh.write ("#   file, exptime, filter, centre_x, centre_y,\n"
              "#   net_flux, net_flux_err, instmag, instmag_err, airmass\n")
    fh.write ("%s %f %s %f %f %f %f %f %f %f\n" \
              % (file_fits, exptime, filter_name, centre_x, centre_y, \
                 net_flux, net_flux_err, instmag, instmag_err, airmass) )
