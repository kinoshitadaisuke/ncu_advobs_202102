#!/usr/pkg/bin/python3.9

#
# Time-stamp: <2021/06/02 03:35:46 (CST) daisuke>
#

# importing argparse module
import argparse

# importing sys module
import sys

# importing numpy module
import numpy

# importing astropy module
import astropy.io.fits
import astropy.modeling

# importing photutils module
import photutils.aperture

# choice of PSF model
list_psf_model = ['2dg', '2dm']

# constructing parser object
desc = "Aperture photometry"
parser = argparse.ArgumentParser (description=desc)

# adding arguments
parser.add_argument ('-m', '--file-match', default='', \
                     help='star-to-star correspondence file')
parser.add_argument ('-o', '--file-output', default='', \
                     help='output file')
parser.add_argument ('-w', '--width', type=int, default=10, \
                     help='half-width of PSF fitting box (default: 10)')
parser.add_argument ('-p', '--psf-model', choices=list_psf_model, \
                     default='2dg', \
                     help='PSF model [2dg=Gaussian, 2dm=Moffat] (default: 2dg)')
parser.add_argument ('-r', '--radius', type=float, default=1.5, \
                     help='radius of aperture in FWHM (default: 1.5)')
parser.add_argument ('-a', '--sky-annulus-inner', type=float, default=5.0, \
                     help='radius of inner sky annulus in FWHM (default: 4)')
parser.add_argument ('-b', '--sky-annulus-outer', type=float, default=8.0, \
                     help='radius of outer sky annulus in FWHM (default: 6)')
parser.add_argument ('-f', '--fwhm-min', type=float, default=2.0, \
                     help='minimum acceptable FWHM in pixel (default: 2)')
parser.add_argument ('-g', '--fwhm-max', type=float, default=8.0, \
                     help='maximum acceptable FWHM in pixel (default: 8)')
parser.add_argument ('fits1', nargs=1, help='FITS file 1')
parser.add_argument ('fits2', nargs=1, help='FITS file 2')

# command-line argument analysis
args = parser.parse_args ()

# input parameters
file_match             = args.file_match
file_output            = args.file_output
file_fits1             = args.fits1[0]
file_fits2             = args.fits2[0]
half_width             = args.width
psf_model              = args.psf_model
aperture_radius_fwhm   = args.radius
sky_annulus_inner_fwhm = args.sky_annulus_inner
sky_annulus_outer_fwhm = args.sky_annulus_outer
fwhm_min               = args.fwhm_min
fwhm_max               = args.fwhm_max

# check of match file
if (file_match == ''):
    print ("Star-to-star correspondence file must be specified.")
    sys.exit ()

# check of output file
if (file_output == ''):
    print ("Output file must be specified.")
    sys.exit ()

# check of FITS file 1 and 2
if not ( (file_fits1[-5:] == '.fits') and (file_fits2[-5:] == '.fits') ):
    print ("FITS files must have extension \".fits\".")
    sys.exit ()

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

# making empty lists for coordinates of stars
list_coord1 = []
list_coord2 = []
    
# opening match file
with open (file_match, 'r') as fh_in:
    # reading match file line-by-line
    for line in fh_in:
        # if the line starts with '#', then skip
        if (line[0] == '#'):
            continue
        # splitting data
        (x1_str, y1_str, x2_str, y2_str) = line.split ()
        # conversion from string into float
        x1 = float (x1_str)
        y1 = float (y1_str)
        x2 = float (x2_str)
        y2 = float (y2_str)
        # appending data to lists
        list_coord1.append ( (x1, y1) )
        list_coord2.append ( (x2, y2) )

# reading FITS files
(header1, image1) = read_fits (file_fits1)
(header2, image2) = read_fits (file_fits2)

# opening file for writing
fh_out = open (file_output, 'w')

# writing header to file
fh_out.write ("# x1, y1, net_flux1, x2, y2, net_flux2\n")

# aperture photometry
for i in range ( len (list_coord1) ):
    # X and Y coordinates
    (x1, y1) = list_coord1[i]
    (x2, y2) = list_coord2[i]

    if not ( (x1 - half_width > 0.0) \
             and (y1 - half_width > 0.0) \
             and (x1 + half_width < header1['NAXIS1']) \
             and (y1 + half_width < header1['NAXIS1']) \
             and (x2 - half_width > 0.0) \
             and (y2 - half_width > 0.0) \
             and (x2 + half_width < header2['NAXIS1']) \
             and (y2 + half_width < header2['NAXIS1']) ):
        continue
    
    # region of calculation
    box1_xmin = int (x1) - half_width
    box1_xmax = int (x1) + half_width
    box1_ymin = int (y1) - half_width
    box1_ymax = int (y1) + half_width
    box2_xmin = int (x2) - half_width
    box2_xmax = int (x2) + half_width
    box2_ymin = int (y2) - half_width
    box2_ymax = int (y2) + half_width
    box1 = image1[box1_ymin:box1_ymax, box1_xmin:box1_xmax]
    box2 = image2[box2_ymin:box2_ymax, box2_xmin:box2_xmax]

    # rough background subtraction
    box1 -= numpy.median (box1)
    box2 -= numpy.median (box2)
    
    # PSF fitting
    box1_y, box1_x = numpy.indices (box1.shape)
    box2_y, box2_x = numpy.indices (box2.shape)
    if (psf_model == '2dg'):
        psf1_init = astropy.modeling.models.Gaussian2D (x_mean=half_width, \
                                                        y_mean=half_width)
        psf2_init = astropy.modeling.models.Gaussian2D (x_mean=half_width, \
                                                        y_mean=half_width)
    elif (psf_model == '2dm'):
        psf1_init = astropy.modeling.models.Moffat2D (x_0=half_width, \
                                                      y_0=half_width)
        psf2_init = astropy.modeling.models.Moffat2D (x_0=half_width, \
                                                      y_0=half_width)
    fit = astropy.modeling.fitting.LevMarLSQFitter ()
    psf1_fitted = fit (psf1_init, box1_x, box1_y, box1, maxiter=1000)
    psf2_fitted = fit (psf2_init, box2_x, box2_y, box2, maxiter=1000)

    # results of fitting
    if (psf_model == '2dg'):
        x1_centre = psf1_fitted.x_mean.value
        y1_centre = psf1_fitted.y_mean.value
        fwhm_x1   = psf1_fitted.x_fwhm
        fwhm_y1   = psf1_fitted.y_fwhm
        fwhm1     = (fwhm_x1 + fwhm_y1) / 2.0
        x2_centre = psf2_fitted.x_mean.value
        y2_centre = psf2_fitted.y_mean.value
        fwhm_x2   = psf2_fitted.x_fwhm
        fwhm_y2   = psf2_fitted.y_fwhm
        fwhm2     = (fwhm_x2 + fwhm_y2) / 2.0
    elif (psf_model == '2dm'):
        x1_centre = psf1_fitted.x_0.value
        y1_centre = psf1_fitted.y_0.value
        fwhm1     = psf1_fitted.fwhm
        x2_centre = psf2_fitted.x_0.value
        y2_centre = psf2_fitted.y_0.value
        fwhm2     = psf2_fitted.fwhm

    # check of measured FWHM
    if not ( (fwhm1 > fwhm_min) and (fwhm1 < fwhm_max) \
             and (fwhm2 > fwhm_min) and (fwhm2 < fwhm_max) ):
        continue

    # aperture radius in pixel
    aperture_radius1_pix = fwhm1 * aperture_radius_fwhm
    aperture_radius2_pix = fwhm2 * aperture_radius_fwhm

    # sky annulus
    sky_annulus_inner1_pix = fwhm1 * sky_annulus_inner_fwhm
    sky_annulus_outer1_pix = fwhm1 * sky_annulus_outer_fwhm
    sky_annulus_inner2_pix = fwhm2 * sky_annulus_inner_fwhm
    sky_annulus_outer2_pix = fwhm2 * sky_annulus_outer_fwhm

    # aperture
    aperture1 \
        = photutils.aperture.CircularAperture ((x1, y1), r=aperture_radius1_pix)
    aperture2 \
        = photutils.aperture.CircularAperture ((x2, y2), r=aperture_radius2_pix)

    # sky annulus
    annulus1 \
        = photutils.aperture.CircularAnnulus ((x1, y1), \
                                              r_in=sky_annulus_inner1_pix, \
                                              r_out=sky_annulus_outer1_pix)
    annulus2 \
        = photutils.aperture.CircularAnnulus ((x2, y2), \
                                              r_in=sky_annulus_inner2_pix, \
                                              r_out=sky_annulus_outer2_pix)

    # masked data for sky annulus
    skyannulus1_data  = annulus1.to_mask (method='center').multiply (image1)
    skyannulus1_mask  = skyannulus1_data <= 0.0
    skyannulus1_mdata = numpy.ma.array (skyannulus1_data, mask=skyannulus1_mask)
    skyannulus2_data  = annulus2.to_mask (method='center').multiply (image2)
    skyannulus2_mask  = skyannulus2_data <= 0.0
    skyannulus2_mdata = numpy.ma.array (skyannulus2_data, mask=skyannulus2_mask)

    # sky background estimate using sigma-clipping algorithm
    skybg1_mean, skybg1_median, skybg1_stddev \
        = astropy.stats.sigma_clipped_stats (skyannulus1_mdata, \
                                             sigma=3.0, maxiters=10, \
                                             cenfunc='median')
    skybg2_mean, skybg2_median, skybg2_stddev \
        = astropy.stats.sigma_clipped_stats (skyannulus2_mdata, \
                                             sigma=3.0, maxiters=10, \
                                             cenfunc='median')

    # sky background per pixel
    skybg1_per_pixel = 3.0 * skybg1_median - 2.0 * skybg1_mean
    skybg2_per_pixel = 3.0 * skybg2_median - 2.0 * skybg2_mean
    
    # aperture photometry
    phot_star1 = photutils.aperture.aperture_photometry (image1, aperture1)
    phot_star2 = photutils.aperture.aperture_photometry (image2, aperture2)

    # net flux
    net_flux1 = phot_star1['aperture_sum'] - skybg1_per_pixel * aperture1.area
    net_flux2 = phot_star2['aperture_sum'] - skybg2_per_pixel * aperture2.area

    # rejecting some data
    if ( (net_flux1 < 0.0) or (net_flux2 < 0.0) ):
        continue
    if ( (net_flux1 / net_flux2 > 2.0) or (net_flux1 / net_flux2 < 0.5) ):
        continue

    # writing data to file
    fh_out.write ("%10.5f %10.5f %12.4f %10.5f %10.5f %12.4f\n" \
                  % (x1, y1, net_flux1, x2, y2, net_flux2) )

# closing file
fh_out.close ()
