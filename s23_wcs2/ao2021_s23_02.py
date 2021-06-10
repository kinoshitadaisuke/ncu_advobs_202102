#!/usr/pkg/bin/python3.9

#
# Time-stamp: <2021/06/10 20:23:38 (CST) daisuke>
#

# importing argparse module
import argparse

# importing datetime module
import datetime

# importing pathlib module
import pathlib

# importing subprocess module
import subprocess

# importing numpy module
import numpy
import numpy.ma

# importing astropy module
import astropy.convolution
import astropy.coordinates
import astropy.io.fits
import astropy.io.votable
import astropy.units
import astropy.stats

# importing astroquery module
import astroquery.gaia

# importing photutils module
import photutils.segmentation

# importing ssl module
import ssl

# allow insecure downloading
ssl._create_default_https_context = ssl._create_unverified_context

# units
u_ha  = astropy.units.hourangle
u_deg = astropy.units.deg

# list of astrometric catalogues
list_astrometric_catalogues = ['gaia_edr3', 'gaia_dr2']

# list of sexagesimal format keywords
list_keyword_sexagesimal = ['HA', 'UT', 'ST', 'RA', 'DEC', \
                            'OBJCTRA', 'OBJCTDEC', 'SITELAT', 'SITELONG']

# list of epoch related keywords
list_keyword_epoch = ['EPOCH', 'EQUINOX']

# list of keywords for pixel scale
list_keyword_pixelscale = ['SECPIX1', 'SECPIX2']

# list of keywords for fundamental WCS
list_keyword_wcs = ['CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', \
                    'CDELT1', 'CDELT2', 'CROTA1', 'CROTA2']

# list of keywords for CD matrix
list_keyword_cdmatrix = ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']

# help messages
help_focallength \
    = 'FITS keyword for focal length in mm (default: FOCALLEN)'
help_pixelsize_x \
    = 'FITS keyword for pixel size in X in micron (default: XPIXSZ)'
help_pixelsize_y \
    = 'FITS keyword for pixel size in Y in micron (default: XPIXSZ)'
help_exptime \
    = 'FITS keyword for integration time in sec (default: EXPTIME)'
help_imwcs \
    = 'location of the command "imwcs" (default: /usr/local/bin/imwcs)'

# constructing parser object
desc = "setting WCS in FITS files"
parser = argparse.ArgumentParser (description=desc)

# adding arguments
parser.add_argument ('-a', '--astrometric-catalogue', \
                     choices=list_astrometric_catalogues, \
                     default='gaia_edr3', \
                     help='choice of astrometric catalogue (default: gaia_edr3')
parser.add_argument ('-e', '--epoch', type=float, default=2000.0, \
                     help='epoch (default: 2000.0)')
parser.add_argument ('-r', '--radius', type=float, default=3.0, \
                     help='search radius for catalogue in FOV (default: 3)')
parser.add_argument ('-s', '--threshold-rejection-sigma', \
                     type=float, default=4.0, \
                     help='sigma-clipping threshold in sigma (default: 4)')
parser.add_argument ('-t', '--threshold-detection-obj', \
                     type=float, default=5.0, \
                     help='detection threshold in sigma (default: 5)')
parser.add_argument ('-u', '--threshold-detection-sky', \
                     type=float, default=2.0, \
                     help='detection threshold for sky estimate (default: 2)')
parser.add_argument ('-n', '--npixels', type=int, default=5, \
                     help='minimum number of pixels for detection (default: 5)')
parser.add_argument ('-d', '--dilate-size', type=int, default=21, \
                     help='dilate size (default: 21)')
parser.add_argument ('-m', '--maxiters', type=int, default=30, \
                     help='maximum number of iterations (default: 30)')
parser.add_argument ('-k', '--gaussian-fwhm', type=float, default=3.0, \
                     help='Gaussian FWHM in pixel for convolution (default: 3)')
parser.add_argument ('-b', '--kernel-size', type=int, default=3, \
                     help='Gaussian kernel array size in pixel (default: 3)')
parser.add_argument ('-w', '--command-imwcs', default='/usr/local/bin/imwcs', \
                     help=help_imwcs)
parser.add_argument ('--keyword-exptime', default='EXPTIME', \
                     help=help_exptime)
parser.add_argument ('--keyword-focallength', default='FOCALLEN', \
                     help=help_focallength)
parser.add_argument ('--keyword-pixelsize-x', default='XPIXSZ', \
                     help=help_pixelsize_x)
parser.add_argument ('--keyword-pixelsize-y', default='YPIXSZ', \
                     help=help_pixelsize_y)
parser.add_argument ('FITS_files', nargs='+', help='FITS files')

# command-line argument analysis
args = parser.parse_args ()

# input parameters
astrometric_catalogue   = args.astrometric_catalogue
epoch                   = args.epoch
radius_fov              = args.radius
threshold_detection_obj = args.threshold_detection_obj
threshold_detection_sky = args.threshold_detection_sky
npixels                 = args.npixels
dilate_size             = args.dilate_size
maxiters                = args.maxiters
threshold_rejection     = args.threshold_rejection_sigma
gaussian_fwhm           = args.gaussian_fwhm
kernel_array_size       = args.kernel_size
keyword_exptime         = args.keyword_exptime
keyword_focallength     = args.keyword_focallength
keyword_pixelsize_x     = args.keyword_pixelsize_x
keyword_pixelsize_y     = args.keyword_pixelsize_y
command_imwcs           = args.command_imwcs
files_fits              = args.FITS_files

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

# function to fetch a subset of astrometric catalogue
def fetch_astrometric_catalogue (catalogue, ra_deg, dec_deg, radius_deg, \
                                 file_output):
    # command for database query
    if (catalogue == 'gaia_edr3'):
        query_p \
            = "POINT('ICRS',gaiaedr3.gaia_source.ra,gaiaedr3.gaia_source.dec)"
        query_c = "CIRCLE('ICRS',%f,%f,%f)" % (ra_deg, dec_deg, radius_deg)
        query \
            = "SELECT * from gaiaedr3.gaia_source WHERE CONTAINS(%s,%s)=1;" \
            % (query_p, query_c)
    elif (catalogue == 'gaia_dr2'):
        query_p \
            = "POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec)"
        query_c = "CIRCLE('ICRS',%f,%f,%f)" % (ra_deg, dec_deg, radius_deg)
        query \
            = "SELECT * from gaiadr2.gaia_source WHERE CONTAINS(%s,%s)=1;" \
            % (query_p, query_c)
    # sending a job to Gaia database
    job = astroquery.gaia.Gaia.launch_job_async (query, dump_to_file=True, \
                                                 output_file=file_output, \
                                                 output_format='votable')
    # returning job object
    return (job)

# function to convert vot file to cat file
def convert_from_vot_to_cat (file_vot, file_cat):
    # reading VOTable file
    table_cat = astropy.io.votable.parse_single_table (file_vot).to_table ()
    
    # source_id, ra, dec, and r'-band mag
    data_id  = numpy.array (table_cat['source_id'])
    data_ra  = numpy.array (table_cat['ra'])
    data_dec = numpy.array (table_cat['dec'])
    data_r   = numpy.array (table_cat['phot_rp_mean_mag'])
    
    # writing data of stars in astrometric catalogue into *.cat file
    with open (file_cat, 'w') as fh_cat:
        # writing a header of a local catalogue for WCSTools
        fh_cat.write ("gaia/d/j\n")
        fh_cat.write ("a tiny portion of astrometric catalogue %s\n" \
                      % astrometric_catalogue)
        # writing information of each star to file
        for i in range ( len (data_id) ):
            if ( numpy.isnan (data_r[i]) ):
                # if there is no information of r'-band mag, then use 99.99
                fh_cat.write ("%s %12.8f %+13.8f %7.3f\n" \
                              % (data_id[i], data_ra[i], \
                                 data_dec[i], 99.99 ) )
            else:
                # if there is information of mag, then use mag
                fh_cat.write ("%s %12.8f %+13.8f %7.3f\n" \
                              % (data_id[i], data_ra[i], \
                                 data_dec[i], data_r[i]) )

# function to calculate sky background level
def calc_skybg (image, threshold_detection_sky, npixels, dilate_size, \
                maxiters, threshold_rejection):
    # masking stars
    source_mask \
        = photutils.segmentation.make_source_mask (image, \
                                                   threshold_detection_sky, \
                                                   npixels=npixels, \
                                                   sigclip_iters=maxiters, \
                                                   dilate_size=dilate_size)
    # making a masked array
    image_masked = numpy.ma.array (image, mask=source_mask)
    # sky background estimate using sigma-clipping
    skybg_mean, skybg_median, skybg_stddev \
        = astropy.stats.sigma_clipped_stats (image_masked, \
                                             sigma=threshold_rejection)
    # mode calculation using empirical formula
    skybg_mode = 3.0 * skybg_median - 2.0 * skybg_mean
    # returning estimated sky background level and its standard deviation
    return (skybg_mode, skybg_stddev)

def do_source_extraction (image, skybg_mode, skybg_stddev, \
                          threshold_detection_obj, \
                          gaussian_fwhm, kernel_array_size):
    # detection threshold in ADU
    threshold_detection_adu \
        = skybg_mode + threshold_detection_obj * skybg_stddev
    # 2D Gaussian kernel for convolution
    gaussian_sigma = gaussian_fwhm * astropy.stats.gaussian_fwhm_to_sigma
    kernel = astropy.convolution.Gaussian2DKernel (gaussian_sigma, \
                                                   x_size=kernel_array_size, \
                                                   y_size=kernel_array_size)
    kernel.normalize ()
    # source detection
    image_segm \
        = photutils.segmentation.detect_sources (image, \
                                                 threshold_detection_adu, \
                                                 npixels=npixels, \
                                                 filter_kernel=kernel)
    # deblending
    image_deblend \
        = photutils.segmentation.deblend_sources (image, image_segm, \
                                                  npixels=npixels, \
                                                  filter_kernel=kernel, \
                                                  nlevels=32, contrast=0.001)
    # making a source catalogue
    catalogue_sources \
        = photutils.segmentation.SourceCatalog (image, image_deblend)
    # making a table
    table_sources = catalogue_sources.to_table ()
    # returning source table
    return (table_sources)

#
# main routine
#

# check of command "imwcs"
path_imwcs = pathlib.Path (command_imwcs)
if not ( path_imwcs.exists () and path_imwcs.is_file () ):
    print ("Cannot find the command \"imwcs\"!")
    sys.exit ()

# processing each FITS file
for file_fits in files_fits:
    # making pathlib object
    path_fits = pathlib.Path (file_fits)
    
    # file names
    file_wcs = path_fits.stem + 'w.fits'
    file_log = path_fits.stem + 'w.log'
    file_vot = path_fits.stem + 'w.vot.gz'
    file_cat = path_fits.stem + 'w.cat'
    file_xym = path_fits.stem + 'w.xym'

    print ('FITS file =', file_fits)
    print ('wcs file  =', file_wcs)
    print ('log file  =', file_log)
    print ('vot file  =', file_vot)
    print ('cat file  =', file_cat)
    print ('xym file  =', file_xym)
    
    # check of FITS file name
    if not ( path_fits.exists () and path_fits.is_file () \
             and (file_fits[-5:] == '.fits') ):
        # writing an error message to log file
        with open (file_log, 'w') as fh_log:
            fh_log.write ("Input file \"%s\" is missing, " % file_fits)
            fh_log.write ("or is not a regular file,\n")
            fh_log.write ("or is not a FITS file (*.fits)\n")
            fh_log.write ("Check the file.\n")
        # skip processing the file
        continue
        
    # reading FITS file
    (header, image) = read_fits (file_fits)

    # important parameters of FITS file
    exptime = header[keyword_exptime]
    
    # reformatting keyword values
    #   e.g. '12 34 56.7' ==> '12:34:56.7'
    for keyword in list_keyword_sexagesimal:
        # if the keyword exists, then process
        if (keyword in header):
            # value
            value = header[keyword]
            # if value is not '12:34:56.7' format, then reformat the string
            if not (':' in value):
                # modifying the keyword value
                header[keyword] = value.replace (' ', ':')

    # calculation of pixel scale
    pixelsize_x_micron  = header[keyword_pixelsize_x]
    pixelsize_y_micron  = header[keyword_pixelsize_y]
    focallength_mm      = header[keyword_focallength]
    pixelsize_x_m       = pixelsize_x_micron * 10**-6
    pixelsize_y_m       = pixelsize_y_micron * 10**-6
    focallength_m       = focallength_mm * 10**-3
    pixelscale_x_deg    = pixelsize_x_m / focallength_m * 180.0 / numpy.pi
    pixelscale_y_deg    = pixelsize_y_m / focallength_m * 180.0 / numpy.pi
    pixelscale_x_arcsec = pixelscale_x_deg * 3600.0
    pixelscale_y_arcsec = pixelscale_y_deg * 3600.0

    # adding SECPIX1 and SECPIX2 keywords to FITS header
    for keyword in list_keyword_pixelscale:
        # SECPIX1
        if (keyword[-1] == '1'):
            header[keyword] = pixelscale_x_arcsec * -1.0
        # SECPIX2
        elif (keyword[-1] == '2'):
            header[keyword] = pixelscale_y_arcsec
            
    # adding EPOCH and EQUINOX keywords to FITS header
    for keyword in list_keyword_epoch:
        # header['EPOCH'] = header['EQUINOX'] = epoch (e.g. 2000.0)
        header[keyword] = epoch

    # number of pixels in X and Y
    if ('NAXIS1' in header):
        npix_x = header['NAXIS1']
    else:
        # writing an error message to log file
        with open (file_log, 'w') as fh_log:
            fh_log.write ("NAXIS1 keyword is missing in %s." % file_fits)
            fh_log.write ("Check the file.\n")
        # skip processing the file
        continue
    if ('NAXIS2' in header):
        npix_y = header['NAXIS2']
    else:
        # writing an error message to log file
        with open (file_log, 'w') as fh_log:
            fh_log.write ("NAXIS2 keyword is missing in %s." % file_fits)
            fh_log.write ("Check the file.\n")
        # skip processing the file
        continue

    # field-of-view
    fov_x_arcsec = pixelscale_x_arcsec * npix_x
    fov_y_arcsec = pixelscale_y_arcsec * npix_y
    if (fov_x_arcsec > fov_y_arcsec):
        fov_arcsec = fov_x_arcsec
    else:
        fov_arcsec = fov_y_arcsec
    fov_deg = fov_arcsec / 3600.0
        
    # RA and Dec
    ra_str  = header['RA']
    dec_str = header['DEC']

    # making SkyCoord object
    coord = astropy.coordinates.SkyCoord (ra_str, dec_str, \
                                          unit=(u_ha, u_deg), frame='icrs')
    coord_ra_deg  = coord.ra.deg
    coord_dec_deg = coord.dec.deg

    # adding initial guess of WCS related keywords
    # a flag for WCS
    have_wcs = 'YES'
    # check of CD matrix keywords
    for keyword in list_keyword_cdmatrix:
        # if there is no CD matrix keyword, then set have_wcs = 'NO'
        if not (keyword in header):
            have_wcs = 'NO'
    # if there is no WCS in FITS header, then add initial guess values
    if (have_wcs == 'NO'):
        # adding WCS related keywords in FITS header
        for keyword in list_keyword_wcs:
            if (keyword == 'CRVAL1'):
                header[keyword] = coord_ra_deg
            elif (keyword == 'CRVAL2'):
                header[keyword] = coord_dec_deg
            elif (keyword == 'CRPIX1'):
                header[keyword] = npix_x / 2.0
            elif (keyword == 'CRPIX2'):
                header[keyword] = npix_y / 2.0
            elif (keyword == 'CDELT1'):
                header[keyword] = pixelscale_x_deg * -1.0
            elif (keyword == 'CDELT2'):
                header[keyword] = pixelscale_y_deg * -1.0
            elif (keyword == 'CROTA1'):
                header[keyword] = 0.0
            elif (keyword == 'CROTA2'):
                header[keyword] = 0.0
    
    # writing a FITS file
    now = datetime.datetime.now ()
    header['comment'] = " "
    header['comment'] = "Update on %s." % now
    header['comment'] = "Some keywords needed for WCS are added."
    header['comment'] = "New file name = %s" % file_wcs
    astropy.io.fits.writeto (file_wcs, image, header=header)

    # retrieving a subset of astrometric catalogue
    radius_deg = fov_deg * radius_fov
    job = fetch_astrometric_catalogue (astrometric_catalogue, \
                                       coord_ra_deg, coord_dec_deg, \
                                       radius_deg, file_vot)

    # reading VOTable file and write *.cat file for WCSTools
    convert_from_vot_to_cat (file_vot, file_cat)

    # sky background level estimate
    skybg_mode, skybg_stddev = calc_skybg (image, threshold_detection_sky, \
                                           npixels, dilate_size, maxiters, \
                                           threshold_rejection)

    # source extraction
    table_sources = do_source_extraction (image, skybg_mode, skybg_stddev, \
                                          threshold_detection_obj, \
                                          gaussian_fwhm, kernel_array_size)

    # writing daofind format file of detected sources
    with open (file_xym, 'w') as fh_xym:
        # for each object
        for i in range ( len (table_sources) ):
            # X and Y coordinate
            x_c = table_sources['xcentroid'][i]
            y_c = table_sources['ycentroid'][i]
            # calculation of instrumental magnitude
            net_flux = table_sources['kron_flux'][i]
            inst_mag = -2.5 * numpy.log10 (net_flux / exptime)
            # writing X, Y, and instrumental magnitude
            fh_xym.write ("%10.3f %10.3f %+7.3f\n" % (x_c, y_c, inst_mag) )

    
    # setting WCS of a FITS file using WCSTools
    # number of execution of imwcs command
    n_exe = 5
    # initial tolerance of fitting in pixel
    tolerance_ini = 9
    # iterations
    for i in range (n_exe):
        # tolerance of fitting in pixel
        tolerance = tolerance_ini - 2 * i
        # command
        command_wcs = "%s -v -o -h 200 -t %d -d %s -c %s %s" \
            % (command_imwcs, tolerance, file_xym, file_cat, file_wcs)
        # executing command
        subprocess.run (command_wcs, shell=True)
