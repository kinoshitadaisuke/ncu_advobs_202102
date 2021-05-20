#!/usr/pkg/bin/python3.9

#
# Time-stamp: <2021/05/21 00:58:41 (CST) daisuke>
#

# importing argparse module
import argparse

# importing sys module
import sys

# importing datetime module
import datetime

# importing numpy module
import numpy.ma

# importing astropy module
import astropy.io.fits

# importing photutils module
import photutils.segmentation

# constructing parser object
desc = 'detecting and masking sources'
parser = argparse.ArgumentParser (description=desc)

# adding arguments
parser.add_argument ('-i', '--input-file', default='', \
                     help='input file name')
parser.add_argument ('-o', '--output-file', default='', \
                     help='output file name')
parser.add_argument ('-t', '--threshold', type=float, default=2.0, \
                     help='detection threshold in sigma (default: 2)')
parser.add_argument ('-n', '--npixels', type=int, default=5, \
                     help='minimum number of pixels for detection (default: 5)')
parser.add_argument ('-s', '--dilate-size', type=int, default=21, \
                     help='dilate size (default: 21)')
parser.add_argument ('-m', '--maxiters', type=int, default=30, \
                     help='maximum number of iterations (default: 30)')

# command-line argument analysis
args = parser.parse_args ()

# file names
file_input  = args.input_file
file_output = args.output_file

# input parameters
threshold   = args.threshold
npixels     = args.npixels
dilate_size = args.dilate_size
maxiters    = args.maxiters

# check of input file name
if not (file_input[-5:] == '.fits'):
    print ("Input file must be a FITS file.")
    sys.exit ()

# check of output file name
if not (file_output[-5:] == '.fits'):
    print ("Output file must be a FITS file.")
    sys.exit ()

# opening FITS file
with astropy.io.fits.open (file_input) as hdu:
    # reading header and image
    header = hdu[0].header
    image  = hdu[0].data
    # if no image in PrimaryHDU, then read next HDU
    if (header['NAXIS'] == 0):
        header = hdu[1].header
        image  = hdu[1].data

# making source mask
source_mask \
    = photutils.segmentation.make_source_mask (image, threshold,
                                               npixels=npixels,
                                               sigclip_iters=maxiters,
                                               dilate_size=dilate_size)

# making masked array
image_masked = numpy.ma.array (image, mask=source_mask)

# now
datetime_now = datetime.datetime.now ()

# adding comments in header
header['comment'] = "Updated on %s" % (datetime_now)
header['comment'] = "Detecting and masking sources"
header['comment'] = "Original file = %s" % (file_input)
header['comment'] = "Options:"
header['comment'] = "  threshold   = %f sigma" % (threshold)
header['comment'] = "  npixels     = %d pixels" % (npixels)
header['comment'] = "  dilate_size = %d pixels" % (dilate_size)
header['comment'] = "  maxiters    = %d" % (maxiters)

# writing a FITS file
astropy.io.fits.writeto (file_output, \
                         numpy.ma.filled (image_masked, fill_value=0.0), \
                         header=header)
