#!/usr/pkg/bin/python3.9

#
# Time-stamp: <2021/05/19 09:12:42 (CST) daisuke>
#

# importing argparse module
import argparse

# importing sys module
import sys

# importing astropy module
import astropy.io.fits
import astropy.stats

# constructing parser object
desc   = "estimating mode"
parser = argparse.ArgumentParser (description=desc)

# adding arguments
parser.add_argument ('-i', '--input', default='', help='input FITS file')
parser.add_argument ('-t', '--threshold', type=float, default=3.0, \
                     help='threshold for sigma-clipping in sigma (default: 3)')
parser.add_argument ('-n', '--maxiters', type=int, default=10, \
                     help='number of maximum iterations (default: 10)')

# parsing arguments
args = parser.parse_args ()

# input parameters
file_input = args.input
threshold  = args.threshold
maxiters   = args.maxiters

# check of input FITS file
if not (file_input[-5:] == '.fits'):
    print ("Input file must be a FITS file.")
    print ("Input file name = %s" % file_input)
    sys.exit ()

# opening FITS file
with astropy.io.fits.open (file_input) as hdu_list:
    # reading image data
    data = hdu_list[0].data

# calculation of mean and median using sigma-clipping
mean, median, stddev \
    = astropy.stats.sigma_clipped_stats (data, sigma=threshold, \
                                         maxiters=maxiters, cenfunc='median')

# calculation of mode using empirical formula
mode = 3 * median - 2 * mean

print ("mean                         = %10.3f ADU" % mean)
print ("median                       = %10.3f ADU" % median)
print ("mode = 3 * median - 2 * mean = %10.3f ADU" % mode)
print ("stddev                       = %10.3f ADU" % stddev)
