#!/usr/pkg/bin/python3.9

#
# Time-stamp: <2021/05/27 15:22:52 (CST) daisuke>
#

# importing argparse module
import argparse

# importing sys module
import sys

# importing numpy module
import numpy.random

# constructing parser object
desc   = 'generating random (x, y) positions of stars'
parser = argparse.ArgumentParser (description=desc)

# adding arguments
parser.add_argument ('-n', '--number', type=int, default=10, \
                     help='number of stars (default: 10)')
parser.add_argument ('-x', '--size-x', type=int, default=1024, \
                     help='image size on x-axis (default: 2048)')
parser.add_argument ('-y', '--size-y', type=int, default=1024, \
                     help='image size on y-axis (default: 2048)')
parser.add_argument ('-a', '--flux-min', type=float, default=100.0, \
                     help='minimum flux of stars (default: 100)')
parser.add_argument ('-b', '--flux-max', type=float, default=100000.0, \
                     help='maximum flux of stars (default: 100000)')
parser.add_argument ('-o', '--file-output', default='', \
                     help='output file name')

# command-line argument analysis
args = parser.parse_args ()

# parameters
nstars      = args.number
size_x      = args.size_x
size_y      = args.size_y
flux_min    = args.flux_min
flux_max    = args.flux_max
file_output = args.file_output

# check of output file name
if (file_output == ''):
    print ("Output file name must be specified.")
    sys.exit ()

# generating random numbers
position_x = numpy.random.default_rng ().uniform (0, size_x, nstars)
position_y = numpy.random.default_rng ().uniform (0, size_y, nstars)
flux       = numpy.random.default_rng ().uniform (flux_min, flux_max, nstars)

# writing data to file
with open (file_output, 'w') as fh_out:
    # for each object
    for i in range ( len (position_x) ):
        # writing x, y, flux to file
        fh_out.write ("%f %f %f\n" % (position_x[i], position_y[i], flux[i]) )
