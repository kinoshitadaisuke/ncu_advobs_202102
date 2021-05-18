#!/usr/pkg/bin/python3.9

# importing argparse module
import argparse

# importing sys module
import sys

# importing astropy module
import astropy.io.fits

# importing matplotlib module
import matplotlib.figure
import matplotlib.backends.backend_agg

# constructing parser object
desc   = "making a histogram"
parser = argparse.ArgumentParser (description=desc)

# adding arguments
parser.add_argument ('-i', '--input', default='', help='input FITS file')
parser.add_argument ('-o', '--output', default='', help='output file')
parser.add_argument ('-a', '--z1', type=float, default=0.0, \
                     help='minimum pixel value (default: 0)')
parser.add_argument ('-b', '--z2', type=float, default=70000.0, \
                     help='maximum pixel value (default: 70000)')
parser.add_argument ('-n', '--nbins', type=int, default=7000, \
                     help='number of bins (default: 7000)')

# parsing arguments
args = parser.parse_args ()

# input parameters
file_input  = args.input
file_output = args.output
z1          = args.z1
z2          = args.z2
nbins       = args.nbins

# check of input FITS file
if not (file_input[-5:] == '.fits'):
    print ("Input file must be a FITS file.")
    print ("Input file name = %s" % file_input)
    sys.exit ()

# check of output file
if not ( (file_output[-4:] == '.eps') or (file_output[-4:] == '.pdf') \
         or (file_output[-4:] == '.png') or (file_output[-3:] == '.ps') ):
    print ("Output file must be either EPS, PDF, PNG, or PS.")
    print ("Output file name = %s" % file_output)
    sys.exit ()

# opening FITS file
with astropy.io.fits.open (file_input) as hdu_list:
    # reading image data
    data = hdu_list[0].data

# flattening data
data_flat = data.flatten ()

# making objects "fig" and "ax"
fig = matplotlib.figure.Figure ()
matplotlib.backends.backend_agg.FigureCanvasAgg (fig)
ax = fig.add_subplot (111)

# axes
ax.set_xlabel ("Pixel Value [ADU]")
ax.set_ylabel ("Number of pixels")

# plotting image
ax.set_xlim (z1, z2)
ax.hist (data_flat, bins=nbins, range=(z1, z2), histtype='bar', align='mid', \
         label=file_input)
ax.legend ()

# saving file
fig.savefig (file_output, dpi=225)
