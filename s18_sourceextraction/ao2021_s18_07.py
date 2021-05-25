#!/usr/pkg/bin/python3.9

#
# Time-stamp: <2021/05/26 00:11:29 (CST) daisuke>
#

# importing argparse module
import argparse

# importing sys module
import sys

# importing astropy module
import astropy.table
import astropy.visualization

# importing photutils module
import photutils.aperture

# importing matplotlib module
import matplotlib.pyplot

# constructing parser object
desc   = 'reading table from a file'
parser = argparse.ArgumentParser (description=desc)

# adding arguments
parser.add_argument ('-c', '--catalogue-file', default='', \
                     help='input catalogue file name')
parser.add_argument ('-i', '--input-file', default='', \
                     help='input FITS file name')
parser.add_argument ('-o', '--output-file', default='', \
                     help='output file name')
parser.add_argument ('-r', '--radius', type=float, default=10.0, \
                     help='radius of aperture in pixel (default: 10)')

# command-line argument analysis
args = parser.parse_args ()

# catalogue file name
file_catalogue = args.catalogue_file
file_input     = args.input_file
file_output    = args.output_file
radius         = args.radius

# check of catalogue file name
if (file_catalogue == ''):
    print ("Catalogue file name must be specified.")
    sys.exit ()

# check of input FITS file name
if not (file_input[-5:] == '.fits'):
    print ("Input file must be a FITS file.")
    sys.exit ()

# check of output file name
if not ( (file_output[-4:] == '.eps') or (file_output[-4:] == '.pdf') \
         or (file_output[-4:] == '.png') or (file_output[-3:] == '.ps')):
    print ("Output file must be either EPS, PDF, PNG, or PS.")
    sys.exit ()

# reading catalogue from a file
table_source = astropy.table.Table.read (file_catalogue, \
                                         format='ascii.commented_header')

# positions of apertures
list_x = list (table_source['xcentroid'])
list_y = list (table_source['ycentroid'])
positions = []
for i in range ( len (list_x) ):
    positions.append ( (list_x[i], list_y[i]) )

# apertures
apertures = photutils.aperture.CircularAperture (positions, r=radius)

# opening FITS file
with astropy.io.fits.open (file_input) as hdu:
    # reading header and image
    header = hdu[0].header
    image  = hdu[0].data
    # if no image in PrimaryHDU, then read next HDU
    if (header['NAXIS'] == 0):
        header = hdu[1].header
        image  = hdu[1].data

# matplotlib
norm \
    = astropy.visualization.mpl_normalize.ImageNormalize \
    ( stretch=astropy.visualization.HistEqStretch (image) )
matplotlib.pyplot.imshow (image, origin='upper', cmap='viridis', norm=norm)
apertures.plot (color='red', lw=1.0, alpha=0.5)
matplotlib.pyplot.title ("Detected sources")
matplotlib.pyplot.savefig (file_output, dpi=225)
