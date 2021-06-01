#!/usr/pkg/bin/python3.9

#
# Time-stamp: <2021/06/01 23:35:42 (CST) daisuke>
#

# importing argparse module
import argparse

# importing sys module
import sys

# importing numpy module
import numpy

# importing astropy module
import astropy.io.fits
import astropy.io.votable
import astropy.visualization
import astropy.wcs

# importing matplotlib module
import matplotlib.figure
import matplotlib.backends.backend_agg

# constructing parser object
desc = "Reading a VOTable file and plotting locations of extracted sources"
parser = argparse.ArgumentParser (description=desc)

# list of colours
list_colours = ['maroon', 'red', 'coral', 'bisque', 'orange', \
                'wheat', 'yellow', 'green', 'lime', 'aqua', \
                'skyblue', 'blue', 'indigo', 'violet', 'pink']

# list of markers
list_markers = ['o', 'v', '^', 's', 'p', 'h', '8']

# list of cmap
list_cmap = [ 'viridis', 'plasma', 'inferno', 'magma', 'cividis', 'gray', \
              'bone', 'pink', 'spring', 'summer', 'autumn', 'winter', \
              'cool', 'hot', 'copper', 'hsv', 'ocean', 'terrain', 'gnuplot', \
              'rainbow', 'turbo'
             ]

# adding arguments
parser.add_argument ('-f', '--file-fits', default='', \
                     help='input FITS file (*.fits)')
parser.add_argument ('-t', '--file-votable', default='', \
                     help='input VOTable file (*.vot)')
parser.add_argument ('-o', '--file-output', default='', \
                     help='output graphics file (EPS, PDF, PNG, or PS file)')
parser.add_argument ('-b', '--cmap', choices=list_cmap, default='bone', \
                     help='cmap for colour bar (default: bone)')
parser.add_argument ('-c', '--colour', choices=list_colours, default='red', \
                     help='colour of marker (default: red)')
parser.add_argument ('-m', '--marker', choices=list_markers, default='o', \
                     help='shape of marker (default: o)')
parser.add_argument ('-r', '--radius', type=float, default=10.0, \
                     help='radius of marker in pixel (default: 10)')

# command-line argument analysis
args = parser.parse_args ()

# input parameters
file_fits    = args.file_fits
file_votable = args.file_votable
file_output  = args.file_output
cmap         = args.cmap
colour       = args.colour
marker       = args.marker
radius       = args.radius

# check of input file
if not (file_fits[-5:] == '.fits'):
    print ("Input file given by -f option must be a FITS file.")
    sys.exit ()
if not (file_votable[-4:] == '.vot'):
    print ("Input file given by -t option must be a VOTable file.")
    sys.exit ()

# check of output file
if not ( (file_output[-4:] == '.eps') or (file_output[-4:] == '.pdf') \
         or (file_output[-4:] == '.png') or (file_output[-3:] == '.ps') ):
    print ("Output file must be either EPS, PDF, PNG, or PS file.")
    sys.exit ()

# function to read header and image of a FITS file
def read_fits (file_fits):
    # opening FITS file
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

# reading VOTable file
votable_sources = astropy.io.votable.parse (file_votable)

# reading the first table in VOTable
table_sources = votable_sources.get_first_table ()

# X and Y coordinates of sources
centroid_x = table_sources.array['xcentroid']
centroid_y = table_sources.array['ycentroid']

# reading FITS file
header, image = read_fits (file_fits)
wcs           = astropy.wcs.WCS (header)

# making objects "fig" and "ax"
fig = matplotlib.figure.Figure ()
matplotlib.backends.backend_agg.FigureCanvasAgg (fig)
ax = fig.add_subplot (111, projection=wcs)

# plotting image
norm \
    = astropy.visualization.mpl_normalize.ImageNormalize \
    ( stretch=astropy.visualization.HistEqStretch (image) )
im = ax.imshow (image, cmap=cmap, norm=norm)
ax.plot (centroid_x, centroid_y, marker=marker, color=colour, \
         markersize=radius, fillstyle='none', linestyle='None')
ax.set_xlabel ('Right Ascension (J2000)')
ax.set_ylabel ('Declination (J2000)')

# saving to file
fig.savefig (file_output, dpi=225)
