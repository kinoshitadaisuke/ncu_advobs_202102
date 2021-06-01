#!/usr/pkg/bin/python3.9

#
# Time-stamp: <2021/06/02 00:39:15 (CST) daisuke>
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
import astropy.table
import astropy.visualization

# importing scikit-image module
import skimage.transform

# importing astroalign module
import astroalign

# importing matplotlib module
import matplotlib.figure
import matplotlib.backends.backend_agg

# date/time
now = datetime.datetime.now ()

# constructing parser object
desc   = 'Finding star-to-star correspondence'
parser = argparse.ArgumentParser (description=desc)

# adding arguments
parser.add_argument ('-m', '--max-points', type=int, default=100, \
                     help='maximum number of control points (default: 100)')
parser.add_argument ('-o', '--file-output', default='', \
                     help='output graphics file')
parser.add_argument ('-l', '--file-log', default='', \
                     help='output log file')
parser.add_argument ('votable1', nargs=1, help='VOTable file 1')
parser.add_argument ('votable2', nargs=1, help='VOTable file 2')
parser.add_argument ('fits1', nargs=1, help='FITS file 1')
parser.add_argument ('fits2', nargs=1, help='FITS file 2')

# command-line argument analysis
args = parser.parse_args ()

# file names
file_vot1  = args.votable1[0]
file_vot2  = args.votable2[0]
file_fits1 = args.fits1[0]
file_fits2 = args.fits2[0]
file_fig   = args.file_output
file_log   = args.file_log
max_points = args.max_points

# check of output file name
if not ( (file_fig[-4:] == '.eps') or (file_fig[-4:] == '.pdf') \
         or (file_fig[-4:] == '.png') or (file_fig[-3:] == '.ps') ):
    print ("Figure file name must be either EPS, PDF, PNG, or PS.")
    sys.exit ()

# check of VOTable file name
if not ( (file_vot1[-4:] == '.vot') and (file_vot2[-4:] == '.vot') ):
    print ("Input file must be a VOTable file (*.vot).")
    print ("vot file 1 = %s" % file_vot1)
    print ("vot file 2 = %s" % file_vot2)
    sys.exit ()

# check of FITS file name
if not ( (file_fits1[-5:] == '.fits') and (file_fits2[-5:] == '.fits') ):
    print ("Input file must be a FITS file (*.fits).")
    print ("FITS file 1 = %s" % file_fits1)
    print ("FITS file 2 = %s" % file_fits2)
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

# reading VOTable files
votable_source1 = astropy.io.votable.parse (file_vot1)
votable_source2 = astropy.io.votable.parse (file_vot2)

# reading the first table in VOTable
table_source1 = votable_source1.get_first_table ()
table_source2 = votable_source2.get_first_table ()

# (x, y) coordinates of sources
list_source1_x = table_source1.array['xcentroid']
list_source1_y = table_source1.array['ycentroid']
list_source2_x = table_source2.array['xcentroid']
list_source2_y = table_source2.array['ycentroid']
position_1 = numpy.transpose ( (list_source1_x, list_source1_y) )
position_2 = numpy.transpose ( (list_source2_x, list_source2_y) )

# finding star-to-star matching
transf, (list_matched_2, list_matched_1) \
    = astroalign.find_transform (position_2, position_1, \
                                 max_control_points=max_points)

# transformation
list_matched_2_aligned \
    = astroalign.matrix_transform (list_matched_2, transf.params)

# writing results of matching to log file
with open (file_log, 'w') as fh:
    fh.write ("#\n")
    fh.write ("# result of image alignment\n")
    fh.write ("#\n")
    fh.write ("#   date/time = %s\n" % now)
    fh.write ("#\n")
    fh.write ("# input files\n")
    fh.write ("#\n")
    fh.write ("#   VOTable file 1 = %s\n" % file_vot1)
    fh.write ("#   VOTable file 2 = %s\n" % file_vot2)
    fh.write ("#   FITS file 1 = %s\n" % file_fits1)
    fh.write ("#   FITS file 2 = %s\n" % file_fits2)
    fh.write ("#\n")
    fh.write ("# transformation matrix\n")
    fh.write ("#\n")
    fh.write ("# [\n")
    fh.write ("#   [%f, %f, %f],\n" \
              % (transf.params[0][0], transf.params[0][1], \
                 transf.params[0][2]) )
    fh.write ("#   [%f, %f, %f],\n" \
              % (transf.params[1][0], transf.params[1][1], \
                 transf.params[1][2]) )
    fh.write ("#   [%f, %f, %f]\n" \
              % (transf.params[2][0], transf.params[2][1], \
                 transf.params[2][2]) )
    fh.write ("# ]\n")
    fh.write ("#\n")
    fh.write ("#\n")
    fh.write ("# list of matched stars\n")
    fh.write ("#\n")
    for i in range ( len (list_matched_1) ):
        fh.write ("%10.4f %10.4f %10.4f %10.4f\n" \
                  % (list_matched_1[i][0], list_matched_1[i][1], \
                     list_matched_2[i][0], list_matched_2[i][1]) )

# reading FITS files
(header1, image1) = read_fits (file_fits1)
(header2, image2) = read_fits (file_fits2)

# byte swap
# data stored in FITS file is network byte-order (big endian).
# Intel/AMD CPUs use little endian.
image1 = image1.byteswap ().newbyteorder ()
image2 = image2.byteswap ().newbyteorder ()

# aligning 2nd image to 1st image
st = skimage.transform.SimilarityTransform (scale=transf.scale, \
                                            rotation=transf.rotation, \
                                            translation=transf.translation)
image2_aligned = skimage.transform.warp (image2, st.inverse)

# marker and colour names for matplotlib
markers = ['o', 'v', '^', 's', 'p', 'h', '8']
colours = ['maroon', 'red', 'coral', 'bisque', 'orange', \
           'wheat', 'yellow', 'green', 'lime', 'aqua', \
           'skyblue', 'blue', 'indigo', 'violet', 'pink']

# making objects "fig" and "ax"
fig = matplotlib.figure.Figure ()
matplotlib.backends.backend_agg.FigureCanvasAgg (fig)
ax1 = fig.add_subplot (121)
ax2 = fig.add_subplot (122)

# plotting first image
norm1 \
    = astropy.visualization.mpl_normalize.ImageNormalize \
    ( stretch=astropy.visualization.HistEqStretch (image1) )
im1 = ax1.imshow (image1, origin='lower', cmap='bone', norm=norm1)
for i in range ( len (list_matched_1) ):
    i_marker = i % len (markers)
    i_colour = i % len (colours)
    ax1.plot (list_matched_1[i][0], list_matched_1[i][1], \
              marker=markers[i_marker], color=colours[i_colour], \
              markersize=8, fillstyle='none')
ax1.set_title ('First Image')

# plotting second image
norm2 \
    = astropy.visualization.mpl_normalize.ImageNormalize \
    ( stretch=astropy.visualization.HistEqStretch (image2_aligned) )
im2 = ax2.imshow (image2_aligned, origin='lower', cmap='bone', norm=norm2)
for i in range ( len (list_matched_2_aligned) ):
    i_marker = i % len (markers)
    i_colour = i % len (colours)
    ax2.plot (list_matched_2_aligned[i][0], list_matched_2_aligned[i][1], \
              marker=markers[i_marker], color=colours[i_colour], \
              markersize=8, fillstyle='none')
ax2.set_title ('Second Image')

# writing to a file
fig.savefig (file_fig, dpi=225)
