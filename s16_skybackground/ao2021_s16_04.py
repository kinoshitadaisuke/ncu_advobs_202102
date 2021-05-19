#!/usr/pkg/bin/python3.9

#
# Time-stamp: <2021/05/19 09:12:13 (CST) daisuke>
#

# importing argparse module
import argparse

# importing sys module
import sys

# importing numpy module
import numpy

# importing astropy module
import astropy.coordinates
import astropy.units

# importing photutils module
import photutils.aperture

# importing matplotlib module
import matplotlib.pyplot

# constructing parser object
desc = 'marking target objects'
parser = argparse.ArgumentParser (description=desc)

# adding argument
parser.add_argument ('-i', '--input', default='', \
                     help='input FITS file name')
parser.add_argument ('-o', '--output', default='', \
                     help='output image file name')
parser.add_argument ('-c', '--catalogue', default='', \
                     help='standard star catalogue file name')
parser.add_argument ('-r', '--radius', type=float, default=10.0, \
                     help='radius of aperture in arcsec (default: 10)')
parser.add_argument ('-a', '--z1', type=float, default=0.0, \
                     help='minimum value to display (default: 0)')
parser.add_argument ('-b', '--z2', type=float, default=100.0, \
                     help='maximum value to display (default: 100)')

# command-line argument analysis
args = parser.parse_args ()

# input parameters
file_fits      = args.input
file_output    = args.output
file_cat       = args.catalogue
radius_arcsec  = args.radius
z1             = args.z1
z2             = args.z2

# unit
u_arcsec = astropy.units.arcsec

# check of FITS file name
if not (file_fits[-5:] == '.fits'):
    print ("The file \"%s\" is not a FITS file!" % file_fits)
    print ("Check the file name.")
    sys.exit ()

# check of output image file
if not ( (file_output[-4:] == '.eps') or (file_output[-4:] == '.pdf') \
         or (file_output[-4:] == '.png') or (file_output[-3:] == '.ps') ):
    print ("Output image file must be either EPS, PDF, PNG, or PS.")
    print ("Given output image file name = %s" % file_output)
    sys.exit ()

# check of catalogue file
if not (file_cat[-13:] == '.txt.clean.v3'):
    print ("Catalogue file must be SDSS photometric standard star catalogue.")
    print ("Given catalogue file name = %s" % file_cat)
    sys.exit ()

# making an empty dictionary to store data for standards
dic_stds = {}

# opening catalogue file
with open (file_cat, 'r') as fh:
    # reading the file line-by-line
    for line in fh:
        # if the line stars with '#', then skip
        if (line[0] == '#'):
            continue
        # splitting the data
        records = line.split ()
        # star ID
        star_id = int (records[0])
        # RA
        ra_deg  = float (records[1])
        # Dec
        dec_deg = float (records[2])
        # adding data to the dictionary
        if not (star_id in dic_stds):
            dic_stds[star_id] = {}
            dic_stds[star_id]['RA_deg'] = ra_deg
            dic_stds[star_id]['Dec_deg'] = dec_deg

# positions of standard stars
positions = []
for star in dic_stds.keys ():
    positions.append ( (dic_stds[star]['RA_deg'], dic_stds[star]['Dec_deg']) )
coords = astropy.coordinates.SkyCoord (positions, unit='deg')
            
# opening FITS file
with astropy.io.fits.open (file_fits) as hdu_list:
    # reading header information
    header = hdu_list[0].header
    # WCS information
    wcs = astropy.wcs.WCS (header)
    # reading image data
    data   = hdu_list[0].data

# making apertures
apertures_sky \
    = photutils.aperture.SkyCircularAperture (coords, \
                                              r=radius_arcsec * u_arcsec)
apertures_pix = apertures_sky.to_pixel (wcs)
    
# making plot
matplotlib.pyplot.imshow (data, origin='upper', vmin=z1, vmax=z2)
apertures_pix.plot (color='red', lw=1.0)
matplotlib.pyplot.colorbar ()
matplotlib.pyplot.savefig (file_output, dpi=225)
