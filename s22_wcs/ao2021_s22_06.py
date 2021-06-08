#!/usr/pkg/bin/python3.9

#
# Time-stamp: <2021/06/08 19:12:31 (CST) daisuke>
#

# importing argparse module
import argparse

# importing sys module
import sys

# importing astropy module
import astropy.io.fits

# command-line argument analysis
desc = 'printing header of a FITS file'
parser = argparse.ArgumentParser (description=desc)
parser.add_argument ('-i', '--file-input', help='input FITS file name')
args = parser.parse_args ()

# FITS file name
file_fits = args.file_input

# check of FITS file name
if not (file_fits[-5:] == '.fits'):
    print ("Input file must be a FITS file")
    sys.exit ()

# function to read header of a FITS file
def read_fits_header (file_fits):
    # opening FITS file
    with astropy.io.fits.open (file_fits) as hdu:
        # reading header and image
        header = hdu[0].header
        # if no image in PrimaryHDU, then read next HDU
        if (header['NAXIS'] == 0):
            header = hdu[1].header
    # returning header and image
    return (header)
  
# reading FITS header
header = read_fits_header (file_fits)

# printing FITS header
print (repr (header))
