#!/usr/pkg/bin/python3.9

#
# Time-stamp: <2021/06/01 17:13:11 (CST) daisuke>
#

# importing argparse module
import argparse

# importing sys module
import sys

# importing pathlib module
import pathlib

# importing astropy module
import astropy.io.fits

# constructing parser object
desc = "a utility like gethead of WCSTools"
parser = argparse.ArgumentParser (description=desc)

# adding arguments
parser.add_argument ('-k', '--keywords', default='', \
                     help='FITS keywords (e.g. "TIME-OBS,EXPTIME,FILTER")')
parser.add_argument ('files', nargs='+', help='FITS files')

# command-line argument analysis
args = parser.parse_args ()

# input parameters
keywords   = args.keywords
files_fits = args.files

# check of keywords
if (keywords == ''):
    print ("No keyword is given. Keywords must be specified.", file=sys.stderr)
    sys.exit ()

# keywords
list_keywords = keywords.split (',')

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

# processing FITS files one-by-one
for file_fits in files_fits:
    # check of file name
    if not (file_fits[-5:] == '.fits'):
        # printing a message to stderr
        print ("Input file \"%s\" is not a FITS file." % file_fits, \
               file=sys.stderr)
        print ("Input file must be \"*.fits\".")
        # skipping
        continue

    # making pathlib object
    path_fits = pathlib.Path (file_fits)
    filename  = path_fits.name

    # if the file does not exist, then skip
    if not (path_fits.exists ()):
        # printing a message to stderr
        print ("The file \"%s\" does not exist." % file_fits, \
               file=sys.stderr)
        # skipping
        continue
    
    # reading FITS header
    header = read_fits_header (file_fits)

    # making an empty list for storing values of keywords
    list_values = []
    
    # values of keywords
    for keyword in list_keywords:
        # if the keyword exists, then read the value from FITS header
        if (keyword in header):
            value = header[keyword]
        # if the keyword does not exists, then use '__NONE__' as a value
        else:
            value = '__NONE__'
        # appending the value to the list
        list_values.append (value)

    # printing information
    data = str (filename)
    for value in list_values:
        data = data + " " + str (value)
    print (data)
