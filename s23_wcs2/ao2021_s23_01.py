#!/usr/pkg/bin/python3.9

#
# Time-stamp: <2021/06/10 20:18:25 (CST) daisuke>
#

# importing argparse module
import argparse

# importing sys module
import sys

# importing pathlib module
import pathlib

# importing shutil module
import shutil

# importing astropy module
import astropy.io.fits

# list of data type
list_datatype = ['LIGHT', 'FLAT', 'DARK', 'BIAS']

# constructing parser object
desc = "copying FITS files"
parser = argparse.ArgumentParser (description=desc)

# adding arguments
parser.add_argument ('-t', '--data-type', choices=list_datatype, \
                     default='LIGHT', help='data type (default: LIGHT)')
parser.add_argument ('-e', '--exptime', type=float, default=0.0, \
                     help='exposure time in sec')
parser.add_argument ('-f', '--filter-name', default='', \
                     help='filter name')
parser.add_argument ('-o', '--object-name', default='', \
                     help='object name')
parser.add_argument ('-d', '--dir-dst', default='.', \
                     help='destination directory (default: .)')
parser.add_argument ('files', nargs='+', help='FITS files')

# command-line argument analysis
args = parser.parse_args ()

# input parameters
data_type   = args.data_type
exptime     = args.exptime
filter_name = args.filter_name
object_name = args.object_name
dir_dst     = args.dir_dst
files_fits  = args.files

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

# check of object name
if (object_name == ''):
    print ("Object name has to be specified.")
    sys.exit ()

# check of filter name
if (filter_name == ''):
    print ("Filter name has to be specified.")
    sys.exit ()

# destination directory
path_dst = pathlib.Path (dir_dst)
path_dst.mkdir (mode=0o755, exist_ok=True)
    
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
    if not ( path_fits.exists () and path_fits.is_file () ):
        # printing a message to stderr
        print ("The file \"%s\" does not exist or is not a regular file." \
               % file_fits, file=sys.stderr)
        # skipping
        continue
    
    # reading FITS header
    header = read_fits_header (file_fits)

    # data type
    if not (data_type == header['IMAGETYP']):
        continue

    # exposure time
    if not (exptime == header['EXPTIME']):
        continue

    # filter name
    if ('FILTER' in header):
        if not (filter_name == header['FILTER']):
            continue

    # object name
    if not (object_name == header['OBJECT']):
        continue

    # printing status
    print ("Now copying the file \"%s\" to \"%s\"..." % (file_fits, dir_dst) )
    
    # copying file
    shutil.copy2 (path_fits, path_dst)
