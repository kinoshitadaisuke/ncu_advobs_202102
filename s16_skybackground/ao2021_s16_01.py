#!/usr/pkg/bin/python3.9

# importing argparse module
import argparse

# importing sys module
import sys

# importing astropy module
import astropy.io.fits

# constructing parser object
desc   = "selecting FITS files"
parser = argparse.ArgumentParser (description=desc)

# adding arguments
parser.add_argument ('-n', '--name', default='', help='target object name')
parser.add_argument ('-e1', '--exptime-min', type=float, default=0.0, \
                     help='minimum exposure time in sec')
parser.add_argument ('-e2', '--exptime-max', type=float, default=3600.0, \
                     help='maximum exposure time in sec')
parser.add_argument ('-a', '--keyword-airmass', default='AIRMASS', \
                     help='FITS keyword for airmass (default: AIRMASS)')
parser.add_argument ('-d', '--keyword-datatype', default='IMAGETYP', \
                     help='FITS keyword for data type (default: IMAGETYP)')
parser.add_argument ('-e', '--keyword-exptime', default='EXPTIME', \
                     help='FITS keyword for exposure time (default: EXPTIME)')
parser.add_argument ('-f', '--keyword-filter', default='FILTER', \
                     help='FITS keyword for filter name (default: FILTER)')
parser.add_argument ('-o', '--keyword-object', default='OBJECT', \
                     help='FITS keyword for object name (default: OBJECT)')
parser.add_argument ('-t', '--keyword-timeobs', default='TIME-OBS', \
                     help='FITS keyword for TIME-OBS (default: TIME-OBS)')
parser.add_argument ('files', nargs='+', help='FITS files')

# parsing arguments
args = parser.parse_args ()

# input parameters
target_name      = args.name
exptime_min      = args.exptime_min
exptime_max      = args.exptime_max
keyword_airmass  = args.keyword_airmass
keyword_datatype = args.keyword_datatype
keyword_exptime  = args.keyword_exptime
keyword_filter   = args.keyword_filter
keyword_object   = args.keyword_object
keyword_timeobs  = args.keyword_timeobs
files_fits       = args.files

# check of target object name
if (target_name == ''):
    print ("Target object name must be specified.")
    sys.exit ()

# making an empty dictionary for storing search results
dic_fits = {}
    
# processing each FITS file
for file_fits in files_fits:
    # check of file name
    if not (file_fits[-5:] == '.fits'):
        print ("The file \"%s\" is not a FITS file!" % file_fits)
        print ("File must be a FITS file. Skipping.")
        continue

    # opening FITS file
    with astropy.io.fits.open (file_fits) as hdu_list:
        # reading header
        header = hdu_list[0].header

    # if target object name does not match with specified target object name,
    # then skip
    if not (header[keyword_object] == target_name):
        continue

    # if exposure time does not match with specified exposure time criteria,
    # then skip
    if not ( (header[keyword_exptime] >= exptime_min) \
             and (header[keyword_exptime] <= exptime_max) ):
        continue

    # adding information to the dictionary
    dic_fits[file_fits] = {}
    dic_fits[file_fits][keyword_airmass]  = header[keyword_airmass]
    dic_fits[file_fits][keyword_datatype] = header[keyword_datatype]
    dic_fits[file_fits][keyword_exptime]  = header[keyword_exptime]
    dic_fits[file_fits][keyword_filter]   = header[keyword_filter]
    dic_fits[file_fits][keyword_object]   = header[keyword_object]
    dic_fits[file_fits][keyword_timeobs]  = header[keyword_timeobs]

# printing results
for file_fits in dic_fits.keys ():
    print (file_fits)
    print ("  %-8s = %-24s    %-8s = %-24s" \
           % (keyword_datatype, dic_fits[file_fits][keyword_datatype], \
              keyword_object, dic_fits[file_fits][keyword_object]) )
    print ("  %-8s = %-24s    %-8s = %f sec" \
           % (keyword_filter, dic_fits[file_fits][keyword_filter], \
              keyword_exptime, dic_fits[file_fits][keyword_exptime]) )
    print ("  %-8s = %-24s    %-8s = %f" \
           % (keyword_timeobs, dic_fits[file_fits][keyword_timeobs], \
              keyword_airmass, dic_fits[file_fits][keyword_airmass]) )
