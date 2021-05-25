#!/usr/pkg/bin/python3.9

#
# Time-stamp: <2021/05/25 23:54:07 (CST) daisuke>
#

# importing argparse module
import argparse

# importing sys module
import sys

# importing astropy module
import astropy.table

# constructing parser object
desc   = 'reading table from a file'
parser = argparse.ArgumentParser (description=desc)

# adding arguments
parser.add_argument ('-c', '--catalogue-file', default='', \
                     help='input catalogue file name')

# command-line argument analysis
args = parser.parse_args ()

# catalogue file name
file_catalogue = args.catalogue_file

# check of catalogue file name
if (file_catalogue == ''):
    print ("Catalogue file name must be specified.")
    sys.exit ()

# reading catalogue from a file
table_source = astropy.table.Table.read (file_catalogue, \
                                         format='ascii.commented_header')

# printing table
print (table_source)
