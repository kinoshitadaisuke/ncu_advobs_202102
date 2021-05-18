#!/usr/pkg/bin/python3.9

# importing argparse module
import argparse

# importing astropy module
import astropy.coordinates

# constructing parser object
desc   = "showing coordinates of photometric standard stars"
parser = argparse.ArgumentParser (description=desc)

# adding arguments
parser.add_argument ('files', nargs='+', help='photometric standard star files')

# parsing arguments
args = parser.parse_args ()

# input parameters
files_cat = args.files

# printing header
print ("# ID, RA (deg), Dec (deg), RA (hex), Dec (hex), g' mag, r' mag, i' mag")

# processing each catalogue file
for file_cat in files_cat:
    # opening the file
    with open (file_cat, 'r') as fh:
        # reading the file line-by-line
        for line in fh:
            # if the line starts with '#', then skip
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
            # g'-band magnitude
            mag_g   = float (records[6])
            # r'-band magnitude
            mag_r   = float (records[9])
            # i'-band magnitude
            mag_i   = float (records[12])
            # coordinates
            coord = astropy.coordinates.SkyCoord (ra_deg, dec_deg, unit='deg')
            # conversion into hmsdms format
            radec_str = coord.to_string ('hmsdms')
            (ra_str, dec_str) = radec_str.split ()
            # printing coordinate
            print ("%03d %9.5f %+8.4f %-16s %-16s %6.3f %6.3f %6.3f" \
                   % (star_id, ra_deg, dec_deg, ra_str, dec_str, \
                      mag_g, mag_r, mag_i) )
