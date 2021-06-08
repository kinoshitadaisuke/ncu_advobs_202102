#!/usr/pkg/bin/python3.9

#
# Time-stamp: <2021/06/08 19:33:13 (CST) daisuke>
#

# importing argparse module
import argparse

# importing astropy module
import astropy.coordinates
import astropy.units

# importing astroquery module
import astroquery.gaia

# importing ssl module
import ssl

# allow insecure downloading
ssl._create_default_https_context = ssl._create_unverified_context

# constructing parser object
desc = "Downloading Gaia EDR3"
parser = argparse.ArgumentParser (description=desc)

# adding arguments
parser.add_argument ('-c', '--coordinate', default='', \
                     help='RA and Dec (hh:mm:ss.ss,+/-dd:mm:ss.s)')
parser.add_argument ('-o', '--file-output', default='', \
                     help='output file name (*.cat)')
parser.add_argument ('-r', '--radius', type=float, default=0.1, \
                     help='radius in degree (default: 0.1)')

# command-line argument analysis
args = parser.parse_args ()

# input parameters
coordinate  = args.coordinate
file_output = args.file_output
radius_deg  = args.radius

# check of coordinate
if not (',' in coordinate):
    print ("Coordinate must be hh:mm:ss.ss,+/-dd:mm:ss.s format.")
    sys.exit ()

# check of output file
if not ( (file_output[-4:] == '.vot') or (file_output[-7:] == '.vot.gz') ):
    print ("Output file must be VOTable file (*.vot or *.vot.gz)")
    sys.exit ()

# units
u_ha = astropy.units.hourangle
u_deg = astropy.units.deg

# coordinate
(coord_ra, coord_dec) = coordinate.split (',')

# making skycoord object
coord = astropy.coordinates.SkyCoord (coord_ra, coord_dec, \
                                      unit=(u_ha, u_deg), frame='icrs')

# RA and Dec in deg
coord_ra_deg  = coord.ra.deg
coord_dec_deg = coord.dec.deg

# command for database query
query_p = "POINT('ICRS',gaiaedr3.gaia_source.ra,gaiaedr3.gaia_source.dec)"
query_c = "CIRCLE('ICRS',%f,%f,%f)" % (coord_ra_deg, coord_dec_deg, radius_deg)
query   = "SELECT * from gaiaedr3.gaia_source WHERE CONTAINS(%s,%s)=1;" \
    % (query_p, query_c)

# sending a job to Gaia database
job = astroquery.gaia.Gaia.launch_job_async (query, dump_to_file=True, \
                                             output_file=file_output, \
                                             output_format='votable')
print (job)

# getting results
results = job.get_results ()

# printing results
print (results)
