#!/usr/pkg/bin/python3.9

#
# Time-stamp: <2021/06/08 17:30:07 (CST) daisuke>
#

# importing astropy module
import astropy.coordinates
import astropy.units

# importing astroquery module
import astroquery.gaia

# importing ssl module
import ssl

# allow insecure downloading
ssl._create_default_https_context = ssl._create_unverified_context

# units
u_ha = astropy.units.hourangle
u_deg = astropy.units.deg

# coordinate
coord_ra  = "18:36:56.3"
coord_dec = "+38:47:01"

# search radius in deg
radius_deg = 0.1

# output file name
file_output = 'sample_gaia_edr3.vot.gz'

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
