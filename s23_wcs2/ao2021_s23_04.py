#!/usr/pkg/bin/python3.9

#
# Time-stamp: <2021/06/10 22:47:15 (CST) daisuke>
#

# importing numpy module
import numpy
import numpy.ma

# importing astropy module
import astropy.coordinates
import astropy.stats
import astropy.units

# data file
file_data = 'v0678vir_astrometry.data'

# units
u_ha  = astropy.units.hourangle
u_deg = astropy.units.deg

# making empty Numpy arrays
array_ra = numpy.array ([])
array_dec = numpy.array ([])

# opening data file
with open (file_data, 'r') as fh:
    # reading data file
    for line in fh:
        # splitting data
        records = line.split ()
        # RA in deg
        ra_deg = float (records[1])
        # Dec in deg
        dec_deg = float (records[2])
        # appending data to Numpy arrays
        array_ra  = numpy.append (array_ra, ra_deg)
        array_dec = numpy.append (array_dec, dec_deg)

# sigma-clipping
array_ra_masked  = astropy.stats.sigma_clip (array_ra, sigma=3.0, \
                                             maxiters=10, cenfunc='median', \
                                             masked=True)
array_dec_masked = astropy.stats.sigma_clip (array_dec, sigma=3.0, \
                                             maxiters=10, cenfunc='median', \
                                             masked=True)

# printing results of sigma-clipping
print (array_ra_masked)
print (array_dec_masked)

# mean and standard deviation
ra_mean    = numpy.ma.mean (array_ra_masked)
ra_stddev  = numpy.ma.std (array_ra_masked)
dec_mean   = numpy.ma.mean (array_dec_masked)
dec_stddev = numpy.ma.std (array_dec_masked)

# coordinate
coord = astropy.coordinates.SkyCoord (ra_mean, dec_mean, \
                                      unit=(u_deg, u_deg), frame='icrs')
coord_str = coord.to_string ('hmsdms')
(coord_ra_str, coord_dec_str) = coord_str.split ()

# printing results
print ("RA  = %12.8f deg = %s" % (ra_mean, coord_ra_str) )
print ("Dec = %12.8f deg = %s" % (dec_mean, coord_dec_str) )
print ("stddev of RA  = %6.3f arcsec" % (ra_stddev * 3600.0) )
print ("stddev of Dec = %6.3f arcsec" % (dec_stddev * 3600.0) )
