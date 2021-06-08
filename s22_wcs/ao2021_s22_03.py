#!/usr/pkg/bin/python3.9

#
# Time-stamp: <2021/06/08 18:38:14 (CST) daisuke>
#

# importing argparse module
import argparse

# importing sys module
import sys

# importing numpy module
import numpy

# importing astropy module
import astropy.io.votable

# command-line argument analysis
parser = argparse.ArgumentParser (description='reading VOTable file')
parser.add_argument ('-i', '--input', help='input VOTable file name')
parser.add_argument ('-o', '--output', help='output catalogue file name')
args = parser.parse_args ()

# VOTable file name
file_votable   = args.input
file_catalogue = args.output

# check of VOTable file name
if not ( (file_votable[-4:] == '.vot') or (file_votable[-7:] == '.vot.gz') ):
    print ("Input file must be VOTable file (*.vot or *.vot.gz)")
    sys.exit ()

# check of catalogue file name
if not (file_catalogue[-4:] == '.cat'):
    print ("Output file must be catalogue file (*.cat)")
    sys.exit ()

# reading VOTable
table = astropy.io.votable.parse_single_table (file_votable).to_table ()

# data
data_id        = numpy.array (table['source_id'])
data_desig     = numpy.array (table['designation'])
data_ra        = numpy.array (table['ra'])
data_dec       = numpy.array (table['dec'])
data_parallax  = numpy.array (table['parallax'])
data_pmra      = numpy.array (table['pmra'])
data_pmdec     = numpy.array (table['pmdec'])
data_b         = numpy.array (table['phot_bp_mean_mag'])
data_g         = numpy.array (table['phot_g_mean_mag'])
data_r         = numpy.array (table['phot_rp_mean_mag'])
data_br        = numpy.array (table['bp_rp'])
data_bg        = numpy.array (table['bp_g'])
data_gr        = numpy.array (table['g_rp'])
data_ra_err    = numpy.array (table['ra_error'])
data_dec_err   = numpy.array (table['dec_error'])
data_pmra_err  = numpy.array (table['pmra_error'])
data_pmdec_err = numpy.array (table['pmdec_error'])
data_p_snr     = numpy.array (table['parallax_over_error'])
data_b_snr     = numpy.array (table['phot_bp_mean_flux_over_error'])
data_g_snr     = numpy.array (table['phot_g_mean_flux_over_error'])
data_r_snr     = numpy.array (table['phot_rp_mean_flux_over_error'])

# writing data to file
with open (file_catalogue, 'w') as fh:
    # printing header
    fh.write ("gaia_edr3/d/j\n")
    fh.write ("a tiny portion of Gaia EDR3\n")
    # writing each object
    for i in range ( len (data_id) ):
        if ( numpy.isnan (data_r[i]) ):
            fh.write ("%012d %12.8f %+13.8f %7.3f\n" \
                      % (i, data_ra[i], data_dec[i], 99.999) )
        else:
            fh.write ("%012d %12.8f %+13.8f %7.3f\n" \
                      % (i, data_ra[i], data_dec[i], data_r[i]) )
