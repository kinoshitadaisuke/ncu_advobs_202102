#!/usr/pkg/bin/python3.9

#
# Time-stamp: <2021/06/02 04:29:00 (CST) daisuke>
#

# importing numpy module
import numpy

# importing matplotlib module
import matplotlib.figure
import matplotlib.backends.backend_agg

# data file
file_phot = 'phot_288_289.data'

# output file
file_output = 'limmag_288_289.pdf'

# magnitude of photometric standard star
std_mag = 14.553

# net flux of photometric standard star
std_flux1 = 258876.2009
std_flux2 = 270233.4415

# making empty list
list_mag_mean = []
list_mag_err  = []

# opening file
with open (file_phot, 'r') as fh_in:
    # reading file
    for line in fh_in:
        # skipping line, if the line start with '#'
        if (line[0] == '#'):
            continue
        # splitting line
        (x1_str, y1_str, flux1_str, x2_str, y2_str, flux2_str) = line.split ()
        # conversion from string to float
        x1    = float (x1_str)
        y1    = float (y1_str)
        flux1 = float (flux1_str)
        x2    = float (x2_str)
        y2    = float (y2_str)
        flux2 = float (flux2_str)
        # skipping standard star
        if ( (flux1 == std_flux1) and (flux2 == std_flux2) ):
            continue
        # calculation of magnitude
        mag1 = std_mag - 2.5 * numpy.log10 (flux1 / std_flux1)
        mag2 = std_mag - 2.5 * numpy.log10 (flux2 / std_flux2)
        mag_mean = (mag1 + mag2) / 2.0
        mag_diff = numpy.absolute (mag1 - mag2)
        mag_err  = mag_diff / numpy.sqrt (2.0)
        # appending data to lists
        list_mag_mean.append (mag_mean)
        list_mag_err.append (mag_err)

# making objects "fig" and "ax"
fig = matplotlib.figure.Figure ()
matplotlib.backends.backend_agg.FigureCanvasAgg (fig)
ax = fig.add_subplot (111)

# plotting data
ax.plot (list_mag_mean, list_mag_err, marker='o', color='blue', \
         markersize=5, linestyle='None')
ax.set_xlabel ('Magnitude')
ax.set_ylabel ('Error on Magnitude')
ax.set_title ('60-sec r-band limiting magnitude')

# saving to file
fig.savefig (file_output, dpi=225)
