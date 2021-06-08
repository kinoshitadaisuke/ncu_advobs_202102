#!/usr/pkg/bin/python3.9

#
# Time-stamp: <2021/06/08 18:58:41 (CST) daisuke>
#

# importing sys module
import sys

# importing pathlib module
import pathlib

# importing shutil module
import shutil

# source directory
dir_src = '0214red'

# destination directory
dir_dst = 'v0678vir'

# making a directory if it does not exist
path_dst = pathlib.Path (dir_dst)
path_dst.mkdir (mode=0o755, exist_ok=True)

# reading data from standard input
for line in sys.stdin:
    # if the line starts with '#", then skip
    if (line[0] == '#'):
        continue
    # splitting data
    record = line.split ()
    # data
    filename   = record[0]
    imagetype  = record[1]
    exptime    = float (record[2])
    filtername = record[3]
    objectname = record[4]
    # rejecting data if criteria do not match
    if not ( (imagetype == 'LIGHT') and (exptime == 60.0) \
             and (filtername == 'rp_Astrodon_2019') \
             and (objectname == 'V0678VIR') ):
        continue
    # copying file
    file_fits = "%s/%s" % (dir_src, filename)
    path_fits = pathlib.Path (file_fits)
    print ("copying file from %s to %s..." % (path_fits, dir_dst) )
    shutil.copy2 (path_fits, dir_dst)
