#!/bin/csh

set url_tar = 'https://www-star.fnal.gov/NorthEqExtension_ugriz/Data/usno40stds.clean.v3.tar.gz'
set url_txt = 'https://www-star.fnal.gov/NorthEqExtension_ugriz/Data/README.txt'

# change following line, if the location of wget is different on your computer.
set wget     = '/usr/pkg/bin/wget'
set wget_opt = '--no-check-certificate'

if (-e ${wget}) then
    ${wget} ${wget_opt} ${url_tar}
    ${wget} ${wget_opt} ${url_txt}
else
    echo "The command \"wget\" is not found."
endif
