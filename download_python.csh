#!/bin/csh

# Time-stamp: <2021/05/21 11:15:55 (CST) daisuke>

# change following line, if you have git command at different location.
set git = '/usr/pkg/bin/git'

# fetching Python scripts
if (-e ${git}) then
    ${git} clone https://github.com/kinoshitadaisuke/ncu_advobs_202102.git
else
    echo "The command \"git\" is not found."
endif
