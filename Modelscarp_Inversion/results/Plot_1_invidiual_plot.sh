#!/bin/bash

# include library of function to plot
. lib_gmt.sh

################################################
########## Main program  #######################
################################################

gmt gmtset PS_MEDIA A1

# maximum number of iterations
itmax=30000;
# height of the scarp (m)
hs=16;
# max age (kyr)
agemax=20;

# chain #0
    chain=0;
    # plot
    function_plot $chain $itmax $hs $agemax

# chain #1
    chain=1;
    # plot
    function_plot $chain $itmax $hs $agemax

