#!/bin/bash

# include library of function to plot
. lib_gmt.sh

################################################
########## Main program  #######################
################################################

gmt gmtset PS_MEDIA A1

# maximum number of iterations
itmax=50000;
# height of the scarp (cm)
hs=1600;
# max age (yr)
agemax=20000;

# plot the density grid
function_grid "merged" $hs $agemax

