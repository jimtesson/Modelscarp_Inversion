#!/bin/bash

# include library of function to plot
. lib_gmt.sh

################################################
########## Main program  #######################
################################################

gmt gmtset PS_MEDIA A1


# height of the scarp (m)
hs=16;
# max age (kyr)
agemax=20;
# max frequency for histogram of the number of event (%)
max_freq_nev=50
# max frequency for histogram of the number of event (%)
max_freq_SR=20
# max frequency for histogram of the number of event (%)
max_freq_RMSw=20

# chain #0
chain='merged';

# plot
function_plot_histo $chain $hs $agemax $max_freq_nev $max_freq_SR $max_freq_RMSw


