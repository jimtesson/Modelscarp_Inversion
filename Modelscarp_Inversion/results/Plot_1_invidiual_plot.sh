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
    burnin0=5000; 
    # plot
    function_plot $chain $burnin0 $itmax $hs $agemax

# chain #1
    chain=1;
    burnin1=10000; 
    # plot
    function_plot $chain $burnin1 $itmax $hs $agemax

# erase the merged_results file
: ' rm merged_results.txt
function_merge 0 $burnin0
function_merge 1 $burnin1
function_merge 2 $burnin2
function_merge 3 $burnin3
function_merge 4 $burnin4
function_merge 5 $burnin5
function_merge 6 $burnin6
function_merge 7 $burnin7
function_merge 10 $burnin10
function_merge 11 $burnin11
function_merge 12 $burnin12
function_merge 13 $burnin13
function_merge 14 $burnin14
function_merge 15 $burnin15'

: '
function_grid 0 $burnin
function_grid 1 $burnin
function_grid 2 $burnin
function_grid 3 $burnin
function_grid 4 $burnin
function_grid 5 $burnin
function_grid 6 $burnin
function_grid 7 $burnin
function_grid 8 $burnin
function_grid 9 $burnin
function_grid 10 $burnin
function_grid 11 $burnin
function_grid 12 $burnin
function_grid 13 $burnin
function_grid 14 $burnin
function_grid 15 $burnin



function_grid "merged" 0'

