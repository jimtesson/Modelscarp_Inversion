#!/bin/bash

# include library of function to plot
. lib_gmt.sh

################################################
########## Main program  #######################
################################################
# merge all desired chains together in "merged_results.txt"

echo ""
echo "Merging chains:"
echo ""

# erase existing output file "the "merged_results.txt"
    rm merged_results.txt

# add chain #0 to merged_results.txt
    chain=0; burnin=5000; 
    function_merge $chain $burnin

# add chain #1 to merged_results.txt
    chain=1; burnin=10000; 
    function_merge $chain $burnin

echo ""
echo " => The merged chains are placed in merged_results.txt"
echo ""
