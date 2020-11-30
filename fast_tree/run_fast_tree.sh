#!/usr/bin/env bash

# run FAST-tree
# This script takes in input directory .. Input directory has ST folders..
# Each ST folder has output folder generated from gubbins
# each of these output folder has tree.filtered_polymorphic_sites.fasta edited file
# Edited: Reference has been removed

module load fasttree/2.1.8

usage() {   #spit out the usage
    cat <<UsageDisplay
#################################################
run_fast_tree -d <input dir> -o <output>

#################################################
UsageDisplay
exit;
}

check_input(){
    if [ "$input_dir" == "" ] || [ ! -d "$input_dir" ]
    then
	printf "Incorrect input directory\n"
	usage;
    fi

    input_dir=$(cd $(dirname $input_dir) ; pwd)/`basename $input_dir` # make full path ..

    if [ "$output_dir" == "" ] || [ ! -d "$output_dir" ]
    then
	printf "Incorrect output directory\n"
	usage;
    fi
    output_dir=$(cd $(dirname $output_dir) ; pwd)/`basename $output_dir` # make full path ..
}

if [ $# -eq 0 ] # exit if no arguments!
then
    usage;
fi

input_dir=""
output_dir=""

while getopts "o:d:h" args # iterate over arguments
do
    case $args in
	d)
	    input_dir=$OPTARG;; #store input direct
	o)
	    output_dir=$OPTARG;;
	h)
	    usage;

    esac # done with case loop
done # completed while loop

check_input;
for i in `ls $input_dir | sort`
do
    dir=$input_dir/$i #temp dir
    input_file=$dir/"tree.filtered_polymorphic_sites.fasta" #FASTA file to run FAST-tree on
    
    if [ -f $input_file ]
    then
	
	mkdir $output_dir/$i
	FastTreeMP -gtr -gamma -log $output_dir/$i/"fast_tree.log" -nt $input_file > $output_dir/$i/"out.tre" 2> $output_dir/$i/"run_error.log"
    else
	printf "$i has some issues\n"
    fi
    
    
done
printf "Done looping! Script run complete\n"
printf "Bye bye\n\n"
