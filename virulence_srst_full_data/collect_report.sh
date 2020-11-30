#!/usr/bin/env bash

# __date__ 20 july 2016
# __author__ Sanjeev Sariya
# __location__ Price Lab 4th Fl, New Hampshire, Cube# 420F

#Loop through output directory, and concatenate all 
#genes report generated from srst virulence

#
# isolate file
# output directory
#

check_input(){
    
    if [ "$isolate_file" == "" ]
    then
	printf "\nisolate Files incorrect. Aborting... \n"
	usage;
    fi
    # -------
    
    if [ "$out_dir" == "" ]
    then
	printf "\nIncorrect directory\n"
	usage;
    fi
    # ------
    
    isolate_file=$(cd $(dirname $isolate_file) ; pwd)/`basename $isolate_file` # make ful path
    out_dir=$(cd $out_dir;pwd)

    if [ ! -f "$isolate_file" ]
    then
	printf "Incorrect isolate fasta file\n"
	usage;
    fi
    #--
    
    if [ ! -d "$out_dir" ]
    then
	printf "\noutput directory incorrect. Aborting...\n"
	usage;
    fi
    #---
}

usage() {   #spit out the usage

cat <<UsageDisplay
#################################################

collect_report.sh -i <isolate file> -o <output direcotry>

-o where indivuldual isolate folders are created
-i isolate file-  contains isolate name - one per line 

#################################################
UsageDisplay
exit;

}
#---------

if [ $# -eq 0 ] # exit if no arguments!
then
    usage;
fi
#-----------

out_dir=""
isolate_file=""

while getopts "i:o:h" args # iterate over arguments
do
    case $args in
	i)
	    isolate_file=$OPTARG ;; # store isolate file
	o)
	    out_dir=$OPTARG ;; # output directory
	h)
	    usage ;;        # help
	*)
	    usage ;;        # any other flag which is not an option
    esac # done with case loop
done # completed while loop

check_input;

for i in `cat $isolate_file | sort`
do
    if [ "$i" != "" ]
    then
	
	temp_file_path="" #concatenate outdir/$i and file path
	result_file=""
	
	#set -x
	result_file=$(ls $out_dir/$i | grep "__fullgenes__ecoli_db__results.txt")
	temp_file_path=$out_dir/$i/$result_file
	
	if [ ! -f "$temp_file_path" ] || [ "$result_file" == "" ]
 	then

	    printf "Issue with isolate \n" >&2 
	    printf $out_dir/$i"\n" >&2 
	    #usage;
	else
	    #cocatenate all __full_gene file
	    #remove out_paired, and other added prefix .. 
	    cat $temp_file_path | sed 's/out_paired_//g' | sed 's/_bcode_lane//g' | sed '1d' >> $out_dir/"final_report.txt"
	fi	

	unset temp_file_path result_file
    fi
done

#add report header .. seprated by Tab
sed -i '1s/^/Sample\tDB\tgene\tallele\tcoveragedepth\tdiffs\tuncertainty\tdivergence\tlength\tmaxMAF\tclusterid\tseqid\tannotation\n/' $out_dir/"final_report.txt"

printf "Check output directory for report\n"
