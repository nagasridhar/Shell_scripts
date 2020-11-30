#!/usr/bin/env bash

# This will take in isolate file. Path to all SRST MLST output/folders.
# Output will be a new file with MLST results concatenated into one single file.
# this has been hard coded as output from MLST SRST __mlst__Escherichia_coli#1__results.txt   file names.

check_files(){ 

    if [ "$isolate_file" == "" ]
    then
	printf "\nIsolate file not provided. Exiting...\n"
	usage;
    fi
    if [ ! -f "$isolate_file" ]
    then
	printf "Incorrect isolate file. Exiting...\n"
	usage;
    fi
}
# function ends 

check_dir(){ #check dir function starts here

    if [ "$mlst_dir" == ""  ]
    then
        printf "\nMLST directory is incorrect. Aborting...\n"
        usage;
    fi

    mlst_dir=$(cd $mlst_dir; pwd)
    # get full path .. that is /home/my_user_name/folder_depth1_/blah1/blah_n/file_R1_001.fastq.gz

    if [ ! -d "$mlst_dir" ] 
    then
        printf "\nmlst directory is incorrect. Aborting...\n"
        usage;
    fi
} # end of check dir function
#

usage() {   #spit out the usage
cat <<UsageDisplay
################################################################################

bash collect_SRST_results.sh -d <SRST output dir> -i <isolate file>

OPTIONS:

-d directory where MLST result txt files are present. 
-i isolate file

Output named full_report_srst.txt will be generated in present working directory!

#################################################################################

UsageDisplay
exit;
} # end of usage function --

if [ $# -eq 0 ] # exit if no arguments!
then
    usage;
fi
#

mlst_dir=""
isolate_file=""

while getopts "i:d:h" args # iterate over arguments
do
    case $args in
       
        i)
            isolate_file=$OPTARG;; 
        d) 
            mlst_dir=$OPTARG ;;            # path to mlst outputs
        h) 
            usage ;; # help
        *) 
            usage ;;  # any other flag which is not an option
    esac # done with case loop
done # completed while loop

check_dir;  # check directories
check_files;
#
header="Sample  ST      adk     fumC    gyrB    icd     mdh     purA    recA    mismatches      uncertainty     depth   maxMAF"
#printf "$header\n" >> full_report_srst.txt

for i in `cat $isolate_file`
do
    if [ "$i" != "" ]
    then
	#mlst__Staphylococcus_aureus__results.txt
        #mlst_file=$i"__mlst__Escherichia_coli#1__results.txt" #hard coded #__scaffolds_mlst_calls.txt 
	mlst_file=$i"__mlst__Staphylococcus_aureus__results.txt"
	seq_type=`sed -n '2p' $mlst_dir/$mlst_file` # store seq type line

	if [[ $seq_type ]]
	then
	    sed -n '2p' $mlst_dir/$mlst_file  >> full_report_srst.txt
	else
	    printf "$i has no sequence type information\n"
	fi
        
        unset  mlst_file

    fi # end if not null
done #end of for loop, iterating on isolate file

printf "Please check `pwd` for full_report_srst.txt\n"
printf "You're done using this script. Now use check_blast python script to check your results\n"
printf "\n<-- TADA -->\n"

