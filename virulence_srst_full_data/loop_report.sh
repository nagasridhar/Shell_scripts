#!/usr/bin/env bash

#Loop through output directory, and concatenate all
#genes report generated from srst virulence


check_input(){
    
    if [ "$st_list" == "" ]
    then
	printf "\nisolate Files incorrect. Aborting... \n"
	usage;
    fi

    if [ "$st_out_dir" == "" ]
    then
	printf "\nIncorrect directory\n"
	usage;
    fi
    st_list=$(cd $(dirname $st_list) ; pwd)/`basename $st_list` # make ful path
    st_out_dir=$(cd $st_out_dir;pwd)

    if [ ! -f "$st_list" ]
    then
	printf "Incorrect isolate fasta file\n"
	usage;
    fi
    #--

    if [ ! -d "$st_out_dir" ]
    then
	printf "\noutput directory incorrect. Aborting...\n"
	usage;
    fi
    
}
#----------------------------------------------------------------->>

usage(){   #spit out the usage

cat <<UsageDisplay
##############################################################
loop_report.sh -s <stlist.txt> -d <output directory>

-s list with ST names

-d output directory where all ST folders are present
##############################################################
UsageDisplay
exit;
}

if [ $# -eq 0 ] # exit if no arguments!
then
    usage;
fi

#----------

st_out_dir=""
st_list=""

while getopts "s:d:h" args # iterate over arguments
do
    case $args in
	s)
	    st_list=$OPTARG ;; # store isolate file
	d)
	    st_out_dir=$OPTARG ;; # output directory
	h)
	    usage ;;        # help
	*)
	    usage ;;        # any other flag which is not an option
    esac # done with case loop
done # completed

check_input;
for x in `cat $st_list | sort`
do
    if [ "$x" != "" ]
    then
	
	#printf "$x\n"
	#printf $st_out_dir"\n"
	
	if [ -d $st_out_dir/$x ]
	then
	    temp_dir=$st_out_dir/$x
	    #printf $temp_dir"\n"
	    
	    if [ -f $temp_dir/"isolate.txt" ]
	    then
		isolate_file=$temp_dir/"isolate.txt"

		cat $isolate_file | xargs -I % echo % >> $temp_dir/"new_isolate.txt"
		new_isolate_file=$temp_dir/"new_isolate.txt"

		bash /home/nagab/scripts/virulence_srst_full_data/collect_report.sh -i $new_isolate_file -o $temp_dir
	    else
		printf "Error with isolate file $x \n"
	    fi
	else
	    printf "$x folder is missing\n"
	fi
    fi
    #----if not null
done
