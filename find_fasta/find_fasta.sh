#!/usr/bin/env bash


# This is made for FUTI isolates
# Find scaffolds files 
# and create symlinks

# Take in isolate file name, find them from the input directory path, and create symlink in the input directory

get_scaffolds_files(){
    
    isolate_name="$1";

    array_isolate=($(find $input_dir -name "*$isolate_name*" 2>/dev/null | grep "__scaffolds.fasta" ))

    if [  ${#array_isolate[@]} -eq 0 ]
    then
	printf "$isolate_name has NO file present  \n" >> log.txt
    fi 
    if [ ${#array_isolate[@]} -gt 1 ] #isolate has more than fasta files
    then
	printf "$isolate_name has multiple file  \n" >> log.txt
    fi
    
    for file in "${array_isolate[@]}"
    do
	printf "$isolate_name\t$file\n" >> path.txt ;
    done

}
#function ends

check_inputs(){ # function to check dirctories
    
    if [ "$input_dir" == "" ] || [ "$output_dir" == "" ]
    then
	usage;
    fi
    output_dir=$(cd $output_dir;pwd)
    input_dir=$(cd $input_dir;pwd)

    if [ ! -d  "$input_dir" ] || [ ! -d "$output_dir" ]
    then
	usage;
    fi
}

usage() {   #spit out the usage
cat <<UsageDisplay
###################
bash find_fasta.sh -i <isolate_file.txt> -d <locatoin where isolates are present> -o <directory for output symlinks>

-i isolate list
-d where scaffolds fasta files are present for the isolate input
-o output where symlinks are to be created

###################
UsageDisplay
exit;
}
# usage function ends

if [ $# -eq 0 ] # exit if no arguments!
then
    usage;
fi

# begin variables
isolate_file="" # isolate list
output_dir="" #where you'd like to have output
input_dir="" #where all input data might be present 

while getopts "i:o:d:h" args # iterate over arguments
do
    case $args in
	d)
	    input_dir=$OPTARG ;;
	i)
	    isolate_file=$OPTARG ;;
	o)
	    output_dir=$OPTARG ;;
	h)
	    usage ;; # print usage
    esac # done with case loop
done # completed while loop

check_inputs;
total_isolates=0 #total isolates found in file provided
isolates_found=0  # total isolates for which scaffolds file was found

for i in `cat $isolate_file | sort`
do
#    printf $i"\n"
    total_isolates=`expr 1 + $total_isolates` # increment total isolates --
    get_scaffolds_files $i;
done
