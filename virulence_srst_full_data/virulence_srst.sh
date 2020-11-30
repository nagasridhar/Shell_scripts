#!/usr/bin/env bash

#
# Sample = isolate and isolate = sample
# We will trim isolates and then run srst2's virulence tool
#

get_fastq_files(){ # get reverse and foward strand of sample isoltae

    # paramter is $1
    trim_forward=""
    trim_reverse="" # variables to store trim reads, if available

    regex_for='(\_L00[0-9]+)?\_R1+(\_00[1-9])?(\_trimmed)?\.fastq\.gz$' # updated on Apr 24 to have ST131 data. which do not have L001 format
    #regex_for='\_L00[0-9]+\_R1+(\_001)?\.fastq\.gz$'

    regex_rev='(\_L00[0-9]+)?\_R2+(\_00[1-9])?(\_trimmed)?\.fastq\.gz$'
    #regex_rev='\_L00[0-9]+\_R2+(\_001)?\.fastq\.gz$'

    #updated find on 24th April to search symlinks.
    array_name=($(find -L $sample_direc -name "*$1_*" 2>/dev/null)) # store find files in an array --

    if [ ${#array_name[@]} -gt 2 ]
    then
	printf "$1 has more than 2 files - Warning\n"  >> "logs_errors.txt"
    fi  # if there're trimmed or some other garbage.. 

    if [ ${#array_name[@]} -gt 0 ]  # if greater than 0 size
    then

	for file in "${array_name[@]}" # iterate over array
	do
	    if [[ ${file: -9} == ".fastq.gz"  ]]  # if fastq.gz extension
	    then
		if [[ "`basename $file`" =~ $regex_for ]] # if R1 or not
		then
		    if [[ "`basename $file`" == *"trimmed"*  ]] # if R1 and trimmed
		    then
			trim_forward=$file # store in trim_forward variable .. -

		    else
			forward_read=$file
		    fi
		elif [[ "`basename $file`" =~ $regex_rev ]]  #else # if it has R2 -- this got updated.. because it was failing on isolate names..
		then

		    if [[ "`basename $file`" == *"trimmed"*  ]]
		    then
			trim_reverse=$file
		    else
			reverse_read=$file # store in reverse+non_trim variable
		    fi
		fi # check for forward and reverse reads
	    else
		if [[ ${file: -4} != ".txt"  ]] # if file doesn't end with .txt & doesn't have .fastq.gz extension --
		then
		    printf "$file\n" >> "logs_errors.txt"
		    printf "$1 doesn't have .fastq.gz files in its folder\n"  >> "logs_errors.txt"

		fi # end of .txt file check if
	    fi # if file has .fastq.gz extension
	done # end of file stored in array
    else
	printf "File not found for $i\n" >> "logs_errors.txt"
    fi
    if [[ "$trim_forward" != "" ]] && [[ "$trim_reverse" != "" ]] # if there are trimmed reads available !!
    then

	forward_read=$trim_forward  # assign trim read file to normal forward reads file
	reverse_read=$trim_reverse
    fi

        unset trim_forward trim_reverse # unset them to NULL

}
#function ends

check_input(){
    # -- function begins    

    if [ "$virulence_fasta" == "" ] || [ "$isolate_file" == "" ] # isolate_file
    then
	printf "\nVirulence fasta/isolate Files incorrect. Aborting... \n"
	usage;
    fi

    if [ "$out_dir" == "" ] || [ "$sample_direc" == "" ]
    then
	printf "\nIncorrect directory\n"
	usage;
    fi
    #----

    virulence_fasta=$(cd $(dirname $virulence_fasta) ; pwd)/`basename $virulence_fasta` # make ful path
    isolate_file=$(cd $(dirname $isolate_file) ; pwd)/`basename $isolate_file` # make ful path
    
    out_dir=$(cd $out_dir;pwd)
    sample_direc=$(cd $sample_direc;pwd)

    if [ ! -f "$isolate_file" ] || [ ! -f "$virulence_fasta"  ]
    then
	printf "Incorrect isolate/virulence fasta file\n"
	usage;
    fi
    
    if [ ! -d "$out_dir" ] ||  [ ! -d "$sample_direc" ]
    then
	printf "\nSample/output directory incorrect. Aborting...\n"
	usage;
    fi

}
#-- function ends

usage() {   #spit out the usage
cat <<UsageDisplay
#################################################

virulence_srst.sh -i <isolate_file> -o <output_Directory> -d <isolate directory> -f <virulence fasta file> 

isolate file - file containing isolate names
output directory - output files to ge written
virulence fasta - fasta after clustering using cdhit
isolate directory - where all symlinked fastq files are present 

#################################################

UsageDisplay
exit;
} # end of usage function --
##

sample_direc=""
forward_read=""
reverse_read=""
adap_file=""
out_dir=""
virulence_fasta=""
isolate_file=""

adap_file="/home/nagab/scripts/illumina_adapters_all.fasta" # adapter file will be used for trimomatic

if [ $# -eq 0 ] # exit if no arguments!
then
    usage;
fi

while getopts "f:d:i:o:h" args # iterate over arguments
do
    case $args in
	f)
	    virulence_fasta=$OPTARG;; #store virulence fasta file
	d)
	    sample_direc=$OPTARG ;; # store sample direc
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


total_isolates=0 #total isolates found in file provided
isolates_found=0  # total isolates which have both forward and reverse reads  

for i in `cat $isolate_file | sort`
do
    
    if [ "$i" != "" ]
    then

	total_isolates=`expr 1 + $total_isolates` # increment total isolates --

	get_fastq_files $i; # function to get forward and reverse read

	if [ "$forward_read" != "" ] && [ "$reverse_read" != "" ]
	then

	    isolates_found=`expr 1 + $isolates_found` # increment isolates found
	    mkdir $out_dir/$i # make direc of isolate name in output direc.

	    forward_read=$(cd $(dirname $forward_read); pwd)/`basename $forward_read` #/`basename "$forward_read"`
	    reverse_read=$(cd $(dirname $reverse_read); pwd)/`basename $reverse_read`

	    ln -s $forward_read $out_dir/$i/`basename "$forward_read"` # make symbolic link in the directory created for forward and reverse read --
	    ln -s $reverse_read $out_dir/$i/`basename "$reverse_read"`

	    forward_link=$out_dir/$i/`basename $forward_read` # save symbolik links
	    reverse_link=$out_dir/$i/`basename $reverse_read`

	    job_id=$(sbatch /home/nagab/scripts/virulence_srst/run_srst.sh -f $forward_link -r $reverse_link -a $adap_file -o $out_dir/$i -v $virulence_fasta | awk '{print $4}')
	    printf "$i submitted for assembly to $job_id\n"  >> "virulence_job_ids.txt"
	    
	    unset forward_link reverse_link job_id
	    
	fi
	#end of not nul if check for forward and reverse reads
	
	unset forward_read reverse_read # make forward and reverse reads back to null -- 
	#
    fi    
done
# for loop of cat file ends

printf "Total Isolates are $total_isolates\n" >> "logs_errors.txt"
printf "Isolates with forward and reverse reads are $isolates_found\n"  >> "logs_errors.txt"

printf "Check log and virulence job text files in working directory\n"
