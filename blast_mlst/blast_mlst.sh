#!/usr/bin/env bash


#SBATCH --time=3-0:0:0  # - time to run          #changed on May 17 from 4 to 7 
#SBATCH -o "blast_main_%j.out"  # output 
#SBATCH -J BLA_ST # job name 
#SBATCH  -p 128gb # partition to run 
#SBATCH  --nodes=1 # number of nodes -- 
#SBATCH -e "blast_log_%j.stderr"  # error log name 

# This script assumes BLAST database has been created.
# This will iterate over isolate files, and will create individual folders for each isolate.
# It will pick scaffold files for each isolate from the directory provided. It's hard coded __scaffolds.fasta string
# It assumes blast_MLST.sh is in PATH. Always get rid of clonal complex

module load blast+/2.2.30+ biopython/1.63 python2.7

usage() {   #spit out the usage
cat <<UsageDisplay
#################################################

bash blast_mlst.sh -d < directory of all scaffolds file> -b <directory of BLAST database> -o <output dir> -i <isolate file> -s <organism name>

OPTIONS:

-d directory where all scaffold files are present. 
   this are suffixed with  __scaffolds.fasta after isolate name.

-o output directory.

-b Path/directory where all blast database has been created
Make sure you do not have clonal complex in profile file.

-i isolate file

-s organism name use all small caps.
   Ex: saureus, ecoli

#################################################

UsageDisplay
exit;
} # end of usage function --
##
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

    if [ "$blast_dir" == ""  ] || [ "$scaffolds_dir" == "" ] || [ "$out_dir" == "" ]
    then
	printf "\nScaffolds/blast/output directory is incorrect. Aborting...\n"
	usage;
    fi

    scaffolds_dir=$(cd $scaffolds_dir;pwd)
    blast_dir=$(cd $blast_dir;pwd)
    out_dir=$(cd $out_dir; pwd)

    # get full path .. that is /home/my_user_name/folder_depth1_/blah1/blah_n/file_R1_001.fastq.gz

    if [ ! -d "$blast_dir" ] || [ ! -d "$scaffolds_dir" ] || [ ! -d "$out_dir" ]
    then
        printf "\nScaffolds/blast/output directory is incorrect. Aborting...\n"
        usage;
    fi
} # end of check dir function
#

if [ $# -eq 0 ] # exit if no arguments!
then
    usage;
fi
#

blast_dir=""
scaffolds_dir=""
out_dir=""
isolate_file=""
organism=""
##

while getopts "s:i:b:d:o:h" args # iterate over arguments
do
    case $args in
	s)
	    organism=$OPTARG ;; #store organism name
	i)
	    isolate_file=$OPTARG;; #store isolate file
        d) 
	    scaffolds_dir=$OPTARG ;;	    # path to scaffolds dir	
	b)
	    blast_dir=$OPTARG ;;	    # path to blast dir
        o)
            out_dir=$OPTARG ;; # output directory
        h) 
            usage ;;             # help
        *) 
            usage ;;             # any other flag which is not an option

    esac # done with case loop
done # completed while loop

## ----  

check_dir;  # check directories
check_files;

if [ "$organism" == "" ]
then
    printf "\nProvide organism name in small caps with -s flag\n"
    printf "Ex: saureus, ecoli...\n"
    printf "Aborting...\n"
    exit;
fi
# <<< ---- Loop over scaffolds directory provided  - -- >>> #

for i in `cat $isolate_file`
do
    if [ "$i" != "" ]
    then
	mkdir $out_dir/$i
	scaf_file=$i"__scaffolds.fasta"
	pushd $out_dir/$i > /dev/null

	bash blast_MLST.sh $blast_dir $scaffolds_dir/$scaf_file $organism 100
	popd > /dev/null
    fi
done #end of for loop, iterating on isolate file

printf "Done BLAST-ing MLST\n"
# <<< --loop over mlst outputs -->>>

for i in `cat $isolate_file`
do
    if [ "$i" != "" ]
    then
	temp_mlst_dir=$out_dir/$i
	mlst_file=$i"__scaffolds_mlst_calls.txt"        # hard coded #__scaffolds_mlst_calls.txt

	if [ -f $temp_mlst_dir/$mlst_file  ] #if file exists
	then
	    sed -n '2p ' $temp_mlst_dir/$mlst_file | sed 's/__scaffolds//g' >> $out_dir/full_report_blast_mlst.txt
	else
	    printf "$i has no MLST BLAST output present\n"
	fi
	unset temp_mlst_dir mlst_file
    fi # end if not null
done #end of for loop, iterating on isolate file

printf "Done collecting reports\n" 
printf "<-- TADA -->\n"
