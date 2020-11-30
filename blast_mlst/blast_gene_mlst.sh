#!/usr/bin/env bash


# This script assumes BLAST database has been created.
# This will iterate over isolate files, and will create individual folders for each isolate.
# It will pick scaffold files for each isolate from the directory provided. It's hard coded __scaffolds.fasta string
# It assumes blast_gene.sh is in PATH. 

module load blast+/2.2.30+ biopython/1.63 python2.7

usage() {   #spit out the usage
cat <<UsageDisplay
#################################################

bash blast_gene_mlst.sh -d < directory of all scaffolds file> -b <directory of BLAST database> -o <output dir> -i <isolate file> -p <percent identity>

OPTIONS:

-p percent identity in BLAST

-d directory where all scaffold files are present. 
   this are suffixed with  __scaffolds.fasta after isolate name.

-o output directory.

-b Path/directory where all blast database has been created
Make sure you do not have clonal complex in profile file.

-i isolate file

#################################################

UsageDisplay
exit;
} # end of usage function --
##

check_input(){
    
    #check_directory
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
    #check Files
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
    if [ "$percent" == "" ]
    then
	printf "\nProvide percent identity\n"
	usage ;
    fi
} #fucntion ends 
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
percent=""
##

while getopts "p:i:b:d:o:h" args # iterate over arguments
do
    case $args in
	p)
	    percent=$OPTARG ;;
	i)
	    isolate_file=$OPTARG;; 
        d) 
	    scaffolds_dir=$OPTARG ;; 	    # path to scaffolds dir	
	b)
	    blast_dir=$OPTARG ;; 	    # path to blast dir
        o)
            out_dir=$OPTARG ;; # output directory
        h) 
            # help
            usage ;;
        *) 
            # any other flag which is not an option
            usage ;;

    esac # done with case loop
done # completed while loop

## ----  

check_input;

# -- ----           -----          ---            ----         Loop over scaffolds directory provided  --- -- -- -- -- -- -- -- --    #

for i in `cat $isolate_file`
do
    if [ "$i" != "" ]
    then
	mkdir $out_dir/$i
	scaf_file=$i"__scaffolds.fasta"
	pushd $out_dir/$i > /dev/null
	
	bash blast_gene.sh $blast_dir $scaffolds_dir/$scaf_file $percent
	
	popd > /dev/null
    fi
done #end of for loop, iterating on isolate file

printf "We are done blast_gene_mlst.sh script\n"
printf "Please go ahead and use collect BLAST gene mlst script\n\n"

bash /home/nagab/scripts/blast_mlst/collect_blast_gene_results.sh -d $out_dir -i $isolate_file -o $out_dir

printf "<-- TADA -->\n"
