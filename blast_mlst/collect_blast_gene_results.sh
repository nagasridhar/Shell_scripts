#!/usr/bin/env bash


# This will take in isolate file. Path to all BLAST MLST output/folders.
# Output will be a new file with MLST results concatenated into one single file.
# this has been hard coded as output from MLST BLAST had suffix __scaffolds_  file names.

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
    if [ "$output_dir" == "" ]
    then
	printf "No output directory provided\n"
	usage;
    fi
    
    output_dir=$(cd $output_dir; pwd) #make full path
    if [ ! -d "$output_dir" ]
    then
	print "Incorrect output directory\n"
	usage ;
    fi
} # end of check dir function
#

usage() {   #spit out the usage
cat <<UsageDisplay
#################################################################################

bash collect_blast_gene_results.sh -d <MLST output dir> -i <isolate file> -o <output_dir>

OPTIONS:

-d directory where all BLAST gene (fimH, gyrA) folders are present. 
-i isolate file
-o output directory
Output named full_report_blast_gene_mlst.txt is generated in the output directory 
where all individual folders are created


##################################################################################

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
output_dir=""
while getopts "o:i:d:h" args # iterate over arguments
do
    case $args in
	o)
	    output_dir=$OPTARG;; # path to output dir 
        i)
            isolate_file=$OPTARG;; 
        d) 
            mlst_dir=$OPTARG ;;            # path to mlst outputs
        h)
            usage ;; # help
        *) 
            usage ;; # any other flag which is not an option

    esac # done with case loop
done # completed while loop

## ----  

check_dir;  # check directories
check_files;

#header="subj_assembly   adk     fumC    gyrB    icd     mdh     purA    recA    MLST_type"
header="qseqid  sseqid  pident  length  mismatch        gapopen qstart  qend    sstart  send    evalue  bitscore"
printf "$header\n">>$output_dir/full_report_blast_gene_mlst.txt
for i in `cat $isolate_file`
do
    if [ "$i" != "" ]
    then
	temp_mlst_dir=$mlst_dir/$i
	mlst_file=$i"__scaffolds_blast_hits.txt"	# hard coded #__scaffolds_
 
	if [ -f $temp_mlst_dir/$mlst_file  ]
	then
	    sed -n '1p' $temp_mlst_dir/$mlst_file  | sed 's/__scaffolds//g' >> $output_dir/full_report_blast_gene_mlst.txt
	else
	    printf "$i has no MLST BLAST output present\n"
	fi
	unset temp_mlst_dir mlst_file
    fi # end if not null
done #end of for loop, iterating on isolate file

printf "Please check $output_dir for full_report_blast_gene_mlst.txt\n"
printf "\n<-- TADA -->\n"
