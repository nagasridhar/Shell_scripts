#!/usr/bin/env bash

# __date__ 19 July 2016
# __author__ Sanjeev Sariya
# __location__ Price Lab 4th Fl, New Hampshire, Cube# 420F

#SBATCH --time=4-0:0:0  # - time to run          #changed on May 17 from 4 to 7
#SBATCH -o "virulence_srst.%j.out"  # output
#SBATCH -J srst_virlnce # job name
#SBATCH  -p 128gb # partition to run
#SBATCH  --nodes=1 # number of nodes --
#SBATCH -e "log_virulence.%j.stderr"  # error log name

# This script is called in from virulence_srst.sh
# All sanity checks have already been performed in there.
# We here do: virulence typing using SRST (trimmomatic reads) 
# trim reads, and then use srst for virulence detection 

# ------       ------          ------                          ------           ------     ------   ------             ------

module load python2.7 samtools/0.1.18 bowtie2/2.2.3
usage() {   #spit out the usage
cat <<UsageDisplay
#################################################

run_srst.sh  -f <forward_read> -r <reverse_read> -o <output_Directory>  -a <illumina_adapter_file>

#################################################
UsageDisplay
exit;
} # end of usage function 

if [ $# -eq 0 ] # exit if no arguments!
then
    usage;
fi

perform_srst(){ # function begins for doing SRST Sequence typing ..
    
    # $1 is dir of trimmed reads # $2 forward trimmed reads  # $3 reverse trimmed reads

    temp_for=$1/$2 #concatenate out_dir and forward read
    temp_rev=$1/$3  # concatenate out_dir and reverse read

    index_read=$(echo $temp_for | grep -bo "_R1" | awk 'BEGIN {FS=":"}{print $1}'  )
    
    if [ "$index_read" == "" ] # _R1 pattern isn't avail..  Exit function
    then
	printf "\nThe reads here do not match our pattern for $isolate. Exiting Sequence typing function. Team Fall back!\n "  >> "logs_srst.txt"
	return 1;
    fi
    
    index_fastq=$(echo $temp_for | grep -bo ".fastq.gz" | awk 'BEGIN {FS=":"} {print $1}')
    for_suff=${temp_for:$index_read:$index_fastq-$index_read} # get forward suffix

    index_read=$(echo $temp_rev | grep -bo "_R2" | awk 'BEGIN {FS=":"}{print $1}' )

    if [ "$index_read" == "" ] # _R2 pattern isn't available.. Exit function
    then
	printf "\nThe reads here do not match our pattern $isolate. Exiting Sequence typing function. Team Fall back! \n "  >> "logs_srst.txt"
	return 1;
    fi

    index_fastq=$(echo $temp_rev | grep -bo ".fastq.gz" | awk 'BEGIN {FS=":"} {print $1}')
    rev_suff=${temp_rev:$index_read:$index_fastq-$index_read} # get reverse suffix

    set -x
    python ~/new_srst2/srst2-master/scripts/srst2.py --input_pe $temp_for $temp_rev --forward $for_suff --reverse $rev_suff --output $out_dir/$isolate  --log --gene_db $virulence_fasta 
    set +x
    
}
# -- function ends

out_dir=""
forward_read=""
reverse_read=""
adap_file=""
out_dir=""
isolate=""
virulence_fasta=""

while getopts "v:r:f:a:i:o:h" args # iterate over arguments
do
    case $args in
	v)
	    virulence_fasta=$OPTARG;;
	f)
	    forward_read=$OPTARG;;
	r)
	    reverse_read=$OPTARG;;
	a)
	    adap_file=$OPTARG;;
	o)
	    out_dir=$OPTARG ;; # output directory
	h)
	    usage ;;        # help
	*)
	    usage ;;        # any other flag which is not an option

    esac # done with case loop
done # completed while loop

temp_forw_read=`basename $forward_read` ; # get only base name of input file provided
temp_rev_read=`basename $reverse_read` ; # --||--
isolate=`basename $out_dir`

printf "we are going to run trimmomatic with phred33 settings\n"

java -jar /c1/apps/trimmomatic/Trimmomatic-0.33/trimmomatic-0.33.jar PE -phred33 -threads 16 -trimlog $out_dir/"trim_log.txt" $forward_read $reverse_read $out_dir/"out_paired_"$temp_forw_read $out_dir/"out_unpaired_"$temp_forw_read $out_dir/"out_paired_"$temp_rev_read $out_dir/"out_unpaired_"$temp_rev_read ILLUMINACLIP:$adap_file:2:25:10 MINLEN:30

rm $out_dir/"out_unpaired_"$temp_forw_read $out_dir/"out_unpaired_"$temp_rev_read

perform_srst $out_dir "out_paired_"$temp_forw_read "out_paired_"$temp_rev_read;


