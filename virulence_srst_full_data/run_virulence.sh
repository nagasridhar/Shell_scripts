#!/usr/bin/env bash
#
# Sample = isolate and isolate = sample
# We will trim isolates and then run srst2's virulence tool
#
usage() {   #spit out the usage
    cat<<UsageDisplay
#################################################################################

bash run_virulence.sh -l st_list.txt -d ST_folders_path -o output_dir -f ecoli_db

  
#################################################################################
UsageDisplay
    exit;
} # end of usage function --
##

if [ $# -eq 0 ] # exit if no arguments!
then
    usage;
fi

check_input(){

    if [ "$virulence_fasta" == "" ] || [ "$st_list" == "" ] # isolate_file
    then
	printf "\nVirulence fasta/isolate Files incorrect. Aborting... \n"
	usage;
    fi
    #--------
    
    if [ "$out_dir" == "" ] || [ "$sample_direc" == "" ]
    then
	printf "\nIncorrect directory\n"
	usage;
    fi
    #--------------
    
    virulence_fasta=$(cd $(dirname $virulence_fasta) ; pwd)/`basename $virulence_fasta` # make ful path
    st_list=$(cd $(dirname $st_list) ; pwd)/`basename $st_list` # make ful path
    
    out_dir=$(cd $out_dir;pwd) #make full path
    sample_direc=$(cd $sample_direc;pwd) #make full path

    #---------------
    if [ ! -f "$st_list" ] || [ ! -f "$virulence_fasta"  ]
    then
	printf "Incorrect isolate/virulence fasta file\n"
	usage;
    fi
    #-----------------
    
    if [ ! -d "$out_dir" ] ||  [ ! -d "$sample_direc" ]
    then
	printf "\nSample/output directory incorrect. Aborting...\n"
	usage;
    fi
    #-------------------

}
#----------------function ends----->>>>>>>

sample_direc=""
forward_read=""
reverse_read=""
adap_file=""
out_dir=""
virulence_fasta=""
st_list=""

#-- variables ends------>>>>

adap_file="/home/nagab/scripts/illumina_adapters_all.fasta" # adapter file will be used for trimomatic

while getopts "f:d:l:o:h" args # iterate over arguments
do
    case $args in
	l)st_list=$OPTARG  ;; #store st list to loop over
	
	f)virulence_fasta=$OPTARG;; #store virulence fasta file ;;

	d)sample_direc=$OPTARG ;; #store all ST directory location
	
	h)usage ;;
	
	o)out_dir=$OPTARG ;;
	
	*);;
    esac
done

check_input;

#-- -create isolate file for each ST folder....

for i in `cat $st_list | sort`
do
    if [ "$i" != "" ]
    then
	
	if [ -d $sample_direc/$i ]
	then

	    mkdir $out_dir/$i #St output directory
	    output_st_dir=$out_dir/$i
	    
	    temp_st_dir=$sample_direc/$i #ST symlinks directroy.......
	    
	    array_fastq_file=`ls $temp_st_dir |rev | cut -c 28- | rev | sort | uniq`
	    
	    for x in "${array_fastq_file[@]}"
	    do
		printf "$x" >>$output_st_dir/"isolate.txt"
		# or do whatever with individual element of the array
	    done
	else
	    printf "$i\n Directory is unavailable"
	fi
	
    fi
done

#-----ST folders created.. isolate.txt list create --->>>>>

for i in `cat $st_list | sort`
do
    if [ "$i" != "" ]
    then
	
	#printf "$i\n"	printf $virulence_fasta"\n"
	
	#need symlinks- fastq direcotry too....
	#need ecoli db
	#need isolate.txt file
	#need output direcotry
	
	temp_isolate_file=$out_dir/$i/"isolate.txt"
	temp_output_dir=$out_dir/$i
	temp_fastq_loc=$sample_direc/$i
	#printf $temp_output_dir,$temp_fastq_loc,$temp_isolate_file"\n"
	
	bash ~/scripts/virulence_srst/virulence_srst.sh -f $virulence_fasta -d $temp_fastq_loc -o $temp_output_dir -i $temp_isolate_file

    fi
done
#---------------------------------------------------

#bash ~/scripts/virulence_srst/collect_report.sh -i isolate.txt -o output_nov_17/
