#!/usr/bin/env bash


#This calls in  blast_gene_mlst.sh  script
# It will create sepaate folder for 100, and 80 percent
#In each 100 and 80 percent 5 folders be created for outputs
#
#  bash automate_gene_blast.sh -l ./ -p test.txt -z /groups/pricelab/FUTI_analysis/FUTI_ecoli_assembly_symlinked/ST10
check_inputs(){

    if [ "$gene_output_location" == "" ] || [ "$isolate_file" == "" ] || [ "$scaffolds_dir" == "" ]
    then
	usage;
    fi
     gene_output_location=$(cd $gene_output_location; pwd)
     scaffolds_dir=$(cd $scaffolds_dir; pwd)
     isolate_file=$(cd $(dirname $isolate_file) ; pwd)/`basename $isolate_file` # make full path .. 

    if [ ! -d  "$gene_output_location" ] || [ ! -f "$isolate_file"  ] || [ ! -d "$scaffolds_dir"  ]
    then
	usage;
    fi
}
#ends
usage() {   #spit out the usage
cat <<UsageDisplay
#################################################################################

bash automate_gene_blast.sh -z <all fasta symlinked> -p <isolate file> - l <output directory>

OPTIONS:

-z directory where all fasta are present

-p isolate file

-l output directory where fimH, gyrA, other will be created


##################################################################################

UsageDisplay
exit;
} # end of usage function 
#

if [ $# -eq 0 ] # exit if no arguments!
then
    usage;
fi

gene_output_location=""
isolate_file=""
scaffolds_dir=""

while getopts "z:p:l:h" args # iterate over arguments
do
    case $args in
       
        z)
            scaffolds_dir=$OPTARG;; 
        p) 
	    isolate_file=$OPTARG;; 
        l)
	    gene_output_location=$OPTARG;;
	h)
            usage ;; # help
        *) 
            usage ;; # any other flag which is not an option

    esac # done with case loop
done # completed while loop

check_inputs;

pushd $gene_output_location > /dev/null
{
    
    mkdir $gene_output_location/80_percent
    pushd $gene_output_location/80_percent > /dev/null
    {
	mkdir parC
	pushd  $gene_output_location/80_percent/parC > /dev/null
	{

	    mkdir $gene_output_location/80_percent/parC/blast_results  
	    pushd $gene_output_location/80_percent/parC/blast_results > /dev/null
	    {
		bash /home/nagab/scripts/blast_mlst/blast_gene_mlst.sh -d $scaffolds_dir -b /lustre/groups/pricelab/_files/st_gene_mlst/parC/db -o $gene_output_location/80_percent/parC/blast_results -i $isolate_file -p 80
	    }
	    popd > /dev/null #exit blast_results
	    bash /home/nagab/scripts/blast_mlst/collect_blast_gene_results.sh -d $gene_output_location/80_percent/parC/blast_results -i $isolate_file
	}
	popd >/dev/null #exit parC 
    }
    popd > /dev/null #exit 80_percent
}
popd  > /dev/null #exit destination


printf "BLAST Done 80 threshold  parC\n"

pushd $gene_output_location > /dev/null
{
    pushd $gene_output_location/80_percent > /dev/null
    {
	mkdir fimH
	pushd  $gene_output_location/80_percent/fimH> /dev/null
	{
	    mkdir blast_results
	    pushd $gene_output_location/80_percent/fimH/blast_results >/dev/null
	    {
		
		bash /home/nagab/scripts/blast_mlst/blast_gene_mlst.sh -d $scaffolds_dir -b /lustre/groups/pricelab/_files/st_gene_mlst/fimH/db -o $gene_output_location/80_percent/fimH/blast_results -i $isolate_file -p 80
	    }
	    popd > /dev/null #exit blast_results
	    bash /home/nagab/scripts/blast_mlst/collect_blast_gene_results.sh -d $gene_output_location/80_percent/fimH/blast_results -i $isolate_file
	}
	popd > /dev/null #exit fimH
    }
    popd > /dev/null #exit 80_percent 
}
popd  > /dev/null #exit destination
printf "Done 80 threshold on fimH gene \n"

#bash /home/nagab/scripts/blast_mlst/blast_gene_mlst.sh -d < directory of all scaffolds file> -b /lustre/groups/pricelab/_files/st_gene_mlst/gyrA/db -o <output dir> -i <isolate file> -p 80

pushd $gene_output_location > /dev/null
{
    pushd $gene_output_location/80_percent > /dev/null
    {
	mkdir gyrA
	pushd  $gene_output_location/80_percent/gyrA> /dev/null
	{
	    mkdir blast_results
	    pushd $gene_output_location/80_percent/gyrA/blast_results >/dev/null
	    {
		
		bash /home/nagab/scripts/blast_mlst/blast_gene_mlst.sh -d $scaffolds_dir -b /lustre/groups/pricelab/_files/st_gene_mlst/gyrA/db -o $gene_output_location/80_percent/gyrA/blast_results -i $isolate_file -p 80
	    }
	    popd > /dev/null #exit blast_results
	    bash /home/nagab/scripts/blast_mlst/collect_blast_gene_results.sh -d $gene_output_location/80_percent/gyrA/blast_results -i $isolate_file
	}
	popd > /dev/null #exit gyrA
    }
    popd > /dev/null #exit 80_percent 
}
popd  > /dev/null #exit destination
printf "Done 80 threshold on gyrA gene \n"

#bash /home/nagab/scripts/blast_mlst/blast_gene_mlst.sh -d < directory of all scaffolds file> -b /lustre/groups/pricelab/_files/st_gene_mlst/CTX-M-15/db -o <output dir> -i <isolate file> -p 80

pushd $gene_output_location > /dev/null
{
    pushd $gene_output_location/80_percent > /dev/null
    {
	mkdir CTX-M-15
	pushd  $gene_output_location/80_percent/CTX-M-15> /dev/null
	{
	    mkdir blast_results
	    pushd $gene_output_location/80_percent/CTX-M-15/blast_results >/dev/null
	    {		
		bash /home/nagab/scripts/blast_mlst/blast_gene_mlst.sh -d $scaffolds_dir -b /lustre/groups/pricelab/_files/st_gene_mlst/gyrA/db -o $gene_output_location/80_percent/CTX-M-15/blast_results -i $isolate_file -p 80
	    }
	    popd > /dev/null #exit blast_results
	    bash /home/nagab/scripts/blast_mlst/collect_blast_gene_results.sh -d $gene_output_location/80_percent/CTX-M-15/blast_results -i $isolate_file
	}
	popd > /dev/null #exit CTX-M-15
    }
    popd > /dev/null #exit 80_percent 
}
popd  > /dev/null #exit destination
printf "Done 80 threshold on CTX-M-15 gene \n"

#bash /home/nagab/scripts/blast_mlst/blast_gene_mlst.sh -d < directory of all scaffolds file> -b /lustre/groups/pricelab/_files/st_gene_mlst/CTX-M_113/db -o <output dir> -i <isolate file> -p 80
#CTX-M_113

pushd $gene_output_location > /dev/null
{
    pushd $gene_output_location/80_percent > /dev/null
    {
	mkdir CTX-M_113
	pushd  $gene_output_location/80_percent/CTX-M_113> /dev/null
	{
	    mkdir blast_results
	    pushd $gene_output_location/80_percent/CTX-M_113/blast_results >/dev/null
	    {		
		bash /home/nagab/scripts/blast_mlst/blast_gene_mlst.sh -d $scaffolds_dir -b /lustre/groups/pricelab/_files/st_gene_mlst/gyrA/db -o $gene_output_location/80_percent/CTX-M_113/blast_results -i $isolate_file -p 80
	    }
	    popd > /dev/null #exit blast_results
	    bash /home/nagab/scripts/blast_mlst/collect_blast_gene_results.sh -d $gene_output_location/80_percent/CTX-M_113/blast_results -i $isolate_file
	}
	popd > /dev/null #exit CTX-M_113
    }
    popd > /dev/null #exit 80_percent 
}
popd  > /dev/null #exit destination
printf "Done 80 threshold on CTX-M_113 gene \n"

# -- 100 percent

pushd $gene_output_location > /dev/null
{
    
    mkdir $gene_output_location/100_percent
    pushd $gene_output_location/100_percent > /dev/null
    {
	mkdir parC
	pushd  $gene_output_location/100_percent/parC > /dev/null
	{

	    mkdir $gene_output_location/100_percent/parC/blast_results  
	    pushd $gene_output_location/100_percent/parC/blast_results > /dev/null
	    {
		bash /home/nagab/scripts/blast_mlst/blast_gene_mlst.sh -d $scaffolds_dir -b /lustre/groups/pricelab/_files/st_gene_mlst/parC/db -o $gene_output_location/100_percent/parC/blast_results -i $isolate_file -p 100
	    }
	    popd > /dev/null #exit blast_results
	    bash /home/nagab/scripts/blast_mlst/collect_blast_gene_results.sh -d $gene_output_location/100_percent/parC/blast_results -i $isolate_file
	}
	popd >/dev/null #exit parC 
    }
    popd > /dev/null #exit 100_percent
}
popd  > /dev/null #exit destination
printf "BLAST Done 100 threshold  parC\n"

pushd $gene_output_location > /dev/null
{
    pushd $gene_output_location/100_percent > /dev/null
    {
	mkdir fimH
	pushd  $gene_output_location/100_percent/fimH> /dev/null
	{
	    mkdir blast_results
	    pushd $gene_output_location/100_percent/fimH/blast_results >/dev/null
	    {
		
		bash /home/nagab/scripts/blast_mlst/blast_gene_mlst.sh -d $scaffolds_dir -b /lustre/groups/pricelab/_files/st_gene_mlst/fimH/db -o $gene_output_location/100_percent/fimH/blast_results -i $isolate_file -p 100
	    }
	    popd > /dev/null #exit blast_results
	    bash /home/nagab/scripts/blast_mlst/collect_blast_gene_results.sh -d $gene_output_location/100_percent/fimH/blast_results -i $isolate_file
	}
	popd > /dev/null #exit fimH
    }
    popd > /dev/null #exit 100_percent 
}
popd  > /dev/null #exit destination
printf "Done 100 threshold on fimH gene \n"

#bash /home/nagab/scripts/blast_mlst/blast_gene_mlst.sh -d < directory of all scaffolds file> -b /lustre/groups/pricelab/_files/st_gene_mlst/gyrA/db -o <output dir> -i <isolate file> -p 80

pushd $gene_output_location > /dev/null
{
    pushd $gene_output_location/100_percent > /dev/null
    {
	mkdir gyrA
	pushd  $gene_output_location/100_percent/gyrA> /dev/null
	{
	    mkdir blast_results
	    pushd $gene_output_location/100_percent/gyrA/blast_results >/dev/null
	    {
		
		bash /home/nagab/scripts/blast_mlst/blast_gene_mlst.sh -d $scaffolds_dir -b /lustre/groups/pricelab/_files/st_gene_mlst/gyrA/db -o $gene_output_location/100_percent/gyrA/blast_results -i $isolate_file -p 100
	    }
	    popd > /dev/null #exit blast_results
	    bash /home/nagab/scripts/blast_mlst/collect_blast_gene_results.sh -d $gene_output_location/100_percent/gyrA/blast_results -i $isolate_file
	}
	popd > /dev/null #exit gyrA
    }
    popd > /dev/null #exit 100_percent 
}
popd  > /dev/null #exit destination
printf "Done 100 threshold on gyrA gene \n"

#bash /home/nagab/scripts/blast_mlst/blast_gene_mlst.sh -d < directory of all scaffolds file> -b /lustre/groups/pricelab/_files/st_gene_mlst/CTX-M-15/db -o <output dir> -i <isolate file> -p 80

pushd $gene_output_location > /dev/null
{
    pushd $gene_output_location/100_percent > /dev/null
    {
	mkdir CTX-M-15
	pushd  $gene_output_location/100_percent/CTX-M-15> /dev/null
	{
	    mkdir blast_results
	    pushd $gene_output_location/100_percent/CTX-M-15/blast_results >/dev/null
	    {		
		bash /home/nagab/scripts/blast_mlst/blast_gene_mlst.sh -d $scaffolds_dir -b /lustre/groups/pricelab/_files/st_gene_mlst/gyrA/db -o $gene_output_location/100_percent/CTX-M-15/blast_results -i $isolate_file -p 100
	    }
	    popd > /dev/null #exit blast_results
	    bash /home/nagab/scripts/blast_mlst/collect_blast_gene_results.sh -d $gene_output_location/100_percent/CTX-M-15/blast_results -i $isolate_file
	}
	popd > /dev/null #exit CTX-M-15
    }
    popd > /dev/null #exit 100_percent 
}
popd  > /dev/null #exit destination
printf "Done 100 threshold on CTX-M-15 gene \n"

#bash /home/nagab/scripts/blast_mlst/blast_gene_mlst.sh -d < directory of all scaffolds file> -b /lustre/groups/pricelab/_files/st_gene_mlst/CTX-M_113/db -o <output dir> -i <isolate file> -p 80
#CTX-M_113

pushd $gene_output_location > /dev/null
{
    pushd $gene_output_location/100_percent > /dev/null
    {
	mkdir CTX-M_113
	pushd  $gene_output_location/100_percent/CTX-M_113> /dev/null
	{
	    mkdir blast_results
	    pushd $gene_output_location/100_percent/CTX-M_113/blast_results >/dev/null
	    {		
		bash /home/nagab/scripts/blast_mlst/blast_gene_mlst.sh -d $scaffolds_dir -b /lustre/groups/pricelab/_files/st_gene_mlst/gyrA/db -o $gene_output_location/100_percent/CTX-M_113/blast_results -i $isolate_file -p 100
	    }
	    popd > /dev/null #exit blast_results
	    bash /home/nagab/scripts/blast_mlst/collect_blast_gene_results.sh -d $gene_output_location/100_percent/CTX-M_113/blast_results -i $isolate_file
	}
	popd > /dev/null #exit CTX-M_113
    }
    popd > /dev/null #exit 100_percent 
}
popd  > /dev/null #exit destination
printf "Done 100 threshold on CTX-M_113 gene \n"


