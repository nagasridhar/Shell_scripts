#!/usr/bin/env python


import os, argparse, re,sys
from Bio import SeqIO
from collections import defaultdict #this is to set int value as 0
"""
this script finds out allele/gene count 
in the virulence database downloaded 
"""
#
if sys.version_info < (2,7):
    print "Python 2.7+ are needed for script"
    exit(1)
##
def check_inputs(**temp_dict):

    virulent_file=""
    if os.path.isfile(temp_dict['fasta_virulent_file']):
        virulent_file=os.path.abspath(temp_dict['fasta_virulent_file'])
    else:
        print "Incorrect virulent file"
        sys.exit(1)
        
    return virulent_file #send complete path
##
def find_colon(seq_id):
    """
    thid function sends only the gene name after trimming at : position
    """
    if (seq_id.find(":") == -1):
        print seq_id, " has issues"
        sys.exit(1)
    else:
        return seq_id[:seq_id.find(":")] #index->seq_id.find(":")

##----->function ends
def find_allele(temp_virul_file):

    allele_count=defaultdict(int)
    # a magic dictionary has been found which sets everythin to 0 by default
    
    seq_dict=SeqIO.parse(open(temp_virul_file,"r"),"fasta") #parse Biopython
    for record in seq_dict:
        """
        loop thorugh the seqs
        """
        trim_seq_id=find_colon(record.id) #seq id after trimming colon in stx2A:87:AY443052:a
        allele_count[trim_seq_id]+=1
    for key in allele_count:
        print key,"\t",allele_count[key]
    
##------> function ends
if __name__=="__main__":

    #make a dict of arguments
    parser=argparse.ArgumentParser("description")
    parser.add_argument ('-f', '--fasta_virulent_file', help='location for all virulence',required=True) # store input directory to fasta files
    args_dict = vars(parser.parse_args()) # make them dict..
    virul_file=check_inputs(**args_dict)
    find_allele(virul_file) #send virulen file
    
    print "Done with the script. Exiting!"
    print "<--Bye-Bye-->"
