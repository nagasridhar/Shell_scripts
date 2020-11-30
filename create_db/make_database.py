#!/usr/bin/env python


import subprocess
import os, argparse, re,sys
from subprocess import call
from Bio import SeqIO
#
if sys.version_info < (2,7):
    print "Python 2.7+ are needed for script"
    exit(1)
##
def check_inputs(**temp_dict):
#function begins

    temp_out_dir=""
    temp_fasta_file=""
    #initiate variables
    
    if os.path.isfile(temp_dict['fasta_file']):
        temp_fasta_file=os.path.abspath(temp_dict['fasta_file'])
    else:
        print "incorrect fasta file"
        sys.exit(1)
    if os.path.exists(temp_dict['out_dir']):
        temp_out_dir=os.path.abspath(temp_dict['out_dir'])
    else:
        print "Bad output directory"
        sys.exit(1)
    return temp_out_dir,temp_fasta_file #send complete path
#function ends
def find_duplicates(temp_fasta):
#{
    duplicates=False
    total=0
    seq_id=[] # list that stores sequence ids
    
    for record in SeqIO.parse(temp_fasta,"fasta"):
        if not record.id in seq_id:
            seq_id.append(record.id)
        else:
            duplicates=True
            print record.id ," is duplicate"
        total+=1
        
    print total ," input sequences "

    if duplicates:
        print "Cannot make db "
        sys.exit()
#}
def create_db(temp_fasta_file):
#function begins
    """
    this function creates database for the FASTA file input
    All database have "fasta" in them

    """
    subprocess.call(["makeblastdb -in %s -dbtype nucl -out %s" %( temp_fasta_file,temp_fasta_file)],
                    shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # subprocess waits until the initiated process ends.. !!
    
    #create database..
#function ends

def parse_fasta(temp_fasta,temp_out):
    #fasta file and output dire are input
#function begins
    
    seq_dict=SeqIO.parse(open(temp_fasta,"r"),"fasta") #parse Biopython
    for record in seq_dict: #iterate over the seqs

        subprocess.call(["mkdir %s" %(temp_out+"/"+record.id)],shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        """
        Individual database folder is created
        """
        fasta_file_location=temp_out+"/"+record.id #directory where fasta will be created
        temp_file_name=fasta_file_location+"/"+record.id+".fasta" #fasta file
        with open(temp_file_name,'w') as out_handle:
            """
            create new fasta file for the sequen id.
            and then create database for it
            """
            out_handle.write(">"+record.id+"\n")
            out_handle.write(str(record.seq)+"\n")
        create_db(temp_file_name)
        
#function ends
if __name__=="__main__":

    #make a dict of arguments
    parser=argparse.ArgumentParser("description")
    parser.add_argument ('-f', '--fasta_file', help='location for fasta file',required=True) # store input fasta file
    parser.add_argument ('-o', '--out_dir', help='output direct',required=True) # store output direc
    args_dict = vars(parser.parse_args()) # make them dict..
    
    #send to verfiy them
    out_dir,fasta_file=check_inputs(**args_dict) #returns name if exists
    find_duplicates(fasta_file)
    """
    Exit script if there are duplicates
    """
    parse_fasta(fasta_file,out_dir)
    print "Done with the script. Exiting!"
    print "<--Bye-Bye-->"
