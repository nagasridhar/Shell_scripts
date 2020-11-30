#!/usr/bin/env python

            
long_description="""
Takes in curated virulence database file as an input
creates unique file per gene (combining all alleles)
finally creates database for them

"""

import subprocess
import os, argparse, re,sys
from subprocess import call
from Bio import SeqIO

#
if sys.version_info < (2,7):
    print "Python 2.7+ are needed for script"
    exit(1)
##
def parse_fasta(temp_fasta,temp_out):

    seq_ident=[]
    for record in SeqIO.parse(temp_fasta,"fasta"):
        
        colon_loc=(record.id).find(":") # get colon's position
        
        gene=record.id[:colon_loc] # get gene name..
        #get rid of allele info - substring
        
        if gene not in seq_ident:
            #new gene../allele
            seq_ident.append(gene)

            subprocess.call(["mkdir %s" %(temp_out+"/"+gene)],shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            #create directory for the gene/allele

        fasta_file_location=temp_out+"/"+gene #directory where fasta will be created
        temp_file_name=fasta_file_location+"/"+gene+".fa" #fasta file

        with open(temp_file_name,'a') as out_handle:
            out_handle.write(record.format("fasta"))
        #with handle closes. write sequence to that file
        
    # for loop end for reading FASTA file
    
    print len(seq_ident)
    print "Done through creating FASTA files"
    return seq_ident
#<<-- function ends-------------->>>

def create_db(temp_list,temp_out):
    #temp_list is gene list
    
    for g in temp_list:
        #g is gene name
        
        gene_location=temp_out+"/"+g
        temp_fasta=gene_location+"/"+g+".fa"

        subprocess.call(["makeblastdb -in %s -dbtype nucl -out %s" %(temp_fasta,temp_fasta)],
                            shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #for loop ends
    print "Done through making BLAST DB"
    
#<<---Function ends---------------->>>>>    
if __name__=="__main__":
    parser=argparse.ArgumentParser("description")
    parser.add_argument ('-f', '--fasta_file', help='location for fasta file',required=True) # store input fasta file
    parser.add_argument ('-o', '--out_dir', help='output direct',required=True) # store output direc
    
    args_dict = vars(parser.parse_args()) # make them dict..

    # os.path.abspath(args_dict['out_dir']),os.path.abspath(args_dict['fasta_file']) - get full path
    abs_fasta_file=os.path.abspath(args_dict['fasta_file'])
    abs_out_dir=os.path.abspath(args_dict['out_dir'])
    gene_list=parse_fasta(abs_fasta_file,abs_out_dir)
    create_db(gene_list,abs_out_dir)
