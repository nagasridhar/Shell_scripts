#!/usr/bin/env python


import sys, os, re, argparse
import os.path
import subprocess as sp
"""
This will take file(s) from SRST and BLAST sequence typing.
Please run using respective flags.
Output will be generated in the designated directory you mention
"""
#
if sys.version_info < (2,7):
    print "Python 2.7+ are needed for script"
    exit(1)
##
def calculate_stats(**temp_dict):
    
    mlst_count={}
    line_header=0
    out_dir,mlst_file=check_inputs(**temp_dict)
    seq_type="" # file showin count per Seq Type
    seq_type_dist=""

    if temp_dict['blast']:
    #{{ blast function  
        seq_type_dist=os.path.join(out_dir,'mlst_dist_blast.txt') # file showin count - for unique, not found, mutliple and anything above 1
        seq_type=os.path.join(out_dir,'mlst_counts_blast.txt') # file showin count per Seq Type
        
        with open(mlst_file,'r') as f:
            for line in f:
                line=line.rstrip()
                line_header+=1

                if line_header >1:                    
                    array=re.split('\s+',line)

                    if not len(array) == 9: # if there's an issue 
                        print "issue with ", line

                    else:                        
                        if array[8] == "(not_found)": # if not found
                            if not "not_found" in mlst_count:
                                mlst_count["not_found"]=1
                            else:
                                mlst_count["not_found"]+=1

                        elif array[8].find(",") > 0: # if has comma
                            if not "multiple" in mlst_count:
                                mlst_count["multiple"]=1
                            else:
                                mlst_count["multiple"]+=1

                        else: # everything else goes and added up
                            if not array[8] in mlst_count:
                                mlst_count[array[8]]=1
                            else:
                                mlst_count[array[8]]+=1
    ## }}end of blast function

    if temp_dict['srst']:
    #{{ SRST function

        seq_type_dist=os.path.join(out_dir,'mlst_dist_srst.txt') # file showin count - for unique, not found, mutliple and anything above 1
        seq_type=os.path.join(out_dir,'mlst_counts_srst.txt') # file showin count per Seq Type

        with open(mlst_file,'r') as f:
            for line in f:
                line=line.rstrip()
                line_header+=1

                if line_header >1:
                    array=re.split('\s+',line)
                    if not len(array) ==13:
                        print "issue with ",line

                    else:
                        st=array[1]
                        
                        if st == "NF" or st == "NF*":
                            if not "not_found" in mlst_count:
                                mlst_count['not_found']=1
                            else:
                                mlst_count['not_found']+=1

                        elif st[-1] == "*": # if last letter is an asterix
                            if not "allele_diff" in mlst_count:
                                mlst_count["multiple"]=1
                            else:
                                mlst_count["multiple"]+=1
                        else:
                            if not st in mlst_count:
                                mlst_count[st]=1
                            else:
                                mlst_count[st]+=1
                            
    #}} end of SRST function           
    
    with open (seq_type,'a') as fh:
 
       for key in mlst_count:            
        fh.write("%s\t%s\n" %(key,mlst_count[key]) ) # print key-value pair .  .  

    get_mlst_dist(mlst_count,seq_type_dist)
        
###
def get_mlst_dist(temp_mlst,file_name):
#{{{function starts
    uniq=0
    not_found=0
    multiple=0
    
    for key in temp_mlst:
        if not key == "multiple" and temp_mlst[key] ==1:
            uniq+=1
        elif not key == "not_found" and temp_mlst[key] ==1:
            uniq+=1
        elif key == "not_found":
            not_found+=1
        elif key == "multiple":
            multiple+=1
        print "Value is mutple ",multiple

    with open(file_name,'a') as fh:
        #if not multiple == 0:
      #  fh.write("%s\t%s\n" %("multiple",str(1) ) )
        #else:
        print "mutliple as ", multiple
        fh.write("%s\t%s\n" %("multiple",str(multiple) ) )
        #if not not_found == 0 :
        #    fh.write("%s\t%s\n" % ("not_found",str(1) ) )
        #else:
        fh.write("%s\t%s\n" % ("not_found",str(not_found) ) )
        fh.write("%s\t%s\n" % ("Unique",str(uniq) ))

        for key in temp_mlst:
           # if key == "multiple":
            #    if temp_mlst["multiple"] ==1:
             #       fh.write("%s\t%s\n" %(key,str(1) ) )

            #elif key == "not_found":
             #   if temp_mlst["not_found"] == 1:
              #      fh.write("%s\t%s\n" %(key,str(1) ) )

            if not temp_mlst[key] ==1 and not key == "not_found" and not key == "multiple" :
                fh.write("%s\t%s\n" %(key,temp_mlst[key] ) )
            #elif not temp_mlst[key] == 1 and not key == "not_found" and not key == "multiple":
             #   fh.write("%s\t%s\n" %("Unique",str(uniq) ) )



#}}}function ends
##
def check_inputs(**temp_dict): # do sanity check on inputs
    blast_file=""
    srst_file=""
    out_dir=""

    if not temp_dict['blast'] and not temp_dict['srst']:                                                                                                                                                           
        print "BLAST or SRST file is required"                                                                                                                                                                     
        sys.exit()   

    if temp_dict['blast'] and temp_dict['srst']:
        print "Either BLAST or SRST output is required. Not both"
        sys.exit()

    if temp_dict['blast']:
        blast_file=os.path.abspath(temp_dict['blast'])
        if not os.path.exists(blast_file):
            print "Incorrect blast file provided"
            sys.exit()
        
    if temp_dict['srst']:
        srst_file=os.path.abspath(temp_dict['srst'])
        if not os.path.exists(srst_file):
            print "Incorrect SRST file provided"
            sys.exit()

    if temp_dict['out']:
        out_dir=os.path.abspath(temp_dict['out'])
        if not os.path.exists(out_dir):
            print "Output directory is incorrect"
            sys.exit()

    if temp_dict['srst']:
        return out_dir,srst_file

    if temp_dict['blast']:
        return out_dir,blast_file
##    
if __name__=="__main__": 
#{{{ main starts
    
    parser=argparse.ArgumentParser("description")
    parser.add_argument ('-b', '--blast', help='BLAST ST file location/path',required=False) # store blast file output
    parser.add_argument ('-s', '--srst', help='SRST ST file location/path',required=False) # store srst file output
    parser.add_argument ('-o', '--out', help='output directory',required=True) # store output directory
    args_dict = vars(parser.parse_args()) # make them dict.. 

    calculate_stats(**args_dict)
    
    
    print "<--Done calculating MLST counts!-->"

#}}} main ends
