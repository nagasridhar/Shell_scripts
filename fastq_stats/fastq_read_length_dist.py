#!/usr/bin/env python


# Get read length numbers:
# mean length, min lenght, max lengh of fastq
# median length, std deviation
# get Read count, and reads with N in them

def get_read(read_format,read_file):
#{
    """
    Return list with read length 
    Also, return dict of Number of Ns found in Reads
    print total number of reads
    print reads with N in them
    """
    import sys
    read_length=[]
    n_count_dict={}
    format_list=[".fastq.gz",".fq",".fastq"]
    # will support these files only
    # won't read fastq.gz files
    
    if format_list.index(read_format) != 0:

        number_reads=0
        number_reads_with_n=0
        """
        Will use only cat to get read lengths
        """
        with open(read_file) as read_h:
            num_row=0 # get only line with ReaDs
            for line in read_h:
                
                num_row+=1
                if(num_row%4 == 2):
                    
                    number_reads+=1 # get number of reads
                    n_count=find_n(line.rstrip())
                    read_length.append(len(line.rstrip()))

                    if n_count>0:
                        number_reads_with_n+=1
                    #get read count with Number of Ns in it
                    
                    if n_count in n_count_dict:
                        n_count_dict[n_count]+=1
                    else:
                        n_count_dict[n_count]=1
                        
                #if ends for read's line
            #-- for loop ends for reading..
        #-- with loop ends ------

    else:
      print "We do not support compressed file dude. Hoom!"
      sys.exit()

    #for key in n_count_dict:
        #print key,"\t",n_count_dict[key]
    print "Total reads ", number_reads
    print "Reads with N in them",number_reads_with_n
    
    return n_count_dict,read_length

#} -------------------------------
def find_n(read):

    return read.count('N')
#}-----------------------------------
def check_input(fastq_file):
    """
    If input is not correct exit script
    If present then get full/absolute path
    """
        
    import os
    fastq_abs_path="" #store abs path
    
    if os.path.isfile(fastq_file):
        fastq_abs_path=os.path.abspath(fastq_file)
        return fastq_abs_path
    
    else:
        """
        Exit script if file not available
        """
        print "Incorrect fastq file"
        sys.exit()
        
#} ------------------------

def check_format(fastq_file):

    fastq_format="" #store file format and return
    
    """
    Check file format If not .gz, fastq or fq exit
    """
    
    if(fastq_file.endswith('.fastq.gz')):
        fastq_format= ".fastq.gz"

    elif (fastq_file.endswith('.fastq')):
        fastq_format=".fastq"

    elif (fastq_file.endswith('.fq')):
        fastq_format=".fq"
    else:
        print "File format unknown"
        sys.exit()
    
    return fastq_format #return the format

#}-----------------------------------
    
if __name__=="__main__":
    import argparse
    import numpy as np
    
    parser=argparse.ArgumentParser("description")
    parser.add_argument ('-f', '--fastq',required=True) # store fastq file
    args_dict = vars(parser.parse_args()) # make them
    
    fastq_file=check_input(args_dict['fastq'])
    fastq_file=check_input(fastq_file)
    fastq_format=check_format(fastq_file)

    n_dict,read_length_list=get_read(fastq_format,fastq_file)
    
    print "Max length is ",max(read_length_list)
    print "Min length is ",min(read_length_list)
    print "Std deviation is ",np.std(read_length_list)
    print "Median length is ",np.median(read_length_list)
    print "Mean length is ",np.mean(read_length_list)
    
