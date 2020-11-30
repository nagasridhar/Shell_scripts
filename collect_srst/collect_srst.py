#!/usr/bin/env python

import sys, os, re, argparse,glob
import subprocess as sp

counter=0
def read_file(isolate, direc):
#isolate file and mlst_directory
    with open(isolate,'r') as fh:
        for line in fh:
            counter=0
            line=line.rstrip()
            
            #filename=direc+line+"__mlst__Klebsiella_pneumoniae__results.txt"
            #filename=direc+line+"__mlst__Escherichia_coli#1__results.txt"
            filename=direc+line+"__mlst__Staphylococcus_aureus__results.txt"
            if os.path.isfile(filename):

                with open (filename,'r') as mlst:
                    for mlst_line in mlst:
                        counter=counter+1
                        mlst_line=mlst_line.rstrip()
                        if counter == 2:
                            array=re.split('\s+',mlst_line)

                            mlst_line=line+"\t"
                            try:
                                for i in range(1,len(array)):
                                    mlst_line+=array[i]+"\t"
                            except IndexError as e:
                                print >> sys.stderr, e
                            finally:
                                print mlst_line
                        #get second line of mlst definition
                #with loop ends of mlst file

        #with loop ends of isolate file
            else:
                print >> sys.stderr, filename, "Have issues with the file. Either missing or some unbeknown stuff"
                #print filename, "Please check the file extension for SRST: ecoli, saureus mlst__Escherichia_coli#1__results.txt"
                #sys.exit()
                
            if counter==1:
                print line+'\t'+"Didn't get any MLST definitions"
                                
###         
if __name__ == "__main__":
    parser=argparse.ArgumentParser("description")
    parser.add_argument ('-i', '--isolate',required=True) # store isolate file
    parser.add_argument ('-d','--dir',required=True) # store directory where all mlst_.txt are present
    args_dict = vars(parser.parse_args()) # make them dict..

    read_file(args_dict['isolate'],args_dict['dir'])
