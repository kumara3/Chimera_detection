'''
Created on Aug 23, 2014

@author: kumara3
This script reads the fasta file
'''
import os
from sys import argv

#script,pathname = argv
#fasta_file = os.path.abspath(pathname)

def readFA(fa_file):
    readfasta = {}
    with open(fa_file) as fp:
        info = fp.readline().rstrip()
        seq=""
        
        while True:
            line=fp.readline()
            if not line or line[0]==">":
                readfasta[info]=seq
                if not line:
                    break
                info = line.rstrip()
                seq = ""
                
            else:
                seq += line.strip()
    return readfasta

#return_readfasta = readFA(fasta_file)

