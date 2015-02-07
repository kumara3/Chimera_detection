'''
Created on Aug 20, 2014

@author: kumara3

Input:meta-velvet.contigs.fa file from MetaVelvet run.
Output:x,y coordinate.x is the bin length and y is number of contigs generated
in metagenomic assembly.
'''

from _bisect import bisect
from sys import argv
import os
import csv


script,pathname = argv
fasta_file = os.path.abspath(pathname)

if os.path.exists(fasta_file):
    print "Good to go........"
else:
    raise TypeError," Wrong Input File"



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

return_readfasta = readFA(fasta_file)

def writeContigs(return_fa,col_width=80):
    for i in range(1000,12000,1000):
        file_name = "bins_%i.fasta"%(i)
        j='greater_than_11000'
        last_file = "bins_%s.fasta"%(j)
        path = '/home/kumara3/workspace/Read_simulator/draw_bar_graph/bins'
        if not os.path.exists(path):
            os.mkdir(path)
    
        
        fp = open(os.path.join(path,last_file),'wb')
        with open(os.path.join(path,file_name),'wb') as fw:
    
            
            for item in return_fa:
                if len(return_fa[item]) < i and len(return_fa[item]) >= i-1000:
                    fw.write((item+"\n"+"\n".join(return_fa[item][i:i+col_width] for i in range(0,len(return_fa[item]),col_width)) +"\n"))
                if len(return_fa[item]) >= 11000:
                    fp.write((item+"\n"+"\n".join(return_fa[item][i:i+col_width] for i in range(0,len(return_fa[item]),col_width)) +"\n"))
                                                     
writeContigs(return_readfasta,col_width=80)
