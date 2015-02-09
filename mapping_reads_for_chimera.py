'''
Created on Nov 13, 2014

@author: kumara3

'''
import re
import read_fasta_file

user_string1 = "/home/kumara3/workspace/Read_simulator/input_data/Metasim_MC_92S_24M/taxon_profile-Empirical.994f05d5.fna"

dict_reads = read_fasta_file.readFA(user_string1)   

map_input_read = {}
def read_mapping(new_dict_reads):
    for each in new_dict_reads:
        my_match = re.match(r'>(r\d+\.\d+)\s+\|.*\|SOURCE_\d+\=\"([a-zA-Z0-9]*).*\".*',each)
        id = my_match.group(1)
        organism=my_match.group(2)
        map_input_read[id]=organism
    return map_input_read
        
def bwa_output(bwa_input_file):
    map_bwa_reads = {}
    with open (bwa_input_file,'r') as fp:
        for each_item in fp:
            my_match = re.match(r'^(r\d+\.\d+).*(NODE.*cov\_\d+?\.\d+).*',each_item)
            if my_match:
                r_id = my_match.group(1)
                contigs = my_match.group(2)
                if contigs not in map_bwa_reads:
                    map_bwa_reads[contigs]=[r_id]
                else:
                    map_bwa_reads[contigs].append(r_id)
        return map_bwa_reads
    
def comparision(new_map_input_read,new_map_bwa_reads):
    for k,v in new_map_bwa_reads.iteritems():
        unmapped_reads=[]
        mapped_reads=[]
        chimeric_reads={}
        clean_reads={}
        chimeric_reads['chimera']={}
        clean_reads['Non_chimeric']={}
        first_organism = new_map_input_read[v[0]]
        clean_reads['Non_chimeric'][k]=[first_organism]
        chimeric_reads['chimera'][k] = [first_organism]
        counter = 0
        number_chimera_reads=0
        unmapped_reads.append(first_organism)
        mapped_reads.append(first_organism)
        for i in range(1,len(v)):  

            if new_map_input_read[v[i]] == first_organism:
                counter += 1
                mapped_reads.append(new_map_input_read[v[i]])

            if new_map_input_read[v[i]] != first_organism:
                unmapped_reads.append(new_map_input_read[v[i]])
        
        if counter == len(v)-1:
            print "The Non chimeric reads is %s"%(k),counter,len(v)-1
                

        if counter < len(v)-1:
            number_of_chimera_reads=len(v)
            print "The chimeric contig is %s"%(k),counter,len(v)-1,len(v)
            
user_string2 = '/home/kumara3/workspace/Read_simulator/input_data/Metasim_MC_92S_24M/hash_length21/first_output.sam'
new_map_input_read=read_mapping(dict_reads)
new_map_bwa_reads = bwa_output(user_string2)
comparision(new_map_input_read,new_map_bwa_reads)




                 
    
    
    
                       
