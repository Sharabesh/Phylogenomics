#!/usr/bin/python


import os
import Bio
from Bio import Phylo
from Bio import AlignIO
import subprocess
from subprocess import Popen, PIPE


recs_file = os.getcwd() + '/masked.fasta' # this needs to be removed; 


"""Uses Bio's AlignIO module to convert fasta into relaxed phylip"""
def fasta_to_phylip():
    
# Phylip seems to be strict about the inputfile format 
# to keep the full accession, relaxed phylip is needed 
# but phylip doesn't accept the input file if the accession has "." 
# so I tried to overwrite the input file ater removing ".#" 
    
    with open(recs_file, "r+") as f:
        input_file = f.readlines()
        for a in range(0,len(input_file)):
            if "." in input_file[a]:
                input_file[a] = input_file[a].split(".", 1)[0]
        f.close()
    with open(recs_file, "w+") as f:
        for item in input_file:
            f.write("%s\n" % item)
        f.close()
                
            
    in_file = recs_file
    out_file = "recs_phylip"
    in_format = "fasta"
    out_format = "phylip-relaxed"
    
    Bio.AlignIO.convert(in_file, in_format, out_file, out_format)



"""Uses PhyML to generate ML tree which is read out by Bio's Phylo module"""
def generate_MLtree():
    directory = os.getcwd()
    recs_phylip = directory + '/recs_phylip'

# I wasn't sure what this line below does; it may need to be fixed if needed
#    os.system("rm {0}/*".format(directory + '/operating_reqs/tree_files'))


# The initial directory and executable: ./PhyML-3.1/PhyML3.1_linux62 
# the directory might need to be changed 
    command = 'PhyML-3.1/PhyML -i {0} -d {1} -m {2}'
    sequence_file_name = recs_phylip # input is in phylip format, not fasta
    data_type = "aa" 
    substitution_model = "LG" 
# the output seems to be created in the current working directory 
    os.system(command.format(sequence_file_name,
        data_type,substitution_model))

# the output name comes out as 'inputname_phyml.tree.txt'       
    tree_file = directory + '/' + 'recs_phylip_phyml_tree.txt'    
    tree = Phylo.read(tree_file, 'newick')
    Phylo.draw_ascii(tree)
    return tree 


"""Uses Phylip to generate NJ tree which is read out by Bio's Phylo module"""
def generate_NJtree():
# python does not wrap Phylip; a string needs to be passed
# to the excuted file     
    
# First, create the distance materix 
# the initial directory and executable: ./phylip-3.696/exe/protdist
    directory = os.getcwd()

# error pops up if there is any existing file that has names that phylip wants to use
# mainly, infile/outifle/outtree but I removed all the exisitng files 
    
    possible_files = ['infile', 'MSA_phylip', 'outfile', 'NJ_tree', 'DistanceMatrix']
    for A in possible_files:
        if os.path.isfile(directory + '/' + A) == 1:
            os.remove(directory + '/' + A)

# it seems Phylip only accets the file name 'infile'         
    os.rename(directory + '/recs_phylip', 'infile')
    sequence_file_name = directory + '/infile'

# The directory might need to be changed 
    p = subprocess.Popen([directory + "/phylip-3.696/exe/protdist", sequence_file_name],
                          stdin = subprocess.PIPE)
    p.communicate("Y")                 
# the distance matrix is created as outfile; below is to construct tree

    os.rename(directory + '/' + 'infile', 'MSA_phylip')
    os.rename(directory + '/' + 'outfile', 'infile')
    
    
# the initial directory and executable: ./phylip-3.696/exe/neighbor
# the directory may need to be changed

    input_matrix = directory + '/infile'
    p = Popen([directory + "/phylip-3.696/exe/neighbor", input_matrix],
                     stdin=PIPE)
    p.communicate("Y")
    
# the output name comes out as 'outfile'
    os.rename(directory + '/' + 'infile', 'DistanceMatrix')
    os.rename(directory + '/' + 'outtree', 'NJ_tree')      
    tree_file = directory + '/' + 'NJ_tree'    
    tree = Phylo.read(tree_file, 'newick')
    Phylo.draw_ascii(tree)
    return tree 




# This part is to test; needs to be removed

fasta_to_phylip()
generate_MLtree()


