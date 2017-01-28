#!/usr/bin/python



from Bio.Phylo.Applications import RaxmlCommandline
import os
from Bio import Phylo



# Parameters need to be adjusted. 
# http://sco.h-its.org/exelixis/resource/download/NewManual.pdf


raxml_cline = RaxmlCommandline(sequences='masked.fasta', model="GTRCAT",
                               name ='ML_out', num_replicates='10', 
                               working_dir=os.getcwd() + os.sep + 'ML_tree')

print (raxml_cline)

ML_tree = Phylo.read('ML_tree/RAxML_bestTree.ML_tree', 'newick')



