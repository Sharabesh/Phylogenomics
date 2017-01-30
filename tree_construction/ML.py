#!/usr/bin/python
import os
from Bio import Phylo



# Parameters need to be adjusted. 
# http://sco.h-its.org/exelixis/resource/download/NewManual.pdf
# s: input file name, m: model name, p: parsimony seed, n:output name
# w: output directiory 

os.system('raxmlHPC -s masked.fasta -m PROTCATDAYHOFF \
            -p 1000 -n ML_out -# 10 -w ' + 'RAxML_bestTree.ML_out')


ML_tree = Phylo.read('RAxML_bestTree.ML_out', 'newick')



