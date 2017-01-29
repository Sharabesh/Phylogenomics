#!/usr/bin/python
import os
from Bio import Phylo



# Parameters need to be adjusted. 
# http://sco.h-its.org/exelixis/resource/download/NewManual.pdf
# s: input file name, m: model name, p: parsimony seed, n:output name
# w: output directiory 

out_dir = os.getcwd() + os.sep + 'ML_tree'
os.mkdir(out_dir)
os.system('.//standard-RAxML-master//raxmlHPC -s masked.fasta -m PROTCATDAYHOFF \
            -p 1000 -n ML_out -# 10 -w ' + out_dir)


ML_tree = Phylo.read(out_dir + '/RAxML_bestTree.ML_out', 'newick')



