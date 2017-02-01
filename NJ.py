#!/usr/bin/env python
import os
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO
from Bio import Phylo 


directory = os.getcwd()


MSA = AlignIO.read(directory + '/' + recs_file, 'fasta')


Distance_model = 'identity' # default model
                            #dna_models, protein_models,models possible

calculator = DistanceCalculator('identity')
distance_matrix = calculator.get_distance(MSA)

constructor = DistanceTreeConstructor()
NJ_tree = constructor.nj(distance_matrix) # perhaps needs to be saved as a file? 
Phylo.draw_ascii(NJ_tree)




