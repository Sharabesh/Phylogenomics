# Phylogenomic Protein Function Prediction
This project seeks to implement a full phylogenetic pipeline to serve as a tool for transmembrane helix prediction. A target sequence and neighboring homologs are annotated by TMHMM and RAxML's maximum likelihood constructor is used to determine contribution scores to each sequence. See Presentation files for a brief overview of core functionality. 

# Implementing a full phylogenomic pipeline 
Given an input sequence this software gathers the closest 100 homologs, generates a multiple sequence alignment, masks that alignment, and then uses RAxML's maximum likelhood estimator to generate a phylogenetic tree. 

# Annotation Transfer Protocol 
This software uses an annotation transfer protocol based on evolutionary distances between proteins to transfer TMH annotations. That is, the more closely related a hit is to the target protein the more it's annotation at a particular site would matter in the determination of the target's true annotation. 
