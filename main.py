import Bio
import sys
from Bio import SwissProt as sp
from  Bio import SeqIO
from urllib.parse import urlparse
import urllib
import re
from Bio.SubsMat.MatrixInfo import blosum62 as bs
import sys
import os

#BLAST XML Reader
from Bio.Blast import NCBIXML

#BLAST Command Line Tools (Python Wrapper)
from Bio.Blast.Applications import NcbiblastpCommandline

#Ensures environment variables are set
assert os.environ.get("BLASTDB")!=None, "Directory of Blast database must be set"

#Set path to BLAST home (Default install location)
if not os.environ["PATH"]:
    os.environ["PATH"] = "$HOME/ncbi-blast-2.2.29+/bin"


#Yields a list of topological annotations for each item in the text file
class Residue:
    def __init__(self,annotation,aa):
        self.annotation = annotation
        self.aa = aa
#Data structure to represent the a sequence of Residues/indels
class Sequence:
    def __init__(self,identifier):
        self.sequence = [] #Linked list may have been better data structure, allows for insertion and deletion in constant time
        self.id = identifier
    def __len__(self):
        return len(self.sequence)
    def __getitem__(self, item):
        return self.sequence[item]
    def __repr__(self):
        output = ""
        for item in self.sequence:
            output += item.aa
        return output

"""Uses BioPython rather than System os"""
def gather_homologs(fasta): #File location passed as input ("fasta.txt")
    blasted = NcbiblastpCommandline(query=fasta,db="nr.61",evalue=1,outfmt=5,out="alignment.xml")
    blasted()


"""Using System os
-qcov_hsp_perc allows for a defined query coverage """
def sys_gather_homologs(fasta):
    generic = "blastp -db {0} -query {1} -evalue {2}" \
              " -outfmt {3} -out {4} -qcov_hsp_perc {5}"
    db="nr.61"
    query = fasta
    evalue = 0.01 #Change if wanted
    outfmt = 5 # Indicates BLAST XML to write to file
    out = "alignment.xml" #Data written to alignment.xml
    query_coverage = 80 #Minimum required overlap
    formatted = generic.format(db,query,evalue,outfmt,out,query_coverage)
    os.system(formatted)


def parse_xml(alignment_file): #This is alignment.xml
    record_handle = open(alignment_file)
    records = NCBIXML.parse(record_handle)
    return records

def generate_msa(records_lst): #Takes output from parse_xml
    pass