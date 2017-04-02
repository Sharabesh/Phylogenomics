import urllib
import urllib.request
from Bio import SwissProt as sp
from Bio import SeqIO
import os
import subprocess
from main import *


#Using Belvu select only those with specified percent identity
THIRTY_PERCENT_FILE = "lessThan30-40.fasta"
FOURTY_PERCENT_FILE = "aligned40-70Percent.fasta"
SEVENTY_PERCENT_FILE = "MoreThan70Percent.fasta"

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)



LINK = "http://www.uniprot.org/uniprot/{0}.txt"

def process(files): #the file is a text file of a swissprot protein
    parsed = sp.parse(files)
    record = next(parsed)
    sequence = record.sequence
    length = len(sequence)
    a = record.features
    finished = False
    output = []
    for item in a:
        if finished:
            break
        else:
            if(not(item[0] == "CHAIN") and not(item[0] == "DOMAIN")): #this requires more refinement to allow for a range of domains
                output.append([item[0],item[1],item[2]])
                if item[2] == length:
                    finished = True
    return output


def get_target_data():
    send_homologs()
    generate_msa()
    textrecord = urllib.request.urlopen(LINK.format("KCNA1_HUMAN"))
    textrecord = process(textrecord)
    helices = [x for x in textrecord if x[0] == "TRANSMEM"]
    for i in range(len(helices)):
        helices[i].append(helices[i][0])
        helices[i].pop(0)
        helices[i] = tuple(helices[i])
    processed = SeqIO.parse(recs_file,"fasta")
    target = next(processed)
    identifier = target.name
    sequence = str(target.seq)
    representation = generate_rep(sequence, helices, identifier)
    return representation


def num_transmembrane(target):
    count = 0
    for i in range(len(target)):
        if target[i].annotation == "H":
            count += 1
    return count



def use_consensus(msaFile=recs_file):
    annotated_target = generate_consensus_withoutTree(get_target_data().id, msaFile)

    correct_target = get_target_data()

    for i in range(len(annotated_target)):
        if annotated_target[i].annotation == 'i' or annotated_target[i].annotation == 'o':
            annotated_target[i].annotation = "H"
    for i in range(len(correct_target)):
        if correct_target[i].annotation == "TRANSMEM":
            correct_target[i].annotation = "H"
    agreements = 0
    for i in range(len(annotated_target)):
        if correct_target[i].annotation == annotated_target[i].annotation:
            agreements += 1

    return (agreements/len(correct_target)) * 100








