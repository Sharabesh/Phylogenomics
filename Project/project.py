import Bio
from Bio import SeqIO
from Bio.SubsMat.MatrixInfo import blosum62 as bs
import numpy as np
import csv

"""PROGRAMMING CONSTANTS"""
recs_file = "homologs_with_functional_sites.fasta"
MAIN = [] # TO hold matrix of aligned sequences

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


def msaprocess(): #Takes an input from the generate_msa function
    processed = SeqIO.parse(recs_file,"fasta")
    for record in processed:
        identifier = record.name
        target = Sequence(identifier)
        gapped_sequence = str(record.seq)
        for item in gapped_sequence:
            target.sequence.append(Residue("",item))
        MAIN.append(target)
    print("Done!")

def scoreAcids(residue1, residue2):
    try:
        return bs[(residue1,residue2)]
    except:
        try:
            return bs[(residue2,residue1)]
        except:
            return 0 #Score of a gap character is set to 0

#BLOSUM Formula
"""Sum of BL62 scores over all residue pairs in a column / #residue pairs"""
def calc_blosum_col_scores():
    total_scores = []
    for column in range(len(MAIN[0])): #Iterate over the columns of a protein
        proteins = []
        for protein in range(len(MAIN)):
            proteins.append(MAIN[protein][column])
        #Iterate over proteins to get the BLOSUM Score
        score = 0
        pairs = 0
        for i in range(len(proteins)-1):
            for j in range(i + 1,len(proteins)):
                score += scoreAcids(proteins[i].aa,proteins[j].aa)
                # if proteins[i].aa != "-" and proteins[j].aa != "-":
                pairs += 1
        total_scores.append(score/pairs)
    return total_scores


#Pairwise Identity Formula
"""#exact matches / number of ungapped alignment columns (no gaps allowed)"""
def calc_pairwise_identity(seq1,seq2): # 'XP_019815713.1/1-495','XP_019737311.1/1-494'
    targets = [x for x in MAIN if x.id==seq1 or x.id==seq2]
    matches = 0
    total = 0
    for column in range(len(targets[0])):
        first = targets[0][column]
        second = targets[1][column]
        if first.aa == second.aa and first.aa != "-":
            matches += 1
            total +=1
        elif first.aa != second.aa and first.aa != "-" and second.aa != "-":
            total += 1
    return matches/total

def calc_pairwise_identity_with_gaps(seq1,seq2):
    targets = [x for x in MAIN if x.id==seq1 or x.id==seq2]
    matches = 0
    total = 0
    for column in range(len(targets[0])):
        first = targets[0][column]
        second = targets[1][column]
        if first.aa == second.aa and first.aa != "-":
            matches +=1 
            total +=1 
        else:
            total +=1 
    return matches/total 


def ent(lst):
    prob_dict = {x:lst.count(x)/len(lst) for x in lst}
    probs = np.array(list(prob_dict.values()))
    return -np.sum(probs * np.log2(probs))


def calc_column_entropy():
    total_entropy = []
    for column in range(len(MAIN[0])):
        acids = []
        for protein in range(len(MAIN)):
            acids.append(MAIN[protein][column].aa)
        total_entropy.append(ent(acids))
    return total_entropy

def calc_average_id():
    total_id = [] 
    for column in range(len(MAIN[0])):
        ids = [] 
        for protein in range(len(MAIN)):
            ids.append(MAIN[protein][column])

        score = 0 
        pairs = 0
        for i in range(len(ids)-1):
            for j in range(i + 1,len(ids)):
                score += 1 if ids[i].aa == ids[j].aa else 0 
                # if proteins[i].aa != "-" and proteins[j].aa != "-":
                pairs += 1
        total_id.append(score/pairs)
    return total_id 




#Calculate Query Coverage
""" (hit_length - gap_open) / query_length"""
def calc_query_coverage(target,seq2):
    main = len([x for x in MAIN if x.id==target][0])
    secondary = [x for x in MAIN if x.id==seq2][0]
    secondary_gaps = 0
    for item in secondary:
        if item.aa == "-":
            secondary_gaps += 1
    return (len(secondary) - secondary_gaps) / main

def calc_coverage():
    coverages = [] 
    for column in range(len(MAIN[0])):
        acids = [] 
        for protein in range(len(MAIN)):
            acids.append(MAIN[protein][column])
        gaps = len([x for x in acids if x.aa == "-"])
        coverages.append(gaps/len(acids))
    return coverages 


def min_percent_identity():
    min_id = ('_','_',100)
    max_id = ('_','_',0)
    all_ids = []
    for i in range(len(MAIN)-1):
        for j in range(i + 1,len(MAIN)):
            x = calc_pairwise_identity(MAIN[i].id,MAIN[j].id) * 100
            if x < min_id[2]:
                min_id = (MAIN[i].id,MAIN[j].id,x)
            if x > max_id[2]:
                max_id = (MAIN[i].id,MAIN[j].id,x)
            all_ids.append(x)
    print("MIN ID: ", min_id)
    print("MAX ID: ", max_id)
    return (np.average(all_ids),min_id,max_id)

def calc_identities_with_gaps():
    min_id = ('_','_',100)
    max_id = ('_','_',0)
    all_ids = [] 
    for i in range(len(MAIN)-1):
        for j in range(i+1,len(MAIN)):
            x = calc_pairwise_identity_with_gaps(MAIN[i].id,MAIN[j].id) * 100 
            if x < min_id[2]:
                min_id = (MAIN[i].id,MAIN[j].id,x)
            if x > max_id[2]:
                max_id = (MAIN[i].id,MAIN[j].id,x)
            all_ids.append(x)
    return (np.average(all_ids),min_id,max_id)


def generate_stats():
    msaprocess()
    with open("stats.txt","w") as file:
        (average_id,min_id,max_id) = min_percent_identity()

        file.writelines("Average ID is: " + str(average_id) + "%\n")
        file.writelines("Min Percent ID BETWEEN: " + str(min_id) + "\n")
        file.writelines("Max Percent ID BETWEEN: " + str(max_id) + "\n")
        file.writelines("The following show the Percent Identities including gaps\n")
        (average_id,min_id,max_id) = calc_identities_with_gaps()
        file.writelines("Average ID is: " + str(average_id) + "%\n")
        file.writelines("Min Percent ID BETWEEN: " + str(min_id) + "\n")
        file.writelines("Max Percent ID BETWEEN: " + str(max_id) + "\n")
    print("Stats are in stats.txt")
    return


def generate_csv():
    count = 0
    msaprocess()
    #Write Blosum Scores
    with open("BL62_score.csv","w") as file:
        for item in range(len(MAIN[0]) - 1):
            file.writelines(str(item) + ',')
        file.writelines(str(len(MAIN[0])-1)) #Avoid comma in last step
        file.writelines("\n") #Finish enumeration
        x = calc_blosum_col_scores()
        for item in x:
            file.writelines(str(item) + ',')
        count += 1

    msaprocess()
    #Write percent identities relative to query sequence
    with open("percent_ids.csv","w") as file:
        for sequence in MAIN:
            file.writelines(sequence.id + ',')
        file.writelines('\n')

        #Assumes target is first sequence
        target = MAIN[0].id
        for sequence in MAIN:
            id = calc_pairwise_identity(target,sequence.id)
            file.writelines(str(id) + ',')
        count += 1

    msaprocess()
    #Write Entropies
    with open("entropy.csv","w") as file:
        for item in range(len(MAIN[0]) - 1):
            file.writelines(str(item) + ',')
        file.writelines(str(len(MAIN[0])-1)) #Avoid comma in last step
        file.writelines("\n") #Finish enumeration
        x = calc_column_entropy()
        for item in x:
            file.writelines(str(item) + ',')
        count += 1

    msaprocess()
    #Write coverages
    with open("coverage.csv","w") as file:
        for sequence in MAIN:
            file.writelines(sequence.id + ',')
        file.writelines('\n')

        #Assumes target is first sequence
        target = MAIN[0].id
        x = calc_coverage()
        for item in x:
            file.writelines(str(item) + ',')

        count += 1
    with open("average_id.csv","w") as file:
        for item in range(len(MAIN[0]) - 1):
            file.writelines(str(item) + ',')
        file.writelines(str(len(MAIN[0])-1)) #Avoid comma in last step
        file.writelines("\n") #Finish enumeration
        x = calc_average_id()
        for item in x:
            file.writelines(str(item) + ',')
        count +=1 

    print("Count written: ",count)


def run():
    generate_stats()
    generate_csv()





