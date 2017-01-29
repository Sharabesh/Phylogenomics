import Bio
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML
from Bio.Blast.NCBIWWW import qblast
from Bio.Align.Applications import MafftCommandline
from Bio.Align import MultipleSeqAlignment
import os


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

E_VAL_THRESH = .05
ALIGN_PERCENT_THRESH = 80
DATABASE = "nr.61"
CONSENSUS_THRESH = .5
INPUT_SEQUENCE = "operating_reqs/fasta.txt"
save_file_name = "operating_reqs/alignment.xml"
recs_file = "operating_reqs/records.fasta"


def length(fasta):
    input = SeqIO.read(fasta,format="fasta")
    return len(input)


input_len = length(INPUT_SEQUENCE)

#
# """Uses BioPython rather than System os"""
# def gather_homologs(fasta): #File location passed as input ("fasta.txt")
#     blasted = NcbiblastpCommandline(query=fasta,db="nr.61",evalue=1,outfmt=5,out="alignment.xml")
#     blasted()


"""Using System os
-qcov_hsp_perc allows for a defined query coverage """
def sys_gather_homologs(fasta):
    generic = "blastp -db {0} -query {1} -evalue {2}" \
              " -outfmt {3} -out {4} -qcov_hsp_perc {5}"
    db=DATABASE
    query = fasta
    evalue = E_VAL_THRESH #Change if wanted
    outfmt = 5 # Indicates BLAST XML to write to file
    out = "alignment.xml" #Data written to alignment.xml
    query_coverage = ALIGN_PERCENT_THRESH #Minimum required overlap
    formatted = generic.format(db,query,evalue,outfmt,out,query_coverage)
    os.system(formatted)



def parse_xml(): #This is alignment.xml
    with open(save_file_name) as result_handle:
        blast_record = NCBIXML.read(result_handle)
    recs = []
    i = 0
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            filtered_seq = Seq(str(hsp.sbjct), IUPAC.protein)
            filtered_id = alignment.title.split("|")[3]
            filtered_name = alignment.title
            filtered_seqrec = SeqRecord(filtered_seq, id=filtered_id, name=filtered_name)
            recs.append(filtered_seqrec)
    _ = SeqIO.write(recs, recs_file, "fasta")
    print("Done")

def generate_msa():
    mafft_cline = MafftCommandline(input=recs_file)
    stdout, stderr = mafft_cline()
    with open(recs_file, "w") as handle:
        handle.write(stdout)

"""Note specific columns can be removed using the argument 
    -- select {{n,l,m-k}}
"""
def mask_msa():
    terminal_source = "trimal -automated1 -in {0} -out {1}"
    os.system(terminal_source.format(recs_file,recs_file))



"""Uses RaxML to generate tree which is read out by Bio's Phylo module"""
def generate_tree():
    pass