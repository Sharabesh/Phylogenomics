from main import *

THIRTY_PERCENT_FILE = ""

# Gather Homologs of various relationships to KNCA1_HUMAN
def gather_homologs():
    sys_gather_homologs()
    parse_xml()
    generate_msa()


#Using Belvu select only those with specified percent identity