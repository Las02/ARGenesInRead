
"""
real    1m52.000s
user    1m46.734s
sys     0m2.438s

"""

import sys
import re
import gzip

def FindKmer(dna, kmer_len):
    '''Find Kmers from dna string and return list with them'''

    # Return None if there is no possible kmer's
    if kmer_len > len(dna):
        return

    kmer_list = []
    from_range = 0
    to_range = len(dna) - kmer_len + 1   #+1 to also add last kmer
    for i in range(from_range, to_range):

        kmer = dna[i: i + kmer_len]
        kmer_list.append(kmer)

    return kmer_list

def CountEachKmer(kmer_list, Data_Structure):
    '''Counts the amount of each kmer in the dict: kmer_count'''
    for kmer in kmer_list:
        if kmer in Data_Structure:
            Data_Structure[kmer]["count"] += 1

def CountKmerPerAR(nKmer_per_AR, Data_Structure):
    '''Counts the amount of all kmers for each AntibioticResistence Gene in the dict: nKmer_per_AR'''

    # Goes trough each kmer
    for data in Data_Structure.values():

        # Go through each AR_gene for which the kmer was found and add the count to it
        for AR_gene in data["AR_genes"]:
            
            if data["count"] == 0:
                pass

            elif AR_gene in nKmer_per_AR:
                nKmer_per_AR[AR_gene] += data["count"]
            
            else:
                nKmer_per_AR[AR_gene] = data["count"]


def AddToDatastructure(Data_Structure, kmer_list, header):
    """Adds the kmer_list and header to the Data_structure """
    for kmer in kmer_list:
        # If the kmer is allready assigned to an AR gene, add the additional
        if kmer in Data_Structure:
            Data_Structure[kmer]["AR_genes"].add(header)
        # Else add the kmer to the datastructure
        else:
            Data_Structure[kmer] = {"count":0, "AR_genes":{header}}

                

# Set the kmer lenght to look for
kmer_length = 5
# Stores {Kmer: {"count":0,{"AR":}}}, The Kmer as key, and then both the number of time it is seen, and the AR genes which has the kmer
Data_Structure = dict()
Could_be_kmers = dict()
nKmer_per_AR = dict()

## Reading in the Antibiotic Resistence (AR) File 

AR_file = open('resistance_genes.fsa.txt', 'r')
filename = "Unknown3_raw_reads_1.txt.gz"
#AR_file = open('ARsmall.txt', 'r')
#filename = "smallfastaseq.txt.gz"

line = 'void'
while line != '' and line[0] != '>':
    line = AR_file.readline()

while line != '':
    line = line.strip()
    header = line
    dna = ''
    line = AR_file.readline()
    while line != "" and line[0] != '>':

        # Quit the program if the found DNA is not DNA
        # TODO Failer hvis DNA string starter med unkown eg. "N"
        line = re.sub("\s", "", line).upper()
        if re.search("[^ATGC]", line) is not None:
            print("Invalid sequence line: " + line)
            sys.exit(1)

        dna += line
        line = AR_file.readline()
    
    # Finding all the posible kmers
    kmer_list = FindKmer(dna, kmer_length)

    # Adding the kmer_list and header to the Data_structure 
    AddToDatastructure(Data_Structure, kmer_list, header)


## Reading in the sequenceing file


sample_file = gzip.open(filename, "r")

last_line = ""
for line in sample_file:
    line = line.decode("utf-8") 

    # Extracting only the dna
    # TODO hver 110% sikre på at den ikke kan være på 2 linjer
    if last_line.startswith("@"):
        dna = line.strip()

        # Find all posible kmers
        kmer_list = FindKmer(dna, kmer_length)
        # Count the foind posible kmers
        CountEachKmer(kmer_list, Data_Structure)
        
    last_line = line
    

CountKmerPerAR(nKmer_per_AR, Data_Structure)

print(nKmer_per_AR)


AR_file.close()
sample_file.close()

