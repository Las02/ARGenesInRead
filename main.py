# XXX Det virker nu men...
# XXX Lige nu kan CountEachKmer blive ekstremt stor hvis input FASTAq er stor
# XXX Vi burde kun adde til den hvis vi har Kmeren i vores dict
# XXX {Kmer: [count,{AR}]} <= evt. :(

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

def CountEachKmer(kmer_list, kmer_count):
    '''Counts the amount of each kmer in the dict: kmer_count'''
    for kmer in kmer_list:
        if kmer in kmer_count:
            kmer_count[kmer] += 1
        else:
            kmer_count[kmer] = 1

def CountKmerPerAR(nKmer_per_AR, AR_to_kmers, kmer_count):
    '''Counts the amount of all kmers for each AntibioticResistence Gene in the dict: nKmer_per_AR'''
    for AR_gene in AR_to_kmers.items():
        kmerlist = AR_gene[1]
        count = 0
        for kmer in kmerlist:
            
            # If the kmer has been found (meaning its in kmer_count)
            # then add the amount of times it has been found to the total count for the AR_gene
            if kmer in kmer_count:
               count += kmer_count[kmer]

        AR_genename = AR_gene[0]
        nKmer_per_AR[AR_genename] = count


# Set the kmer lenght to look for
kmer_length = 6
# Store which Antibiotic Resistence (AR) gene has which kmers
AR_to_kmers = dict()
Could_be_kmers = dict()

## Reading in the Antibiotic Resistence (AR) File 

AR_file = open('ARsmall.txt', 'r')

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
    # Adding all posible kmers to the relevant AR gene
    AR_to_kmers[header] = set(kmer_list)  

# Dict used to store each found kmer, 
kmer_count = dict() 
nKmer_per_AR = dict()

## Reading in the sequenceing file

filename = "smallfastaseq.txt.gz"
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
        CountEachKmer(kmer_list, kmer_count)
        
    last_line = line


CountKmerPerAR(nKmer_per_AR, AR_to_kmers, kmer_count)
print(nKmer_per_AR)



AR_file.close()
sample_file.close()

