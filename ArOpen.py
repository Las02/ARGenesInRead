ResistanceGenes = open('resistance_genes.fsa.txt', 'r')

import sys
import re

from findkmer import FindKmer

def Count(Kmerlist,KmCount,header,ARCount):
    for Kmer in Kmerlist:

        if Kmer in KmCount:
            KmCount[Kmer] += 1
        else:
            KmCount[Kmer] = 1

    ARCount[header] = Kmerlist
    
KmerList = list()
KmCount = dict()
ARCount = dict()

#Reading the resistance genes file
line = 'void'

while line != '' and line[0] != '>':
    line = ResistanceGenes.readline()
while line != '':
    line = line.strip()
    header = line
    #print(line)
    dna = ''
    line = ResistanceGenes.readline()
    while line != "" and line[0] != '>':
        line = re.sub("\s", "", line).upper()
        if re.search("[^ATGC]", line) is not None:
            print("Invalid sequence line: " + line)
            sys.exit(1)
        dna += line
        line = ResistanceGenes.readline()
    KmerList = FindKmer(dna,19)
    Count(KmerList,KmCount,header,ARCount)

print(ARCount)
