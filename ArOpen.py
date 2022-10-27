ResistanceGenes = open('resistance_genes.fsa.txt', 'r')

import sys
import re

from findkmer import FindKmer




def FinalCounter(finalcount, ARCount, KmCount):
    for AR in ARCount.items():
        kmerlist = AR[1]
        count = 0
        for kmer in kmerlist:
            count += KmCount[kmer]
        
        finalcount[AR[0]] = count


    
def Count(list,KmCount,header,ARCount):
    for Kmer in list:

        # Den her skal over ift til den anden fil
        ARCount[header] = []
        if Kmer in KmCount:
            KmCount[Kmer] += 1
        else:
            KmCount[Kmer] = 1
        
        # Men ikke den her del som skal blive
    ARCount[header] = list




finalcount = dict()
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

#print(ARCount)
FinalCounter(finalcount, ARCount, KmCount)
print(finalcount)