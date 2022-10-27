ResistanceGenes = open('resistance_genes.fsa.txt', 'r')

import sys
import re

from findkmer import FindKmer




def FinalCounter(finalcount, ARCount, KmCount):
    for AR in ARCount.items():
        kmerlist = AR[1]
        count = 0
        for kmer in kmerlist:
            if not kmer in KmCount:
                print(kmer, "is not in Kmcount")
            else:
                count += KmCount[kmer]
        
        finalcount[AR[0]] = count


    
def Count(list,KmCount):
    for Kmer in list:
        if Kmer in KmCount:
            KmCount[Kmer] += 1
        else:
            KmCount[Kmer] = 1

def addToARcount(ARCount, list, header):
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
    dna = ''
    line = ResistanceGenes.readline()
    while line != "" and line[0] != '>':

        line = re.sub("\s", "", line).upper()
        if re.search("[^ATGC]", line) is not None:
            print("Invalid sequence line: " + line)
            sys.exit(1)

        dna += line
        line = ResistanceGenes.readline()
    Kmer_list = FindKmer(dna,19)
    addToARcount(ARCount, Kmer_list, header)

#########################################
#FASTAQ

import gzip

filename = "Unknown3_raw_reads_1.txt.gz"
sample_file = gzip.open(filename, "r")

last_line = ""
i=0
for line in sample_file:
    line = line.decode("utf-8") 

    if last_line.startswith("@"):
        dna = line.strip()

        Kmer_list = FindKmer(dna,19)
        Count(Kmer_list, KmCount)

    # Break early for testing
    i+=1
    if i == 10:
        break

    last_line = line


FinalCounter(finalcount, ARCount, KmCount)

















# ARCount is now done





#FinalCounter(finalcount, ARCount, KmCount)
#print(finalcount)
