ResistanceGenes = open('resistance_genes.fsa.txt', 'r')

import sys
import re


line = 'void'


while line != '' and line[0] != '>':
    line = ResistanceGenes.readline()
while line != '':
    line = line.strip()
    print(line)
    dna = ''
    line = ResistanceGenes.readline()
    while line != "" and line[0] != '>':
        line = re.sub("\s", "", line).upper()
        if re.search("[^ATGC]", line) is not None:
            print("Invalid sequence line: " + line)
            sys.exit(1)
        dna += line
        line = ResistanceGenes.readline()
        print(dna)
