def read_fasta(filename):
    ## Reading in the Antibiotic Resistence (AR) File 

    oldheader = None
    file = open(filename, 'r')
    for line in file:
        line = line.strip()
        if line.startswith(">"):
            newheader = line
            if oldheader is not None:
                yield dna, oldheader
            dna = ""
            oldheader = newheader
        else:
            dna += line
    yield dna, oldheader


for dna, header in read_fasta("ARsmall.txt"):
    print(header,dna)
    