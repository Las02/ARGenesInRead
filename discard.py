def read_fasta(filename):
    ## Reading in the Antibiotic Resistence (AR) File 
    AR_file = open(filename, 'r')
    dna =""

    line = AR_file.readline()
    print(line)
    while line != "":
        if line[0] == (">"):
            header = line
        else: 
            dna += line
        line = AR_file.readline()
    yield dna, header
    AR_file.close()

for dna, header in read_fasta("ARsmall.txt"):
    print(header,dna.strip())
    