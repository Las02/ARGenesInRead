def read_fasta(filename):
    ## Reading in the Antibiotic Resistence (AR) File 
    AR_file = open(filename, 'r')
    dna =""

    for line in AR_file:
        if line[0] == ">":
            header = line
            dna = ""
            print(header)
        else:
            dna += line[:-1]
            if line == header:
                break
        print(dna)



    return dna
print(read_fasta("ARsmall.txt"))