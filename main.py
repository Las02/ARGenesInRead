import gzip

def CompareDNA(dna):
    print(dna)

filename = "Unknown3_raw_reads_1.txt.gz"
sample_file = gzip.open(filename, "r")

in_dna = False
for line in sample_file:
    line = line.decode("utf-8")

    if in_dna:
        dna = line.strip()
        CompareDNA(dna)
        in_dna = False

    if line.startswith("@"):
        in_dna = True

    
