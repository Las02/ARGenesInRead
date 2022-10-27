import gzip

filename = "Unknown3_raw_reads_1.txt.gz"
sample_file = gzip.open(filename, "r")

last_line = ""
for line in sample_file:
    line = line.decode("utf-8") 

    if last_line.startswith("@"):
        dna = line.strip()
        print(line)

    last_line = line
 