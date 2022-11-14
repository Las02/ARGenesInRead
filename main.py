import sys
import re
import gzip
import argparse

#gene_filename = 'ARsmall.txt'
#read_filename = "smallfastaseq.txt.gz"
#gene_filename = 'resistance_genes.fsa'
#read_filename = "Unknown3_raw_reads_1.txt.gz"
#read_filename ="mini.txt.gz"


def read_fasta(filename):
    '''Reading in several fasta files'''

    oldheader = None
    file = open(filename, 'r')
    for line in file:
        line = line.strip()
        if line.startswith(">"):    #if line is header
            newheader = line
            # if not the first line
            if oldheader is not None:
                yield dna, oldheader
            dna = ""
            oldheader = newheader
        else:
            dna += line
    # To yield the last fastaformatted dna
    yield dna, oldheader
    file.close()
    
def read_qfasta(filename):
    '''Extract the dna from a gzippedfastaQ file '''

    sample_file = gzip.open(filename, "r")

    last_line = ""
    for line in sample_file:
        line = line.decode("utf-8")

        # Extracting the dna
        if last_line.startswith("@"):
            dna = line.strip()
            yield dna

        last_line = line

    sample_file.close()

def find_all_kmers(dna, kmer_len):
    '''Find Kmers from dna string and return list with them
       in additon to list with their positions'''

    # Return None if there is no possible kmer's
    if kmer_len > len(dna):
        return

    kmer_list = []
    range_list = []
    from_range = 0
    to_range = len(dna) - kmer_len + 1   #+1 to also add last kmer
    for i in range(from_range, to_range):

        kmer = dna[i: i + kmer_len]
        kmer_list.append(kmer)
        range_list.append(i)

    return kmer_list, range_list

def read_is_valid(gene, read):
    '''
    return True if the read is covering enough of the gene,
    else returns False
    '''
    # Various parameters
    max_space = 1   # The maximum space in the read covering the gene
    side_bonus = int(0.5 * len(read))   # How much of the read can be outside the gene
    threshold_score = int(0.5 *len(read))   #The percent of the read which needs to cover the gene

    # Various init
    maxcount = 0
    count = 0
    exitcount = 0
    in_kmer = False
    len_dna = len(gene)

    # Score the reads coverage of the gene
    for pos,value in enumerate(gene):

        # Start counting the coverage at first [1] in gene
        if value == 1: in_kmer = True
        if in_kmer:
            if value == 0:
                exitcount += 1
            else:
                count += 1
                # If at either side position, add a bonus to the score
                if pos in [0, len_dna-1]:
                    count += side_bonus

            # If a better fitting coverage of the read to the gene is found, use it
            if count > maxcount:
                maxcount = count
            
            # Stop counting the coverage if there is too much spacing
            if exitcount >= max_space:
                count = 0
                exitcount = 0
                in_kmer = False

    return maxcount >= threshold_score

def coverage_stats(dna):
    '''from depht array return coverage, avg depht and min_depht'''

    # Init
    count = 0
    total_depht = 0
    min_depth = 99999

    # Find total_depht, coverage and avg_depht
    for base in dna:
        total_depht += base
        if base != 0:
            count += 1
        if base < min_depth:
            min_depth = base
    coverage = count / len(dna)
    avg_depht = total_depht / len(dna)
    return coverage, avg_depht, min_depth

### Setting up the parser for the program
parser = argparse.ArgumentParser(
    prog = "coverage_counter",
    description =  """
    Compares genes given in fasta format to reads in fastaQ format
    returns the genes with best coverage and high enough depht
    """
    )
# Adding arguments to the parser and defining them in the program
parser.add_argument("-g", dest ="gene_filename", type=str)
parser.add_argument("-r", dest ="read_filename", type=str)
parser.add_argument("-k", dest = "kmer_length", type=int, default = 19)
args = parser.parse_args()
read_filename = args.read_filename
kmer_length = args.kmer_length

### Reading in the genefile and saving it in the datastructure: gene_data
kmer2gene2kmerpos = dict()
for dna, header in read_fasta(args.gene_filename):
    (kmer_list, kmer_positions) = find_all_kmers(dna, args.kmer_length)
    # Make gene_data datastructure which consists of: {Kmer: {"gene_name":(kmer_position_in_gene,length_gene)}}
    for kmer, kmer_pos_in_gene in zip(kmer_list, kmer_positions):
        if kmer in kmer2gene2kmerpos:
            kmer2gene2kmerpos[kmer][header] = (kmer_pos_in_gene, len(dna))
        else:
            kmer2gene2kmerpos[kmer] = {header : (kmer_pos_in_gene, len(dna))}

### Reading in the fastaq file and ----
Could_be_kmers = dict()
gene_count = dict()
for dna_read in read_qfasta(read_filename):
    (kmer_list, _) = find_all_kmers(dna_read, args.kmer_length)

    # Save the depht for all found kmers in the read in "gene2depht_count"
    gene2depht_count = dict()
    for kmer in kmer_list:
        # If the kmer is equal to a kmer in the genes
        if kmer in kmer2gene2kmerpos:
            # Go through each gene which has the kmer and add depht to it.
            for genename in kmer2gene2kmerpos[kmer]:
                len_gene = kmer2gene2kmerpos[kmer][genename][1]
                kmer_pos = kmer2gene2kmerpos[kmer][genename][0]

                if genename not in gene2depht_count:
                    # Make vector of [0] to represent depht of each nt corresponding to length of gene
                    gene2depht_count[genename] = [0] * len_gene

                # Add depht to it, corresponding the kmer found
                for i in range(kmer_pos, kmer_pos + kmer_length):
                    gene2depht_count[genename][i] = 1

    # If the read is valid, add it to the depht count
    for genename, gene in gene2depht_count.items():
        if read_is_valid(gene, dna):
            # Add den til final count
            if genename in gene_count:
                # elementwise addition
                for i in range(len(gene_count[genename])):
                    gene_count[genename][i] += gene[i]     
            else :
                gene_count[genename] = gene



###################################################################
print(gene_count)   # adder ikke kestra på
# Find the actual gene coverage
coverage_depht = dict()
for genename, dna in gene_count.items():
    (coverage, avg_depht, min_depht) = coverage_stats(dna)
    if coverage > 0.95 and avg_depht > 10:
        coverage_depht[genename] = (coverage, avg_depht)
        #hvis der ikke er nogen

sorted_coverage_depht = sorted(coverage_depht, key= coverage_depht.get, reverse = True)

for genename in sorted_coverage_depht:
    coverage = coverage_depht[genename][0]
    avg_depht = coverage_depht[genename][1]
    (name,AR) = genename.split(maxsplit = 1)
    print("Gene name:", name[1:])
    print("Antibiotic restistance gene:", AR)
    print("Coverage:",  coverage*100,"%")
    print("Average depht:",avg_depht)
