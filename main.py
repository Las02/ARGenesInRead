import sys
import re
import gzip
def read_fasta(filename):
    ## Reading in the Antibiotic Resistence (AR) File 
    AR_file = open(filename, 'r')
    line = 'void'
    while line != '' and not line.startswith('>'):
        line = AR_file.readline()

    while line != '':
        line = line.strip()
        header = line
        dna = ''
        line = AR_file.readline()
        while line != "" and line[0] != '>':

            dna += line
            line = AR_file.readline()
        yield dna, header
    AR_file.close()
    
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

def FindKmer(dna, kmer_len):
    '''Find Kmers from dna string and return list with them'''

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
        range_list.append([i, len(dna)])

    return kmer_list, range_list

def CountEachKmer(kmer, gene_data, found_kmer_in_gene,kmer_length):
    '''if input kmer is equal to kmer in gene_data 
       make array with depht and store in found_kmer_in_gene'''

    def AddDepth(kmer_pos, found_kmer_in_gene, gene,kmer_length):
        """Add the count for the kmer to position from which it was found"""
        from_range = kmer_pos
        to_range = kmer_pos + kmer_length     
        for i in range(from_range, to_range):
            found_kmer_in_gene[gene][i] = 1

    # If the kmer is equal to a kmer in the genes
    if kmer in gene_data:
        # Go through each gene which has the kmer and add depht to it.
        for genename in gene_data[kmer]:
            len_gene = gene_data[kmer][genename][1]
            kmer_pos = gene_data[kmer][genename][0]

            if genename in found_kmer_in_gene:
                # Add depht to it, corresponding the kmer found
                AddDepth(kmer_pos, found_kmer_in_gene, genename, kmer_length)
            else:
                # Make vector of [0] to represent depht of each nt corresponding to length of gene
                found_kmer_in_gene[genename] = [0] * len_gene

                # Add depht to it, corresponding the kmer found
                AddDepth(kmer_pos, found_kmer_in_gene, genename ,kmer_length)
        


def AddToDatastructure(gene_data, kmer_list, header, range_list):
    """Adds the kmer_list and header to the gene_data """
    for i in range(len(kmer_list)):
        kmer = kmer_list[i]
        kmer_range = range_list[i]

        # If the kmer is allready assigned to an AR gene, add the additional
        if kmer in gene_data:
            gene_data[kmer][header] = (kmer_range)
        # Else add the kmer to the datastructure
        else:
            gene_data[kmer] = {header:kmer_range}

def find_coverage_depht(dna):

    count = 0
    total_depht = 0

    for base in dna:
        total_depht += base

        if base != 0:
            count += 1
        
    coverage = count / len(dna)
    avg_depht = total_depht / len(dna)
    if coverage > 0.95 and avg_depht > 0.5:
        return coverage, avg_depht
    else:
        return None, None
def judge(gene, read):

    max_space = 1
    side_bonus = int(0.5 * len(read))
    threshold_score = int(0.9 *len(read))

    maxcount = 0
    count = 0
    exitcount = 0
    in_kmer = False
    len_dna = len(gene)

    for pos,value in enumerate(gene):
        if value == 1: in_kmer = True

        if in_kmer:
            if value == 0:
                exitcount += 1
            else:
                count += 1
                # If at either side position, add a bonus to the score
                if pos in [0, len_dna-1]:
                    count += side_bonus

            if count > maxcount:
                maxcount = count
            
            if exitcount >= max_space:
                count = 0
                exitcount = 0
                in_kmer = False
    
    return maxcount >= threshold_score   




# Set the kmer lenght to look for
kmer_length = 4

# Stores {Kmer: {"AR":kmer_range}}, The Kmer as key. 
#then there is a inner dict with the AR gene as key and the kmer_ range as values
gene_data = dict()
Could_be_kmers = dict()

#filename = "Unknown3_raw_reads_1.txt.gz"
#read_filename = "test.read.txt.gz"
gene_filename = 'ARsmall.txt'
read_filename = "smallfastaseq.txt.gz"

for dna, header in read_fasta(gene_filename):
    # Finding all the posible kmers and its positions
    (kmer_list, range_list) = FindKmer(dna, kmer_length)

    # Adding the kmer_list and header to the gene_data 
    AddToDatastructure(gene_data, kmer_list, header, range_list)





## Reading in the sequenceing file

gene_count = dict()

for dna in read_qfasta(read_filename):
    # Find all posible kmers
    (kmer_list, range_list) = FindKmer(dna, kmer_length)
    # Count the found posible kmers
    found_kmer_in_gene = dict()
    for kmer in kmer_list:
        CountEachKmer(kmer, gene_data, found_kmer_in_gene,kmer_length)

    # Do some check with it
    for genename,gene in found_kmer_in_gene.items():
        if judge(gene, dna):
            # Add den til final count
            if genename in gene_count:
                for i in range(len(gene_count[genename])):
                    gene_count[genename][i] += gene[i]    # Skal plusse element pr element  
            else :
                gene_count[genename] = gene 
    print(found_kmer_in_gene)

coverage_depht = dict()
for genename, dna in gene_count.items():
    (coverage, avg_depht) = find_coverage_depht(dna)
    if coverage != None:
        coverage_depht[genename] = (coverage, avg_depht)             

sorted_coverage_depht = sorted(coverage_depht, key= coverage_depht.get, reverse = True)


for genename in sorted_coverage_depht:
    coverage = coverage_depht[genename][0]
    avg_depht = coverage_depht[genename][1]
    (name,AR) = genename.split(maxsplit = 1)
    print("Name:", name,"ARgene:", AR ,"coverage:",  coverage, "avg_depght: ",avg_depht)


