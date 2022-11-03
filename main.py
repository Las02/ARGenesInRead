import sys
import re
import gzip


def judge(dna):

    max_space = 2
    side_bonus = 1
    threshold_score = 2

    maxcount = 0
    count = 0
    exitcount = 0
    in_kmer = False
    len_dna = len(dna)

    for pos,value in enumerate(dna):
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

def CountEachKmer(kmer_list, Data_Structure, nKmer_per_AR):
    '''Counts the amount of each kmer in the dict: kmer_count'''

    # Goes through each found kmer
    for kmer in kmer_list:
        # If the kmer is equal to an AR gene
        if kmer in Data_Structure:
            AR_to_range = Data_Structure[kmer]
    
            # Go through each AR_gene for which the kmer was found and add the count to it
            for AR_gene in AR_to_range:
                
                if AR_gene in nKmer_per_AR:
                    from_to_len = AR_to_range[AR_gene]
                    AddDepth(from_to_len, nKmer_per_AR, AR_gene)
                
                else:
                    # Extract [from, to, len(dna)]
                    from_to_len = AR_to_range[AR_gene]
                    # Make vector of [0] to represent depht of each nt
                    length_of_gene = from_to_len[1]
                    nKmer_per_AR[AR_gene] = [0] * length_of_gene

                    # Add the count for the kmer to the specific place
                    AddDepth(from_to_len, nKmer_per_AR, AR_gene)
            

def AddDepth(from_to_len, nKmer_per_AR, AR_gene):
    """Add the count for the kmer to the specific place"""
    from_range = from_to_len[0]
    to_range = from_range + kmer_length
    add = 1
    for i in range(from_range, to_range):
        nKmer_per_AR[AR_gene][i] = 1


def AddToDatastructure(Data_Structure, kmer_list, header, range_list):
    """Adds the kmer_list and header to the Data_structure """
    for i in range(len(kmer_list)):
        kmer = kmer_list[i]
        kmer_range = range_list[i]

        # If the kmer is allready assigned to an AR gene, add the additional
        if kmer in Data_Structure:
            Data_Structure[kmer][header] = (kmer_range)
        # Else add the kmer to the datastructure
        else:
            Data_Structure[kmer] = {header:kmer_range}
                

# Set the kmer lenght to look for
kmer_length = 5

# Stores {Kmer: {"AR":kmer_range}}, The Kmer as key. 
#then there is a inner dict with the AR gene as key and the kmer_ range as values
Data_Structure = dict()
Could_be_kmers = dict()


## Reading in the Antibiotic Resistence (AR) File 

#AR_file = open('resistance_genes.fsa.txt', 'r')
#filename = "Unknown3_raw_reads_1.txt.gz"
AR_file = open('ARsmall.txt', 'r')
filename = "smallfastaseq.txt.gz"

line = 'void'
while line != '' and line[0] != '>':
    line = AR_file.readline()

while line != '':
    line = line.strip()
    header = line
    dna = ''
    line = AR_file.readline()
    while line != "" and line[0] != '>':

        # Quit the program if the found DNA is not DNA
        # TODO Failer hvis DNA string starter med unkown eg. "N"
        line = re.sub("\s", "", line).upper()
        if re.search("[^ATGC]", line) is not None:
            print("Invalid sequence line: " + line)
            sys.exit(1)

        dna += line
        line = AR_file.readline()
    
    # Finding all the posible kmers and its positions
    (kmer_list, range_list) = FindKmer(dna, kmer_length)

    # Adding the kmer_list and header to the Data_structure 
    AddToDatastructure(Data_Structure, kmer_list, header, range_list)


## Reading in the sequenceing file


sample_file = gzip.open(filename, "r")

last_line = ""
for line in sample_file:
    line = line.decode("utf-8") 

    # Extracting only the dna
    # TODO hver 110% sikre på at den ikke kan være på 2 linjer
    if last_line.startswith("@"):
        dna = line.strip()

        # Find all posible kmers
        (kmer_list, range_list) = FindKmer(dna, kmer_length)
        # Count the found posible kmers
        nKmer_per_AR = dict()
        CountEachKmer(kmer_list, Data_Structure, nKmer_per_AR)
        
        gene_count = dict()

        # Do some check with it
        for genename,gene in nKmer_per_AR.items():
            if judge(gene):
                # Add den til final count
                if genename in gene_count:
                    for i in range(len(gene_count[genename])):
                        gene_count[genename][i] += gene[i]    # Skal plusse element pr element  
                else :
                    gene_count[genename] = gene 
                
    last_line = line


"""

# Quick print the found values
# And the depht of each pp
for item in gene_count.items():
    count = sum(item[1])
    print("*"*50)
    print("AR genename:", item[0], end="\t")
    print("coverage:", count, "out of", len(item[1]), "positons")
    #print("The depht of each position:")
    #print(item[1])
    print("*"*50)
print("kmers of size:", kmer_length)
#Count found to be: 129.0, 118.0


"""


AR_file.close()
sample_file.close()

