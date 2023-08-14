from functions import argument_parser, read_fasta, read_fastq, get_kmers_and_pos, read_is_valid, coverage_stats
import sys
import gzip

# Setting up the argument parser
(kmer_length,gene_filename,read_filename) = argument_parser(sys.argv, 19, "resistance_genes.fsa", "Unknown3_raw_reads_1.txt.gz")


# Read in the file with the antibiotic resistence genes and save it in the datastructure kmer2gene2kmerpos
kmer2gene2kmerpos = dict()
trans = str.maketrans("ACTG","TGAC")
for non_complement_dna, header in read_fasta(gene_filename):
    # Complement the dna strand (non_complement_dna) given in the gene files
    complement_dna= non_complement_dna.translate(trans)

    # Go through all possible variations of the gene dna
    for dna in [non_complement_dna, complement_dna,non_complement_dna[::-1], complement_dna[::-1]]:
        # Make "gene_data" datastructure which consists of: {Kmer: {"gene_name":(kmer_position_in_gene,length_gene)}}
        (kmer_list, kmer_positions) = get_kmers_and_pos(dna, kmer_length)
        for kmer, kmer_pos_in_gene in zip(kmer_list, kmer_positions):
            if kmer in kmer2gene2kmerpos:
                kmer2gene2kmerpos[kmer][header] = (kmer_pos_in_gene, len(dna))
            else:
                kmer2gene2kmerpos[kmer] = {header : (kmer_pos_in_gene, len(dna))}


# Read in the fastaq file, for each read evaluate if read is valid.
# If valid, add it to total depht count for each gene saved in TOTALgene2depht_count
TOTALgene2depht_count = dict()
for dna_read in read_fastq(read_filename):
    (kmer_list, _) = get_kmers_and_pos(dna_read, kmer_length)

    # For all kmers check if they match a kmer in the gene
    # If they match save where they covered the gene in the datastructure: "READgene2depht_count"
    # where the depht_count is a list where 1 is matching position and 0 is non matching
    READgene2depht_count = dict()
    for kmer in kmer_list:
        # If the kmer is equal to a kmer in the genes
        if kmer in kmer2gene2kmerpos:
            # Go through each gene which has the kmer and add depht to it.
            for genename in kmer2gene2kmerpos[kmer]:
                len_gene = kmer2gene2kmerpos[kmer][genename][1]
                kmer_pos = kmer2gene2kmerpos[kmer][genename][0]
                if genename not in READgene2depht_count:
                    # Make vector of [0] to represent depht of each nt corresponding to length of gene
                    READgene2depht_count[genename] = [0] * len_gene
                # Add depht to it, corresponding the kmer found
                for i in range(kmer_pos, kmer_pos + kmer_length):
                    READgene2depht_count[genename][i] = 1
    
    # If the read is valid, add the gene to the final depht_count for each gene (FINALgene2depht_count)
    for genename, depht_count in READgene2depht_count.items():
        if read_is_valid(depht_count, dna_read):
            if genename in TOTALgene2depht_count:
                # elementwise addition
                for i in range(len(TOTALgene2depht_count[genename])):
                    TOTALgene2depht_count[genename][i] += depht_count[i]    
            else :
                TOTALgene2depht_count[genename] = depht_count


# Calculate coverage, min depht and avg coverage based on depht_counts from TOTALgene2depht_count
# Add genes with coverage > 95 percent and min_depht >10 gene2coverage_depht
gene2coverage_depht = dict()
for genename, dna in TOTALgene2depht_count.items():
    (coverage, avg_depht, min_depht) = coverage_stats(dna)
    if coverage > 0.95 and min_depht > 10:
        gene2coverage_depht[genename] = (coverage, avg_depht)

# sort the genes based on coverage then depht
sorted_gene_coverage_depht = sorted(
            gene2coverage_depht, 
            key= gene2coverage_depht.get, 
            reverse = True)

# output and format sorted_gene_coverage_depht to tab seperated format
print("gene\tresistence\tcoverage\tavg_depht")
for genename in sorted_gene_coverage_depht:
    coverage = round(gene2coverage_depht[genename][0],2)
    avg_depht = round(gene2coverage_depht[genename][1],2)
    (gene,resistence) = genename.split(maxsplit = 1)
    gene = gene[1:]
    print(f"{gene}\t{resistence}\t{coverage}\t{avg_depht}")
