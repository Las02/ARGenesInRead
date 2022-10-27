

def FindKmer(dna, kmer_len = 3):
    '''Find Kmers from dna string and return list with them'''

    # Return None if there is no kmer's
    if kmer_len > len(dna):
        return

    kmer_list = []
    from_range = 0
    to_range = len(dna) - kmer_len + 1   #+1 to also add last kmer
    for i in range(from_range, to_range):

        kmer = dna[i: i + kmer_len]
        kmer_list.append(kmer)

    return kmer_list

print(FindKmer("ATGAT", 3))