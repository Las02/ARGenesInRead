

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



read = "ATCG"

read_length = len(read)
readcount = [0] * read_length

kmerlist = FindKmer(read,2)[0]

print(kmerlist)
