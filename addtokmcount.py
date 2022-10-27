def Count(list,dict1):
    for Kmer in KmerList:
        if Kmer in KmCount:
            KmCount[Kmer] += 1
        else:
            KmCount[Kmer] = 1
