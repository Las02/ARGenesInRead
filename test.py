
from findkmer import FindKmer

def Count(Kmerlist,KmCount,header,ARCount):
    for Kmer in Kmerlist:
        
        if Kmer in KmCount:
            KmCount[Kmer] += 1
        else:
            KmCount[Kmer] = 1
            
    ARCount[header] = Kmerlist

Kmerlist = FindKmer("ATCGY", 3)
KmCount = dict()
ARCount = dict()
header = "AR"
Count(Kmerlist,KmCount,header,ARCount)


#print(Kmerlist)