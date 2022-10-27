
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
for header in ["AR1", "AR2"]:
    Count(Kmerlist,KmCount,header,ARCount)



finalcount = dict()

def FinalCounter(finalcount, ARCount, KmCount):
    for AR in ARCount.items():
        kmerlist = AR[1]
        count = 0
        for kmer in kmerlist:
            count += KmCount[kmer]
        
        finalcount[AR[0]] = count


FinalCounter(finalcount, ARCount, KmCount)
print(finalcount)