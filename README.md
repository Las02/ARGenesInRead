# Goal
The goal of this project has been to develop a commandline program to find resistance genes in reads from metagenomic samples.  
The program runs without any dependencies apart from base Python in O(N) time corresponding to the size of the metagenomic sample.  

# How does it work?
The program makes a database structure which store all k-mers based on a FASTA file of antibiotic resistance genes.
It then looks at each k-mer from each read in an input FASTAQ file to see if they match.
In case of a match the read might be from bacterial dna which has a resistance gene.
The output of the program is then antibiotic genes with coverage above 95% and a depht above 10.  

# Program manual
The program is named "ARGenesInRead.py," It takes 3 arguments:  
```
-g : a fasta file containg the resistence genes
-r : a gzipped FASTAQ file containg read data
-k : an argument giving the kmer size to look for
```
The following is an example of using the program:
```
python ./ARGenesInRead.py -g resistance_genes.fsa -r Unknown3_raw_reads_1.txt.gz -k 19
```
The output of the program is printed, and can be used as if it was std.out,
It is in the format of a tab separated file. Here is an example:

```
    gene    resistence      coverage        avg_depht
fosA_3_ACWO01000079     Fosfomycin resistance:  1.0     47.15     
aac(6')Ib-cr_1_DQ303918 Fluoroquinolone and aminoglycoside resistance:  1.0   23.95     
blaOXA-1_1_J02967       Beta-lactam resistance: 1.0     22.82   
strB_1_M96392   Aminoglycoside resistance:Alternate name; aph(6)-Id     1.0     22.32     
```

If the program is not giving any arguments or lacking some it will instead use the defaults which are:

```
-g : "Unknown3_raw_reads_1.txt.gz" 
-r : "resistance_genes.fsa" 
-k : 19
```
Example files to run the program can be found below:  
Link to "resistance_genes.fsa":  
teaching.healthtech.dtu.dk/material/36610/resistance_genes.fsa  
Link to "Unknown3_raw_reads_1.txt.gz":  
teaching.healthtech.dtu.dk/material/36610/Unknown3_raw_reads_1.txt.gz   
Running the program with these files takes about 4 minutes.  

# Program structure
The program structure, described in pseudocode can be seen below:
```
Read in the antibiotic genes file:
    for each gene, their complement and their reverse
        Save the kmers, their positions and in which gene they were found

Read in the fastaQ file:
    Go throught each kmer:
        If the kmer was present in the antibiotic genes file:
            Save its coverage
    
    Go through each gene present in the read
        if its valid based on it's coverage:
            add its coverage to a variable saving the total coverage of the gene

Go through each gene
    If its total coverage is above 95 percent and depth above 10
        print it to std.out sorted based on coverage then average depth
```
