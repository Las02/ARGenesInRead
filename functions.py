import sys
import gzip

def argument_parser(argv, kmer_length, gene_file, read_file):
    """
    Sets up a argument parser from the argv vector. 
    In addition takes in default values for variables, if not given 
    arguments in argv
    """
    last_arg = ""
    # Go throug argv. If the last argument was a "-" argument (eg -k), save the argument in correct variable
    for arg in argv:
        # Exit if missing an argument
        if last_arg.startswith("-") and arg.startswith("-"):
            sys.exit(f"Missing argument in {last_arg}")
        elif arg.startswith("-") and argv[-1] == arg:
            sys.exit(f"Missing argument in {arg}")

        elif last_arg == "-k":
            try:
                kmer_length = int(arg)
            except ValueError as e:
                sys.exit(f"The kmer needs to be an integer, error:\n{e}")

        elif last_arg == "-g":
            gene_file = arg

        elif last_arg == "-r":
            read_file = arg

        last_arg = arg

    return (kmer_length, gene_file, read_file)

def read_fasta(filename):
    '''Reading in several fasta files
       yield each pair of header and dna seperate
    '''
    # Try to open the file, error if file does not exsist
    try: file = open(filename, "r")
    except FileNotFoundError as errormessage:
        sys.exit(f"The file '{filename}' could not be found, error: {errormessage}")

    # Extract the dna, and the headers.
    oldheader = "FIRST HEADER"
    for line in file:
        line = line.strip()
        #if line is header yield dna and header except FIRST HEADER
        if line.startswith(">"):    
            newheader = line
            if oldheader != "FIRST HEADER":
                yield dna, oldheader
            dna = ""
            oldheader = newheader
        else:
            dna += line
    # Yield the last header and dna
    yield dna, oldheader
    file.close()
    
def read_fastq(filename):
    '''Extract dna from a gzipped fastQ file '''
    # Try to open the file, error if file does not exsist
    try: sample_file = gzip.open(filename, "r")
    except FileNotFoundError as errormessage:
        sys.exit(f"The file '{filename}' could not be found, error: {errormessage}")

    last_line = ""
    try:
        for line in sample_file:
            line = line.decode("utf-8")

            # Extracting the dna
            if last_line.startswith("@"):
                dna = line.strip()
                yield dna

            last_line = line
    
    # Error message if the gzipped file does not work
    except gzip.BadGzipFile as errormessage:
        sys.exit(f"The file '{filename}' could not be gzipped, error: {errormessage}")
    sample_file.close()

def get_kmers_and_pos(dna, kmer_len):
    '''Find Kmers from a dna string and return list with them
       in additon to list with their positions in the dna string
       '''
    # Return empty lists if there is no possible kmer's
    if kmer_len > len(dna):
        return [],[]

    #  Find all possible kmers
    kmer_list = []
    range_list = []
    to_range = len(dna) - kmer_len + 1   #+1 to also add last kmer
    for i in range(0, to_range):
        kmer = dna[i: i + kmer_len]
        kmer_list.append(kmer)
        range_list.append(i)

    return kmer_list, range_list

def read_is_valid(depht_count, dna_read):
    '''
    return True if the read is covering enough of the gene,
    else returns False
    The input is a list of ints representing the depth of each nucleotide. 
    '''
    # Set parameters
    MAX_SPACE = 1   # The maximum space allowed in the read covering the gene
    SIDE_BONUS = int(0.70 * len(dna_read))   # How much of the read can be outside the gene and still be considered valid
    THRESHOLD_SCORE = int(0.95 *len(dna_read))   # The percent of the read which needs to cover the gene

    # Find the reads largest continuous coverage of the gene.
    # a coverage is considered continuous if the spacing does not exceed the MAX_SPACE 
    END_OF_GENE = len(depht_count) -1
    (maxcount, score, exitcount) = (0,0,0)
    in_coverage = False
    for pos,value in enumerate(depht_count):
        # Start counting the coverage when a depht != 0
        if value == 1: in_coverage = True
        # Start counting the coverage
        if in_coverage:
            # If a position is not covered, start counting to exit
            if value == 0: exitcount += 1
            # Else add to the coverage score for each position
            else: score += 1
            # If at either side position of the gene, add a bonus to the score
            # The positions are one step from either ends to still get SIDE_BONUS when pos 0 is not covered but pos 1 is.
            if pos in [0+1, END_OF_GENE -1]:  
                score += SIDE_BONUS
            # update the maxscore for the read
            if score > maxcount: maxcount = score
            # Stop counting the coverage if spacing exceeds MAX_SPACE
            if exitcount >= MAX_SPACE:
                score = 0
                exitcount = 0
                in_coverage = False

    return maxcount >= THRESHOLD_SCORE

def coverage_stats(depht_count):
    '''from depht array return coverage, avg depht and min_depht
    '''

    # Find total_depht, coverage and avg_depht
    count = 0
    total_depht = 0
    min_depth = 99999
    for base in depht_count:
        total_depht += base
        if base != 0:
            count += 1
        if base < min_depth:
            min_depth = base
    coverage = count / len(depht_count)
    avg_depht = total_depht / len(depht_count)
    return coverage, avg_depht, min_depth
