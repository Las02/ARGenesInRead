

import sys

def argument_parser(argv,kmer_length, gene_file, read_file):

    last_arg = ""

    # Go throug argv. If the last argument was a "-" argument, save the argument in correct variable
    for arg in argv:

        # Exit if missing an argument
        if last_arg.startswith("-") and arg.startswith("-"):
            sys.exit(f"Missing argument in {last_arg}")

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




(k,g,r) = argument_parser(sys.argv,"1a9","file1","file2")

print(k,g,r)

