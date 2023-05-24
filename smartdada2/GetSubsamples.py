import argparse as arg
from smartdada2.reader import reader

def subsample(input_fp:str, subsample:int, output_fp:str):
    """Takes one FASTQ file and outputs subsets of the reads. 

    Args:
        input_fp (str): FASTQ file path 
        subsample (int): number of reads in subset 
        output_fp (str): FASTQ file subset output path 
    
    Returns: 
        A fastq formatted file containing a subset of the original
        reads.
    """
     # make sure input is of right type
    if not isinstance(input_fp, str):
        raise TypeError("must be a file path of type string")
    if not isinstance(subsample, int):
        raise TypeError("subsampling must be of integer type")
    if not isinstance(output_fp, str):
        raise TypeError("must be a file path of type string")

    # read in data
    fqe = reader.FastqReader(input_fp)

    # make sure data is correct type
    if not isinstance(fqe, reader.FastqReader):
        raise TypeError("must be of FastqEntry type")

    # perform reservoir sampling 
    samples = fqe.reservoir_sampling(n_samples=subsample,seed=5)
    
    # create list of samples 
    sample_list = []
    for sample in samples: 
        sample_list.extend((sample.header, sample.seq, "+", sample.scores))

    # write into a fastq formatted file 
    with open(output_fp, "+w") as f:
        for line in sample_list:
            line = str(line) + "\n"
            f.write(line)
    
    return None

# created parser to run through terminal
def main():
    parser = arg.ArgumentParser(
        description="This script will take a fastq formatted file"
        + "and output a subset of it"
    )

    parser.add_argument(
        "--fq",
        "-fastq",
        type=str,
        required=True,
        help="file path to fastq files",
    )

    parser.add_argument(
        "--n", 
        "-n_sample", 
        type=int, 
        required=True, 
        help="Number of samples in subset"
    )

    parser.add_argument(
        "--of", 
        "-output_file",
        type=str, 
        required=True, 
        help="output file path"
    )



    args = parser.parse_args()

    subsample(args.fq, args.n, args.of)

    exit()

if __name__ == "__main__":
    main()