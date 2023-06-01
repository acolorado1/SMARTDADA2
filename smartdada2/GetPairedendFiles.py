import argparse as arg

from smartdada2.reader import reader


def getPairedFiles(main_fastq: str, output_f: str, output_r: str):
    """If paired end sequencing was performed get
    two files containing reverse and forward sequences.

    Args:
        main_fastq (str): one FASTq formatted file containing all
        reads
        output_f (str): output file path for forward sequences
        output_r (str): output file path for reverse sequences

    Returns:
        files: two fastq files will be written containing
        forward and reverse reads separately.
    """
    # make sure input is of right type
    if not isinstance(main_fastq, str):
        raise TypeError("must be a file path of type string")

    # read in data
    fqe = reader.FastqReader(main_fastq)

    # make sure data is correct type
    if not isinstance(fqe, reader.FastqReader):
        raise TypeError("must be of FastqEntry type")

    # get forward reads in list format
    forward = []
    for reads in fqe.get_forward_reads():
        forward.extend((reads.header, reads.seq, "+", reads.scores))

    # write forward file
    with open(output_f, "+w") as f:
        for line in forward:
            line = str(line) + "\n"
            f.write(line)

    # get reverse reads in list format
    reverse = []
    for reads in fqe.get_reverse_reads():
        reverse.extend((reads.header, reads.seq, "+", reads.scores))

    # write reverse file
    with open(output_r, "+w") as f:
        for line in reverse:
            line = str(line) + "\n"
            f.write(line)

    return None


# created parser to run through terminal
def main():
    parser = arg.ArgumentParser(
        description="This script will separate forward and reverse "
        + "reads for either paired or single end sequencing"
    )

    parser.add_argument(
        "--fq",
        "-fastq",
        type=str,
        required=True,
        help="file path to fastq files",
    )

    parser.add_argument(
        "--of_f",
        "-output_file_forward",
        type=str,
        required=True,
        help="output file path for forward sequences",
    )

    parser.add_argument(
        "--of_r",
        "-output_file_reverse",
        type=str,
        required=True,
        help="output file path for reverse sequences",
    )

    args = parser.parse_args()

    getPairedFiles(args.fq, args.of_f, args.of_r)

    exit()


if __name__ == "__main__":
    main()
