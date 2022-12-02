import argparse as arg
import warnings
from smartdada2.reader import reader
import GetTrimParameters as GTP
import GetMaxEE as GME


def main():

    parser = arg.ArgumentParser(
        description="This script will output a dataframe containing"
        + "trimming information read length, average expected error by"
        + "position and the number of reads below and over a threshold"
    )

    parser.add_argument(
        "--fq",
        "-fastq",
        type=str,
        required=True,
        help="file path to fastq files",
    )
    parser.add_argument(
        "--of",
        "-output_file",
        type=str,
        required=False,
        help="name of the output file",
        default="parameter_info.tsv",
    )
    parser.add_argument(
        "--th",
        "-threshold",
        type=int,
        required=False,
        help="Phred score threshold for obvious trimming",
        default=30,
    )
    parser.add_argument(
        "--o_mtp",
        "-obvs_max_trim_percentage",
        type=float,
        required=False,
        help="max percent of the length of the original trim to be"
        + "trimmed on either end of the read in obvious trimming",
        default=0.1,
    )
    parser.add_argument(
        "--a_mtp",
        "-additional_max_trim_percentage",
        type=float,
        required=False,
        help="when performing further trimming the max percentage of"
        + "the read to be trimmed at either end",
        default=0.2,
    )
    parser.add_argument(
        "--mEE",
        "-maxEE",
        type=float,
        required=False,
        help="max expected error default is set to same value"
        + "as default DADA2 parameters",
        default=2.0,
    )

    args = parser.parse_args()

    # ensure o_mtp is less than a_mtp
    if args.o_mtp > args.a_mtp:
        warnings.warn(
            "max percentage of trimming done in obvious trimming is larger"
            + "than the max in following steps which is not recommended"
        )

    # read in data
    fqe = reader.FastqReader(args.fq, scale=239)

    # get average scores per position
    avg_scores_df = fqe.get_average_score()
    avg_scores_list = avg_scores_df["AverageQualityScore"].values.tolist()

    # perform obvious trimming
    left, right = GTP.trim_ends_less_than_threshold(
        avg_scores_list, args.th, args.o_mtp
    )

    # get average EE by position for different read sizes
    EE_by_size_df = GTP.read_size_by_avg_EE(fqe, left, right, args.a_mtp)

    # add number of reads over and below the threshold to the dataframe
    reads_over_under = GME.read_size_by_maxEE(fqe, EE_by_size_df, args.mEE)

    # convert final dataframe to TSV
    reads_over_under.to_csv(args.of, sep="\t", index=False)

    exit()


if __name__ == "__main__":
    main()
