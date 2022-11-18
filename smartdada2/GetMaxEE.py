import pandas as pd
import numpy as np
import smartdada2.GetTrimParameters as GTP
from smartdada2.reader import reader


def read_size_by_maxEE(FastqEntries, TrimInfo: pd.DataFrame, maxEE=2.0):
    """Take FastqEntries and table output from the read_size_by_avg_EE
    and count number of reads that are above and below the maxEE threshold.

    Args:
        FastqEntries (FastEntry): reads taken from Fastq files
        TrimInfo (pd.DataFrame): dataframe containing 3 columns
        maxEE (float, optional): Max sum of the expected errors a read
                                contains. Defaults to 2.0 because
                                DADA2 uses this as a default.

    Returns:
        pd.DataFrame: dataframe containing 5 columns
    """
    if not isinstance(maxEE, float):
        raise TypeError("maxEE must be of type float")

    # ensure input TSV is of the correct format 
    expected_cols = ["Indexes", "ReadLength", "AvgEEPerPosition"]
    if (TrimInfo.columns != expected_cols).all():
        raise ValueError("dataframe does not contain the right column names")

    # get list of expected errors per read
    raw_EE = []

    # for each entry
    for entry in FastqEntries.unpack_entries():
        phred_scores = list(entry.scores)

        # convert scores to numeric format to integer 0-42
        scores = np.array(phred_scores).view(np.int32) - 33

        # converts phred scores to expected error 
        EE_entry = [10 ** (-score / 10) for score in scores]

        # create list of lists of expected error
        raw_EE.append(EE_entry)

    # get list of posible indexes
    indexes = TrimInfo["Indexes"].values.tolist()

    # for trim values get sum of EE for read count if under or over threshold
    under_maxEE = []
    over_maxEE = []

    # for each trim values
    for index in indexes:
        index = index.split(":")
        under_count = 0
        over_count = 0

        # get the sum of EE
        for EE_list in raw_EE:
            indexed_EE = EE_list[int(index[0]):int(index[1])]
            sum_EE = sum(indexed_EE)

            # if under threshold
            if sum_EE < maxEE:
                under_count += 1
            # if over threshold
            else:
                over_count += 1

        under_maxEE.append(under_count)
        over_maxEE.append(over_count)

    # add lists as columns to dataframe
    TrimInfo["ReadsUnderMaxEE"] = under_maxEE
    TrimInfo["ReadsOverMaxEE"] = over_maxEE

    expected_colnames = [
        "Indexes",
        "ReadLength",
        "AvgEEPerPosition",
        "ReadsUnderMaxEE",
        "ReadsOverMaxEE",
    ]
    if (TrimInfo.columns != expected_colnames).all():
        raise ValueError("output dataframe not as expected")

    return TrimInfo
