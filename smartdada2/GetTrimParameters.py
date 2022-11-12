import warnings
import pandas as pd
from smartdada2.reader import reader


def trim_ends_less_than_threshold(avg_qs_list, threshold=30, max_trim_per=0.1):
    """Perform obvious trimming: moving from the middle index to the left
    and then to the right, the first index with a value lower than the
    threshold will be the potential trim index

    :param avg_qs_list: list of floats 
        list of average quality scores per position
    :param threshold: integer
        quality score between 0 and 42
    :param max_trim_per: float
        maximum percentage of trimming performed at either end
    :return: return tuple of two integers
        potential left and right trim values
    """
    if not isinstance(avg_qs_list, list):
        raise TypeError("average list of quality scores must be of type list")
    if not isinstance(threshold, int):
        raise TypeError("threshold must be of type integer")
    if not isinstance(max_trim_per, float):
        raise TypeError("max trim percentage must be of type float")
    if threshold < 0 or threshold > 42:
        raise ValueError("threshold must be between 0 and 42")
    for avg_score in avg_qs_list:
        if not isinstance(avg_score, float):
            raise TypeError("quality scores must be of type float")
        if avg_score < 0 or avg_score > 42:
            raise ValueError("average scores must be between 0 and 42")

    # find the middle index
    len_list = len(avg_qs_list)
    mid_index = len_list // 2

    trim_left_index = 0
    trim_right_index = 0

    # get left index before value is below threshold
    current_index = mid_index

    # get left index before value is below threshold
    while current_index >= 0:
        # if value at current index is below threshold
        if avg_qs_list[current_index] < threshold:
            # get the prior index
            trim_left_index = current_index + 1
            break
        else:
            current_index -= 1

    # get right index before value is below threshold
    current_index = mid_index
    while current_index <= len_list:

        # check for last index
        if current_index >= len_list:
            break

        # if value at current index is below threshold
        if avg_qs_list[current_index] < threshold:
            # get prior checked index
            trim_right_index = current_index
            break
        else:
            current_index += 1

    # if not obvious trimming was performed at all 
    if trim_left_index == 0 and trim_right_index == 0:
        warnings.warn("No obvious trimming performed")

    # if no obvious right trimming was performed
    if trim_right_index == 0:
        # right index is last index
        trim_right_index = len_list 

    # warn if trim values create a short read
    perc_10 = round(len_list * max_trim_per)
    if trim_left_index > perc_10:
        warnings.warn("trim left value might be too high")

    if trim_right_index < len_list - perc_10:
        warnings.warn("trim right value might be too low")

    return trim_left_index, trim_right_index


def get_trim_length_avgEE(avgEE_list, left, right):
    """This function takes a list of average expected error per position and indexing 
    sites to return a list of trim sites, read length, and average EE 

    Args:
        avgEE_list (list): average expected error per position
        left (int): left index
        right (int): right index 

    Returns:
        list: string of indexes, length of read, average EE per position
    """
    if not isinstance(avgEE_list, list): 
        raise TypeError('average EE must be of type list')
    if not isinstance(left, int):
        raise TypeError('left index must be of type int')
    if not isinstance(right, int):
        raise TypeError('right index must be of type index')
    for EE in avgEE_list:
        if not isinstance(EE, float):
            raise TypeError('all values in avgEE must be of type float')

    # to get average EE per position of read 
    current_read_len = len(avgEE_list[left:right])
    current_sumEE = sum(avgEE_list[left:right])

    # get indexes, read length, and avgEE 
    trim_indexes = str(left) + ":" + str(right)
    read_len = current_read_len
    avgEE_position = current_sumEE / current_read_len

    return [trim_indexes, read_len , avgEE_position]



def read_size_by_avg_EE(FastqEntries, left, right, max_trim_perc = 0.20): 
    """Creates a dataframe containing average expected error per position for each 
    length of read calculated within certain ranges of positions (e.g., not taking 
    out more than 20% off of either end)

    Args:
        FastqEntries (FastqEntry): reads taken from fastq files 
        left (integer): left index 
        right (integer): right index
        max_trim_perc (float, optional): max percentage of the read trimmed on 
        either end. Defaults to 0.2.

    Returns:
        pandas dataframe: dataframe containing 3 columns: trim positions, read 
        length, and average EE per position 
    """
    if not isinstance(FastqEntries, reader.FastqReader):
        raise TypeError('must be of FastqEntry type')
    if not isinstance(left, int):
        raise TypeError('left must be of type int')
    if not isinstance(right, int):
        raise TypeError('right must be of type int')
    if not isinstance(max_trim_perc, float):
        raise TypeError('max trim percentage must be of type float')

    # get average expected error per position
    avg_EE_df = FastqEntries.get_expected_error()
    avg_EE_list = list(avg_EE_df["AverageExpectedError"])

    for avgEE in avg_EE_list:
        if not isinstance(avgEE, float):
            raise TypeError('avg EEs must be of type float')

    # get length of reads
    read_len = len(avg_EE_list)

    # create lists with trim sites, read length, avg EE per position
    trim_indexes = []
    read_len_list = []
    avgEE_position = []

    # start trimming from the left
    current_left = left 
    while current_left < read_len * max_trim_perc:
        trim_readLen_avgEE = get_trim_length_avgEE(avg_EE_list, current_left, right)
        trim_indexes.append(trim_readLen_avgEE[0])
        read_len_list.append(trim_readLen_avgEE[1])
        avgEE_position.append(trim_readLen_avgEE[2])

        current_left += 1

    # start trimming from the right 
    current_right = right
    while current_right > read_len - read_len*max_trim_perc:
        trim_readLen_avgEE = get_trim_length_avgEE(avg_EE_list, left, current_right)
        trim_indexes.append(trim_readLen_avgEE[0])
        read_len_list.append(trim_readLen_avgEE[1])
        avgEE_position.append(trim_readLen_avgEE[2])

        current_right -= 1

    # start trimming on either end 
    while left < read_len*max_trim_perc and right > read_len - read_len*max_trim_perc: 
        trim_readLen_avgEE = get_trim_length_avgEE(avg_EE_list, left, right)
        trim_indexes.append(trim_readLen_avgEE[0])
        read_len_list.append(trim_readLen_avgEE[1])
        avgEE_position.append(trim_readLen_avgEE[2])

        left += 1
        right -= 1

    # create three column pandas data frame
    df = pd.DataFrame(list(zip(trim_indexes, read_len_list, avgEE_position)),
                    columns=['Indexes', 'ReadLength', 'AvgEEPerPosition'])
    return df 


    # TODO 
    # CREATE FINAL FUNCTION FOR CALCULATION OF TRIM PARAMS 
