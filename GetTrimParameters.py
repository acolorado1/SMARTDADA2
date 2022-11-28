import warnings


def trim_ends_less_than_threshold(avg_qs_list, threshold=20, max_trim_per=0.1):
    """Perform obvious trimming: moving from the middle index to the left
    and then to the right, the first index with a value lower than the
    threshold will be the potential trim index

    :param avg_qs_list: list of integers
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
        raise ValueError("threshold must be between 0 and 40")
    for avg_score in avg_qs_list:
        if not isinstance(avg_score, int):
            raise TypeError("quality scores must be of type int")
        if avg_score < 0 or avg_score > 42:
            raise ValueError("average scores must be between 0 and 40")

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

        # if value at current index is below threshold
        if avg_qs_list[current_index] < threshold:
            # get prior checked index
            trim_right_index = current_index
            break
        else:
            current_index += 1

    # if no obvious right trimming was performed
    if trim_right_index == 0:
        trim_right_index = len_list

    # exit if trim values create a short amplicon
    perc_10 = round(len_list * max_trim_per)
    if trim_left_index > perc_10:
        warnings.warn("trim left value might be too high")

    if trim_right_index < len_list - perc_10:
        warnings.warn("trim right value might be too low")

    return trim_left_index, trim_right_index
