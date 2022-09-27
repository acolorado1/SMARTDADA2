"""
search.py are modules  

* binary_search() - searches if a specific integer value exists in an array
* key_binary_search() - searches if a specific key (string value) exists 
"""
import numpy as np
from typing import Union
from typing import Optional
from typing import Sequence

from smartdada2.common.errors import SearchExceededError


def binary_search(
    target: int,
    arr: Union[list, np.ndarray[int]],
    max_iter: Optional[int] = 10000,
    sorted: Optional[bool] = False,
):
    """Applies binary search to find specific value in an array.

    Parameters
    ----------
    target : int
        target value to find in array.
    arr : Union[int, np.ndarray]
        array like object to search in.
    max_iter : Optional[int], optional
        Max amount of iteration when searching. Raises StopIteration error if
        exceeds number of searches , by default 10000
    sorted : Optional[bool], optional
        indicates whether array is sorted. False will automatically sort the list
        . If set to True, list will not be sorted, default is False

    Returns
    -------
    bool
        Returns True if the key is found. raises Value error if key is not
        found.

    Raises
    ------
    TypeError
        Raised if a non array object is passed
    StopIteration
        Raised if search exceed number of iterations
    ValueError
        Raised if the target value cannot be found in given array.
    """

    #  type checking
    if isinstance(arr, Sequence):
        arr = np.array(arr)
    elif not isinstance(arr, np.ndarray):
        raise TypeError("Requires array to list type object")

    # sort array
    if sorted is False:
        arr = np.sort(arr)

    # set index positions
    low_idx = -1
    high_idx = len(arr)

    # start search
    search_count = 0
    while high_idx - low_idx > 1:

        # count search
        search_count += 1

        # if iteration is maxed, raise value error
        if search_count == max_iter:
            raise SearchExceededError("Search iteration exceeded")

        # set mid index position
        mid_idx = (high_idx + low_idx) // 2

        # if mid idx points to target value, return True
        if target == arr[mid_idx]:
            return True

        # update mid index position
        if target < arr[mid_idx]:
            high_idx = mid_idx
        else:
            low_idx = mid_idx

    raise ValueError("Unable to find target value in provided array")
