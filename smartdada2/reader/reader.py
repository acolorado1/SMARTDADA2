import warnings
import itertools
from pathlib import Path
from dataclasses import dataclass
from typing import Iterable
from typing import Optional
from typing import Union
from typing import List

import pandas as pd
import numpy as np

from smartdada2.seq_utils.search import binary_search
from smartdada2.common.errors import FastqFormatError


# constants
DNA = list("ATCGU")
AMB_DNA = list("NRYKMSWBDHV")
ALL_DNA = DNA + AMB_DNA
ASCII_SCORES = [chr(i) for i in range(33, 70 + 1)]
NUMERICAL_SCORES = [i for i in range(33, 70 + 1)]


@dataclass(slots=True)
class FastqEntry:
    """Contains the contents of a single read as a FastqEntry"""

    header: str
    seq: str
    scores: str
    length: int

    rseq: Union[None, bool] = None

    def __post_init__(self):

        # checking header (1 == forward, 2 == rev)
        self.check_seq_dir()

    def check_seq_dir(self):
        """Checks if the entry is a forward or reverse sequence"""

        # parse the header
        header_id = self.header.split()[0]
        if not header_id.startswith("@"):
            raise FastqFormatError("Unable to find header id")

        # find the direction of the sequence
        direction = header_id.split(".")[-1]

        # convert string into integer type
        if not direction.isdigit():
            raise FastqFormatError("Unable to find the sequence direction")
        else:
            direction = int(direction)

        # setting up sequence direction into FastqEntry
        if direction == 1:
            self.rseq = False
        elif direction == 2:
            self.rseq = True
        else:
            raise FastqFormatError("Unable to find the sequence direction")

        # making sequences to be capital letters
        self.seq = self.seq.upper()


class FastqReader:
    """
    FastqReader is an efficient Fastq file reader that optimizes for memory
    usage and speed.

    Every entry within the Fastq file is stored within a FastqEntry object,
    which contains header, sequence, sequence length, score, and direction of
    each read.

    FastqReader contains methods that are meant for gathering information
    within fastq files.
    """

    def __init__(
        self,
        fpath: str,
        reverse_seq: Optional[bool] = False,
        technology: Optional[str] = "illumina",
    ):

        # FastqReader accessible parameters
        self.fpath: Path = Path(fpath)
        if not self.fpath.is_file():
            raise FileNotFoundError(f"{self.fpath} does not exist")

        self.reversed: bool = reverse_seq
        self.technology: str = technology

        # states
        self.__n_entries: Optional[int] = 0
        self.__counted: bool = False

    def get_quality_scores(self) -> pd.DataFrame:
        """Returns scores in a per sequence bases

        Returns
        -------
        pd.DataFrame
            DataFrame containing quality scores per sequence

        """

        # converting to np.array
        all_scores = []
        for entry in self.iter_reads():
            phred_scores = list(entry.scores)
            scores = np.array(phred_scores).view(np.int32) - 33
            all_scores.append(scores)

        all_scores = pd.DataFrame(data=np.array(all_scores))
        return all_scores

    def get_average_score(self) -> pd.Series:
        """Returns average score of all sequences. Returns a a pd.Series object

        Returns
        -------
        pd.Series
            average score of all sequences

        """

        # get all scores df and take the average per column basis
        scores_df = self.get_quality_scores()
        average_score = scores_df.mean().reset_index()
        average_score.columns = ["Position", "AverageQualityScore"]
        return average_score

    def get_expected_error(self) -> pd.DataFrame:
        """Calculates average expected error per nucleotide

        Returns
        -------
        pd.DataFrame
            position and expected error in a dataframe
        """
        # loading average quality scores
        scores_df = self.get_average_score()

        # convert quality score values into expected errors values
        scores_df["AverageExpectedError"] = scores_df[
            "AverageQualityScore"
        ].apply(lambda score: 10 ** (-score / 10))
        expected_error_df = scores_df.drop(columns="AverageQualityScore")

        return expected_error_df

    # TODO: finish documentation
    def get_max_ee(self) -> pd.Series:
        """Creates a pandas dataframe object that contains the sum of expected
        error of all sequences. The index will represent the sequence. The two
        columns "direction" and "total_error"

        Returns
        -------
        pd.DataFrame
            A two column dataframe that contains "direction" and "total_error"
        """

        # create empty df
        expected_error_df = pd.DataFrame()

        # quality scores
        qual_scores_df = self.get_quality_scores()
        raw_scores = qual_scores_df.apply(
            lambda row: 10 ** (-row / 10), axis=1
        )

        # extract read direction
        _dir = []
        for read in self.iter_reads():
            if read.rseq is True:
                _dir.append("forward")
            else:
                _dir.append("reverse")

        expected_error_df["max_ee"] = raw_scores.apply(
            lambda row: np.sum(row), axis=1
        )

        expected_error_df["length"] = [
            entry.length for entry in self.iter_reads()
        ]
        expected_error_df["direction"] = _dir

        # rearranging columns
        expected_error_df = expected_error_df[
            ["length", "direction", "max_ee"]
        ]

        return expected_error_df

    # TODO: DRY practice found here, need to create a function that collects scores
    def get_avg_ee(self) -> pd.DataFrame:
        """Get average expected error per sequence. Generates a pandas
        dataframe that contains sequence length, direction, and average error
        score

        Returns
        -------
        pd.DataFrame
            average expected error per sequence dataframe
        """
        # create empty df
        avg_expected_error_df = pd.DataFrame()

        # quality scores
        qual_scores_df = self.get_quality_scores()
        raw_scores = qual_scores_df.apply(
            lambda row: 10 ** (-row / 10), axis=1
        )

        # extract read direction
        _dir = []
        for read in self.iter_reads():
            if read.rseq is True:
                _dir.append("forward")
            else:
                _dir.append("reverse")

        avg_expected_error_df["avg_ee"] = raw_scores.apply(
            lambda row: np.mean(row), axis=1
        )

        avg_expected_error_df["length"] = [
            entry.length for entry in self.iter_reads()
        ]
        avg_expected_error_df["direction"] = _dir

        # rearranging columns
        avg_expected_error_df = avg_expected_error_df[
            ["length", "direction", "avg_ee"]
        ]

        return avg_expected_error_df

    def sequence_df(self) -> pd.DataFrame:
        """
        Returns a Dataframe structure of sequence reads. Row represents a
        sequence and the columns represents the individual nucleotides
        """
        seq_df = pd.DataFrame(
            data=(list(entry.seq) for entry in self.iter_reads())
        )
        return seq_df

    def ambiguous_nucleotide_counts(self) -> pd.DataFrame:
        """Counts all ambiguous nucleotides in all reads.

        Returns
        -------
        pd.DataFrame
            DataFrame that contains the nucleotide position and the number
            of ambiguous nucleotides
        """
        seq_df = self.sequence_df()
        ambi_nuc_counts = (
            seq_df.apply(search_ambiguous_nucleotide).to_frame().reset_index()
        )
        ambi_nuc_counts.columns = ["Position", "AmbiguousCounts"]
        return ambi_nuc_counts

    def total_reads(self) -> int:
        """Returns total number of reads. Changes the value of self.n_entries
        to prevent any recalculation

        Returns
        -------
        int
            total number of reads
        """
        # checks if the n_entries have been counted before
        # -- if not, count and return value
        # -- if yes, return stored count
        if self.__n_entries == 0 and self.__counted is False:
            warnings.warn("No reads have been counted. Counting reads now...")
            self.__n_entries = sum(1 for _ in self.iter_reads())
            self.__counted = True
            return self.__n_entries
        else:
            return self.__n_entries

    def get_forward_reads(self) -> Iterable[FastqEntry]:
        """Obtains all forward reads within fastq files

        Parameters
        ----------
        tolist : bool, optional
            Saves the forward reads in memory (list), by default False

        Returns
        -------
        Iterable[FastqEntry]
            returns a generator if tolist is False. if tolist is True, then the
            reads will be saved into memory and will return a python list
            object.
        """

        # checking if the user wants to load reads in memory
        for read in self.__loader():
            if read.rseq is False:
                yield read

    def get_reverse_reads(self, tolist=False) -> Iterable[FastqEntry]:
        """Obtains all reverse reads within fastq files

        Parameters
        ----------
        tolist : bool, optional
            Saves the reverse reads in memory (list), by default False

        Returns
        -------
        Iterable[FastqEntry]

        """
        for read in self.__loader():
            if read.rseq is True:
                yield read

    def sample(
        self,
        frac: Optional[float] = 0.3,
        seed: Optional[int] = None,
    ) -> Iterable[FastqEntry]:
        """Sub samples sequences randomly selected. If a seed value is provided,
        the randomness is controlled and always produces the same output order
        associated with a specific seed value.

        Parameters
        ----------
        frac : Optional[float], optional
            subsample size, by default 0.3
        seed : Optional[int], optional
            Value that retains the state randomness, by default None

        Yields
        ------
        Iterator[Iterable[FastqEntry]]
            Generator object that contains FastqEntry

        Raises
        ------
        ValueError
            Raised when either frac or prob is larger than 1.0
        """
        # type checks
        if not isinstance(frac, float):
            frac = float(frac)

        # selecting sampling type

        # checking if fraction subsample is larger than the whole dataset
        if frac > 1.0:
            raise ValueError("frac cannot be larger than 1.0")

        # checking if reads have been counted
        # this will update the `self.__n_entries` variable
        if self.__counted is False:
            self.total_reads()

        # initializing random seed if Seed is not None
        # -- setting the seed will have a controlled randomness
        # -- good for reproducibility (share seed value)
        if seed is not None:
            np.random.seed(seed)

        # calculate subsample_size (how many entries we want)
        subsample_size = int(np.round(self.__n_entries * frac))

        # selecting random indices without replacement
        # range selection is (0, total_n_reads)
        # -- done without replacement (no repeated indices selected)
        random_indices = np.sort(
            np.random.choice(
                np.arange(0, self.__n_entries),
                size=subsample_size,
                replace=False,
            )
        )

        # to stop further computation once all indices are found
        stop_idx = np.max(random_indices)

        # loading the selected reads
        for idx, read_entry in enumerate(self.iter_reads()):

            # stop if next read idx is larger than stop idx
            if idx > stop_idx:
                break

            # if index exists within randomly selected indices, yield read
            try:
                binary_search(idx, random_indices, sorted=True)
                yield read_entry

            # if it does not exists, go to the next
            except ValueError:
                continue

    def reservoir_sampling(
        self, n_samples: int, seed: Optional[bool] = None
    ) -> list[FastqEntry]:
        """Generates a uniform random subset sample of FastqEntry.

        Parameters
        ----------
        n_samples : int
            size of sub sample
        seed : Optional[int], optional
            Value that retains the state randomness, by default None

        Returns
        -------
        list[FastqEntry]
            subset of FastqEntry objects
        """

        # type checking
        if seed is not None and not isinstance(seed, int):
            raise ValueError(
                f"seed value must be an integer type, you provided {type(seed)}"
            )

        # set a results array and have elements inside as place holders
        # -- seen is equals to n_samples because those elements make up the results array
        seen = n_samples
        subset_reads = []

        # populating results array with FastqEntry objects as place holders

        # -- checking if number of entries have counted before
        # -- -- if yes, set slice of FastqEntries as placeholders
        if self.__counted is True:
            if self.__n_entries < n_samples:
                raise ValueError(
                    "requested sample size is larger than total number of entries"
                )
            else:
                subset_reads = list(
                    itertools.islice(self.iter_reads(), n_samples)
                )

        # -- if size is unknown, iterate through loader and populate results
        # with FastqEntry placeholders
        else:
            try:
                entries = self.iter_reads()
                for _ in range(n_samples):
                    subset_reads.append(next(entries))

            # raise an exception if the subsample size is larger than the loader population
            except StopIteration:
                raise StopIteration(
                    "number of requested samples exceeded number of entries"
                )

        # next is to update the results array
        # -- checking if seed has been added
        if seed is not None:
            np.random.seed(seed)

        # -- now replacing placeholders by randomly selecting values as ind
        for entry in self.iter_reads():

            # select random int
            seen += 1
            rand_int = np.random.randint(0, seen - 1)

            # replace placeholders via index reassignment
            if rand_int < n_samples:
                subset_reads[rand_int] = entry

        return subset_reads

    @staticmethod
    def unpack_reads(iterator: Iterable[FastqEntry]) -> List[FastqEntry]:
        """
        Takes the iterator object produced from FastqReader and saves them into
        memory as a python list.
        """
        return [read for read in iterator]

    def iter_reads(self) -> Iterable[FastqEntry]:
        """Returns a python generator containing FastqEntries

        Returns
        -------
        Sequence[FastqEntry]
            Generator object containing FastqEntries
        """
        return self.__loader()

    def slice_reads(
        self, range_idx: Union[tuple[int, int], list[int, int]], to_list=False
    ) -> Union[list[FastqEntry], Iterable[FastqEntry]]:
        """Allows index slicing for reads. If to_list is True, the method will
        return a list else it will return a generator

        Parameters
        ----------
        range_idx : tuple
            start and ending position (exclusive), to which reads to select.
        to_list : bool, optional
            convert to list object, If False, returns a generator.
            by default False

        Returns
        -------
        Union[list[FastqEntry], Iterable[FastqEntry]]
        return a list or iterator of FastqEntry objects

        Yields
        ------
        Iterator[Union[list[FastqEntry], Iterable[FastqEntry]]]
            if to_list is false, returns a generator of reads, else it will
            return a list of reads

        Raises
        ------
        TypeError
            Raised when a tuple or list is not provided
        ValueError
            Raised if range_idx does not contains 2 values
        TypeError
            Raised if values in range_idx are not integers
        ValueError
            raised if starting position value is larger than the ending
            position
        """
        if not isinstance(range_idx, tuple) and not isinstance(
            range_idx, list
        ):
            raise TypeError(
                "Please provide a tuple or lists with starting and ending idx"
            )
        elif len(range_idx) != 2:
            raise ValueError("'range_idx' only takes two value (start, end)")
        elif not all(isinstance(value, int) for value in range_idx):
            raise TypeError(
                "Values must be integers. Not floats, strings or booleans"
            )
        elif range_idx[0] > range_idx[1]:
            raise ValueError(
                "starting position cannot be larger than ending position"
            )

        if to_list == True:
            return [entry for entry in self.__slice(range_idx)]
        else:
            return self.__slice(range_idx)

    # ----------------------------------------
    # private functions: users do not interact with this
    # ----------------------------------------
    def __loader(self) -> Iterable[FastqEntry]:
        """Creates a generator object that streams sequence reads as
        FastqEntry objects.

        While iterating, a count is also being conducted.

        Returns
        -------
        Sequence[FastqEntry]
            Generator object with FastqEntries
        """

        if self.fpath.suffix.lower() != ".fastq":
            raise ValueError(
                "FastqReader only takes files '.fastq' or '.FASTQ' files"
            )
        elif self.fpath.stat().st_size == 0:
            raise FastqFormatError("Fastq file contains no contents")

        # iterate all row contents in fastq file and collect entries
        entry_count = 0
        with open(self.fpath, "r") as fastq_file:

            contents_chunk = []
            for row_entry in fastq_file:

                # cleaning entries
                if row_entry == "":
                    continue

                content = row_entry.rstrip("\n")
                contents_chunk.append(content)

                # checking if there are 4 elements in the list
                # -- 4 lines = 1 entry
                if len(contents_chunk) == 4:

                    # check for valid sequences
                    seq_check = set(contents_chunk[1]) - set(ALL_DNA)
                    if len(seq_check) > 0:
                        raise FastqFormatError(
                            "File contains invalid sequence characters"
                        )

                    # check for valid ascii scores
                    score_check = set(contents_chunk[3]) - set(ASCII_SCORES)
                    if len(score_check) > 0:
                        raise FastqFormatError(
                            "File contains invalid score characters"
                        )

                    # convert into FastqEntry
                    entry_count += 1
                    fastq_entry = FastqEntry(
                        header=contents_chunk[0],
                        seq=contents_chunk[1],
                        scores=contents_chunk[3],
                        length=len(contents_chunk[1]),
                    )

                    # clear list
                    contents_chunk = []

                    # yield entry
                    yield fastq_entry

        self.__n_entries = entry_count
        self.__counted = True

    def __slice(self, range_idx: tuple[int, int]) -> Iterable[FastqEntry]:
        """ "Creates a generator of py FastqEntry object by a given range

        Parameters
        ----------
        idx_range : tuple
            start, end range of selecting FastqEntries

        Returns
        -------
        Iterable[FastqEntry]
            Generator with FastqEntry objects

        Yields
        ------
        Iterator[Iterable[FastqEntry]]
            FastqEntry objects
        """

        # unpacking range
        start, end = range_idx

        # iterating reads
        for idx, read in enumerate(self.iter_reads()):

            # break the iteration if the current idx == ending idx
            if idx == end:
                break

            # only selecting reads between start and end positions
            elif idx >= start and idx < end:
                yield read


# NOTE: Should this be in the FastqReader class or a separate function?
def search_ambiguous_nucleotide(nucleotides: pd.Series) -> int:
    """
    Parameters
    ----------
    nucleotides : pd.Series
        Series of nucleotides characters.

    Returns
    -------
    int
        total number of all ambiguous nucleotides found.
    """

    # counting all nucleotides
    ambiguous_nucleotide = list("NRYKMSWBDHV")
    count_series = nucleotides.value_counts()

    # searching for  nucleotides
    found_ambiguous_nucleotide = []
    for ambi_nuc in ambiguous_nucleotide:
        if ambi_nuc in count_series.index.tolist():
            found_ambiguous_nucleotide.append(ambi_nuc)

    total_count = count_series[found_ambiguous_nucleotide].sum()

    return total_count
