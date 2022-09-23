from pathlib import Path
from dataclasses import dataclass
from typing import Sequence
from typing import Optional
from typing import Union
from typing import Optional

import pandas as pd
import numpy as np


@dataclass(slots=True)
class FastqEntry:
    """Contains the contents of a single read as a FastqEntry"""

    header: str
    seq: str
    scores: str
    length: int

    # set to none unless user declares there is reversed sequences
    rseq: Union[None, int] = None


class FastqReader:
    """Memory efficient FastqReader"""

    def __init__(
        self,
        fpath: str,
        reverse_seq: Optional[bool] = False,
        technology: Optional[str] = "illumina",
    ):

        # FastqReader accessible parameters
        self.fpath: Path = Path(fpath)
        self.reversed: bool = reverse_seq
        self.technology: str = technology
        self.n_entries: Optional[int] = 0

        # states
        self.__counted: bool = False

    def get_quality_scores(self) -> pd.DataFrame:
        """Returns scores in a per sequence bases

        Returns
        -------
        pd.DataFrame
            Dataframe containing quality scores per sequence

        """

        # converting to np.array
        all_scores = []
        for entry in self.__loader():
            phred_scores = list(entry.scores)
            scores = [ord(phred_score) - 33 for phred_score in phred_scores]
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
        average_score = scores_df.mean()
        return average_score

    def iter_reads(self) -> Sequence[FastqEntry]:
        """Returns a python generator containing FastqEntries

        Returns
        -------
        Sequence[FastqEntry]
            Generator object containing FastqEntries
        """
        return self.__loader()

    def to_list(self) -> list[FastqEntry]:
        """Saves all lists into memory. Warning, large file sizes will use more
        memory but increase iteration performance.

        Returns
        -------
        list[FastqEntry]
            List of FastqEntries
        """
        return [read for read in self.__loader()]

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
        if self.n_entries == 0 and self.__counted is False:
            self.n_entries = sum(1 for _ in self.__loader())
            self.__counted = True
            return self.n_entries
        else:
            return self.n_entries

    # ----------------------------------------
    # private functions: users do not interact with this
    # ----------------------------------------
    def __loader(self) -> Sequence[FastqEntry]:
        """Creates a generator object that contains sequence read data as
        FastqEntry object.

        While iterating, a count is also being conducted

        Returns
        -------
        Sequence[FastqEntry]
            Generator object with FastqEntries
        """

        # iterate all row contents in fastq file
        entry_count = 1
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

                    # convert into FastqEntry
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


