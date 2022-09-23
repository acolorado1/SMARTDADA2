import pandas as pd
from pathlib import Path
from typing import Union
from typing import Optional
from dataclasses import dataclass


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

    def get_quality_scores(self) -> pd.DataFrame:
        """Returns scores in a per sequence bases"""

        # converting to np.array
        all_scores = []
        for entry in self.__loader():
            phred_scores = list(entry.scores)
            scores = [ord(phred_score) - 33 for phred_score in phred_scores]
            all_scores.append(scores)

        all_scores = pd.DataFrame(data=np.array(all_scores))
        return all_scores

    def get_average_score(self) -> pd.Series:
        """Returns average score of all sequences. Returns a a pd.Series object"""
        # get all scores df and take the average per column basis
        scores_df = self.get_quality_scores()
        average_score = scores_df.mean()
        return average_score

    def iter_reads(self):
        """Returns a python generator containing FastqEntries"""
        return self.__loader()

    def to_list(self):
        """Saves all lists into memory. Warning, large file sizes will use more
        memory but increase iteration performance.
        """
        return [read for read in self.__loader()]

    # ----------------------------------------
    # private functions: users do not interact with this
    # ----------------------------------------
    def __loader(self):
        """Creates a generator object that contains sequence read data as
        FastqEntry object.

        Parameters
        ----------
        """

        # iterate all row contents in fastq file
        with open(self.fpath, "r") as fastq_file:

            contents_chunk = []
            for idx, row_entry in enumerate(fastq_file):

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
