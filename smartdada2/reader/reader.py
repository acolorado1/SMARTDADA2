import sy
from pathlib import Path
from typing import Union
from typing import Optional
from typing import Iterable
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
        reversed: Optional[bool] = False,
        technology: Optional[str] = "illumina",
    ):

        # FastqReader accessible parameters
        self.fpath: str = fpath
        self.reversed: bool = reversed
        self.technology: str = technology

        # parameters that are not set
        self.reads: Union[None, Iterable[FastqEntry]] = None

        # executing setter function when initializing FastqReader
        self.__set_reads()

    def get_quality_scores(self) -> pd.DataFrame:
        """Returns scores in a per sequence bases"""

        # converting to np.array
        all_scores = []
        for entry in self.reads:
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

    # setter function
    def __set_reads(self) -> None:
        """Reads fastq file and converts each entry into FastqEntry object,
        which stores the header, sequence

        Returns
        -------
        None
            Sets FastqReader attributes


        Raises
        ------
        FileNotFoundError
            raised if the provided fpath
        ValueError
            raised invalid files are passed
        """
        # creating path object
        fastq_path = Path(self.fpath)

        # error handling
        try:
            if not fastq_path.is_file():
                raise FileNotFoundError
            elif not fastq_path.suffix == ".fastq":
                raise ValueError
        except FileNotFoundError as e:
            e_name = f"{e.__class__.__name__}:"
            e_msg = f"File path provided is invalid"
            print(e_name, e_msg)
            sys.exit(1)
        except ValueError as e:
            e_name = f"{e.__class__.__name__}:"
            e_msg = f"Invalid file format provided. requires .fastq files"
            print(e_name, e_msg)
            sys.exit(1)

        # reading in file
        contents_chunk = []
        fastq_entries = []

        # reading the fastq file
        with open(self.fpath, "r") as fastq_file:
            for idx, line_cont in enumerate(fastq_file):
                content = line_cont.rstrip("\n")
                contents_chunk.append(content)

                # checking if there are 4 elements in the list
                # -- 4 lines = 1 entry
                if len(contents_chunk) == 4:

                    # convert into FastqEntry and store it
                    entry = FastqEntry(
                        header=contents_chunk[0],
                        seq=contents_chunk[1],
                        scores=contents_chunk[3],
                        length=len(contents_chunk[1]),
                    )
                    fastq_entries.append(entry)

                    # empty out list container
                    contents_chunk.clear()

        # setting class attributes
        self.reads = fastq_entries
