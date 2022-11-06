"""
Testing module tests all random, positive, and negative cases that FastqReader
may encounter
"""
import os
import unittest

# smartdada2 imports
from smartdada2.testing.help_test_funcs import toy_sequencer
from smartdada2.reader.reader import FastqReader, FastqEntry


class TestFastqEntry(unittest.TestCase):
    """FastqEntry is a datatype that is used within FastqReader. Each entry
    that is read within a FastqFile is stored in a FastqEntry data type.

    This datatype is important for all the methods to work within the
    FastqReader
    """

    # ------------------------------
    # Setup class methods
    # -- setups up files
    # -- removes files
    # ------------------------------
    @classmethod
    def set_up(cls) -> None:
        """Sets up files positive, negative, and random"""

    @classmethod
    def tear_down(cls) -> None:
        """removes all files"""
        os.remove(cls.upper_fastq)
        os.remove(cls.lower_fastq)
        pass

    # ------------------------------
    # Unit Testing
    # ------------------------------
    def test_entry(self):
        """Tests a positive case of entry"""

        # set entries
        entries = [
            [
                "test_fastq.test.1 1 length=50",
                "KDUHDNKSBARKSMKVNRBKWYRWDVBDYDRCUDYDVKANAAKKGBSBKG",
                "test_fastq.test.1 1 length=50",
                "(#A!.&+,2@\"$&#&2('-231%-!+!$,!,;6%)',(>#;>20:\"$%#E",
            ],
            [
                "test_fastq.test.2 1 length=50",
                "KWTBBNNDAVHWMYRRSBCRDSNGMDYYNSKTKNBVRDMUSHNRRSWDHK",
                "test_fastq.test.2 1 length=50",
                "&\"@*%3\"&D((3'*!2#2D%%'#<&!2.,2&A)0,')+&<.'1+*--!+1",
            ],
        ]
        pass

    def test_entry_lowercases(self):
        """Tests if fastq files contains lower case sequences"""

        # set entries with lower case sequences
        raw_entries = [
            [
                "@test_fastq.test.1 1 length=50",
                "kduhdnksbarksmkvnrbkwyrwdvbdydrcudydvkanaakkgbsbkg",
                "@test_fastq.test.1 1 length=50",
                "(#A!.&+,2@\"$&#&2('-231%-!+!$,!,;6%)',(>#;>20:\"$%#E",
            ],
            [
                "@test_fastq.test.2 1 length=50",
                "kwtbbnndavhwmyrrsbcrdsngmdyynsktknbvrdmushnrrswdhk",
                "@test_fastq.test.2 1 length=50",
                "&\"@*%3\"&D((3'*!2#2D%%'#<&!2.,2&A)0,')+&<.'1+*--!+1",
            ],
        ]

        # create entry classes
        raw_reads_1 = raw_entries[0]
        test_entry_1 = FastqEntry(
            header=raw_reads_1[0],
            seq=raw_reads_1[1],
            scores=raw_reads_1[3],
            length=len(raw_reads_1[1]),
        )

        raw_reads_2 = raw_entries[1]
        test_entry_2 = FastqEntry(
            header=raw_reads_2[0],
            seq=raw_reads_2[1],
            scores=raw_reads_2[3],
            length=len(raw_reads_2[1]),
        )

        # setting expected variables
        expected_header_1 = "@test_fastq.test.1 1 length=50"
        expected_seq_1 = "KDUHDNKSBARKSMKVNRBKWYRWDVBDYDRCUDYDVKANAAKKGBSBKG"
        expected_score_1 = "(#A!.&+,2@\"$&#&2('-231%-!+!$,!,;6%)',(>#;>20:\"$%#E"
        expected_length_1 = 50
        expected_r_seq_1 = False

        expected_header_2 = "@test_fastq.test.2 1 length=50"
        expected_seq_2 = "KWTBBNNDAVHWMYRRSBCRDSNGMDYYNSKTKNBVRDMUSHNRRSWDHK"
        expected_score_2 = "&\"@*%3\"&D((3'*!2#2D%%'#<&!2.,2&A)0,')+&<.'1+*--!+1"
        expected_length_2= 50
        expected_r_seq_2 = True

        # setting tests
        # -- testing first entry
        self.assertEqual(expected_header_1, test_entry_1.header)
        self.assertEqual(expected_seq_1, test_entry_1.seq)
        self.assertEqual(expected_score_1, test_entry_1.scores)
        self.assertEqual(expected_length_1, test_entry_1.length)
        self.assertEqual(expected_r_seq_1, test_entry_1.rseq)

        # -- testing second entry
        self.assertEqual(expected_header_2, test_entry_2.header)
        self.assertEqual(expected_seq_2, test_entry_2.seq)
        self.assertEqual(expected_score_2, test_entry_2.scores)
        self.assertEqual(expected_length_2, test_entry_2.length)
        self.assertEqual(expected_r_seq_2, test_entry_2.rseq)


    def test_amplicon_direction(self):
        """Tests if can detect amplicon direction"""
        pass


class TestFastqReader(unittest.TestCase):
    """FastqReader is the main engine that attempts to extract as much
    information from fastq sequences."""

    # ------------------------------
    # Setup class methods
    # -- setups up files
    # -- removes files
    # ------------------------------
    @classmethod
    def set_up(cls) -> None:
        """Sets up files positive, negative, and random"""

        # setting vars
        cls.upper_fastq = "upper_seq.fastq"
        cls.lower_fastq = "lower_seq.fastq"

        # generating regular fastq fule
        with open(cls.upper_fastq, "w") as f:
            reads = toy_sequencer(200, 30, rev_seq=True, seed=42)
            for read in reads:
                for read_data in read:
                    f.write(f"{read_data}\n")

        # generating sequence with lower case sequence
        with open(cls.lower_fastq, "w") as f:
            reads = toy_sequencer(200, 30, rev_seq=True, seed=42)
            for read in reads:
                read[1] = read[1].lower()
                for read_data in read:
                    f.write(f"{read_data}\n")

    @classmethod
    def tear_down(cls) -> None:
        """removes all files"""
        pass

    # ------------------------------
    # Unit testing
    # ------------------------------
    def test_reader_instance_1(self) -> None:
        """Positive case of a successful entry"""
        pass

    def test_reader_instance_2(self) -> None:
        """Tests if can create fastq files with lower case sequences"""
        pass

    def test_reader_instance_3(self) -> None:
        """Tests if Fastq files contains unwanted types. For example, if digits
        or booleans or unknown characters are captured within the sequence"""
        pass

    def test_reader_instance_4(self) -> None:
        pass
