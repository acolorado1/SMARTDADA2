import unittest

import smartdada2.GetPairedendFiles as GPF


class MyTestCase(unittest.TestCase):
    def test_get_paired_files(self):
        # make sure it does not run with this data
        self.assertRaises(
            TypeError, GPF.getPairedFiles, "./test_data/SRR1591840_tunc.fastq", False
        )

        # TODO
        # run tests to make sure it works on proper test data


if __name__ == "__main__":
    unittest.main()
