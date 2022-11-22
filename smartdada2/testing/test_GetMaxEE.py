import unittest
import pandas as pd
import numpy as np
import smartdada2.GetMaxEE as GME
import smartdada2.GetTrimParameters as GTP
from smartdada2.reader import reader


class MyTestCase(unittest.TestCase):
    def test_read_size_by_maxEE(self):
        self.assertRaises(
            TypeError, GME.read_size_by_maxEE, self.fastqentry, self.EE_by_size_df
        )
        self.assertRaises(
            ValueError, GME.read_size_by_maxEE, self.fastqentry, self.rand_df
        )
        df = GME.read_size_by_maxEE(self.fastqentry, self.EE_by_size_df)
        self.assertEqual(isinstance(df, pd.DataFrame), True)


@classmethod
def setUpClass(cls):
    cls.fastqentry = reader.FastqReader("./test_data/SRR1591840_tunc.fastq")
    cls.avg_scores = cls.fastqentry.get_average_score()
    cls.avg_scores_list = cls.avg_scores["AverageQualityScore"].values.tolist()
    left, right = GTP.trim_ends_less_than_threshold(cls.avg_scores_list)
    cls.EE_by_size_df = GTP.read_size_by_avg_EE(fqe, left, right)
    cls.rand_array = np.random.randint(low=2, high=10, size=(2, 3))
    cls.rand_df = pd.DataFrame(
        cls.rand_array, columns=["A", "B", "C"], index=["a", "b"]
    )


@classmethod
def tearDownClass(cls):
    del cls.fastqentry
    del cls.rand_array
    del cls.rand_df
    del cls.avg_scores
    del cls.avg_scores_list
    del cls.EE_by_size_df


if __name__ == "__main__":
    unittest.main()
