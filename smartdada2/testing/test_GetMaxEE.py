import unittest
import smartdada2.GetMaxEE as GME
from smartdada2.reader import reader


class MyTestCase(unittest.TestCase):
    def test_read_size_by_maxEE(self):
        # must be able to run with the testdata
        return None

if __name__ == '__main__':
    unittest.main()




fqe = reader.FastqReader(
    "/Users/burkhang/PythonProjs/DADA2ParameterExploration/data/SRR1591840_tunc.fastq"
)

# perform obvious trimming
avg_scores_df = fqe.get_average_score()
avg_scores_list = avg_scores_df["AverageQualityScore"].values.tolist()

left, right = GTP.trim_ends_less_than_threshold(avg_scores_list)

EE_by_size_df = GTP.read_size_by_avg_EE(fqe, left, right)

print(read_size_by_maxEE(fqe, EE_by_size_df) )