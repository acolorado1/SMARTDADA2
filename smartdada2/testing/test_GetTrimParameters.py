import unittest
import warnings
import smartdada2.testing.GetTrimParameters as GetTrimParameters


class MyTestCase(unittest.TestCase):
    def test_trim_ends_less_than_threshold(self):
        known_list = [15, 20, 22, 30, 30, 40, 30, 30, 25, 20, 15]
        result = GetTrimParameters.trim_ends_less_than_threshold(known_list)
        self.assertEqual(result, (1, 10))

        # check raised errors
        string_list = ['30', '35', '36']
        self.assertRaises(TypeError,
                          GetTrimParameters.trim_ends_less_than_threshold,
                          string_list)

        impossible_values = [20, 50, 2]
        self.assertRaises(ValueError,
                          GetTrimParameters.trim_ends_less_than_threshold,
                          impossible_values)

        # check return values for short return sizes
        trim_L_short = [15, 18, 18, 30, 30, 40, 30, 30, 25, 20, 15]
        trim_R_short = [15, 20, 22, 30, 30, 40, 30, 30, 18, 18, 15]

        with warnings.catch_warnings(record=True) as w:
            GetTrimParameters.trim_ends_less_than_threshold(trim_L_short)
            GetTrimParameters.trim_ends_less_than_threshold(trim_R_short)
            print(w[0].category)
            assert issubclass(w[0].category, UserWarning)
            assert issubclass(w[1].category, UserWarning)


if __name__ == '__main__':
    unittest.main()
