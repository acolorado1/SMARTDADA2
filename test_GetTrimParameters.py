import unittest
import GetTrimParameters


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
        Ls_res = GetTrimParameters.trim_ends_less_than_threshold(trim_L_short)
        self.assertEqual(Ls_res,
                         'trim left value might be too high: 3')

        trim_R_short = [15, 20, 22, 30, 30, 40, 30, 30, 18, 18, 15]
        Rs_res = GetTrimParameters.trim_ends_less_than_threshold(trim_R_short)
        self.assertEqual(Rs_res,
                         'trim right value might be too low: 8')


if __name__ == '__main__':
    unittest.main()
