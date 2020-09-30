import unittest

from pathlib import Path
import numpy as np
import ACES_DataSetCalc_v031 as sy

ILLUMINANT = (0.29902, 0.31485)
CCT = 5400
TINT = 1
COLORSPACE = 'AP0'


class MyTestCase(unittest.TestCase):

    def test_identical_output(self):
        reference_output = np.genfromtxt('data/reference_output.csv', delimiter=',', usecols = [1, 2, 3])
        dataSet = sy.calc_data_set(ILLUMINANT, CCT, TINT, sy.macbeth, COLORSPACE)
        csv_pathname = Path('/tmp', sy.unique_filename((ILLUMINANT[0], ILLUMINANT[1]),
                                                       CCT, TINT, COLORSPACE))
        sy.output_to_csv_file(str(csv_pathname), dataSet)
        test_output = np.genfromtxt(str(csv_pathname), delimiter=',', usecols=[1, 2, 3])
        self.assertIsNone(np.testing.assert_array_equal(test_output, reference_output))


if __name__ == '__main__':
    unittest.main()
