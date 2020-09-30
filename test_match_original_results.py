import unittest

from pathlib import Path
import numpy as np
import ACES_DataSetCalc_v031 as sy
from tty_ui import parse_cct, parse_tint

USER_ILLUMINANT = [0.29902, 0.31485]
CCT = 5400
D_UV = 1
COLORSPACE = 'AP0'


class MyTestCase(unittest.TestCase):

    def test_identical_output(self):
        reference_output = np.genfromtxt('data/reference_output.csv', delimiter=',', usecols = [1, 2, 3])
        target_cap_xyz_values = sy.calc_target_cap_xyz_values(USER_ILLUMINANT, CCT, D_UV, sy.macbeth, COLORSPACE)
        csv_pathname = Path('/tmp', sy.unique_filename((USER_ILLUMINANT[0], USER_ILLUMINANT[1]),
                                                       CCT, D_UV, COLORSPACE))
        sy.write_csv(str(csv_pathname), target_cap_xyz_values)
        test_output = np.genfromtxt(str(csv_pathname), delimiter=',', usecols=[1, 2, 3])
        self.assertIsNone(np.testing.assert_array_equal(test_output, reference_output))

    def test_cct_parser(self):
        self.assertEqual(5600.0, parse_cct('5600'))
        self.assertEqual(5600.0, parse_cct('  5600 K  '))
        self.assertEqual(5600.0, parse_cct('5600.'))
        self.assertEqual(5600.0, parse_cct('5600.0'))
        self.assertEqual(5600.1, parse_cct('5600.1'))
        self.assertEqual(5600.0, parse_cct('5600K'))
        self.assertEqual(5600.0, parse_cct('5600.K'))
        self.assertEqual(5600.0, parse_cct('5600.0K'))
        self.assertEqual(5600.1, parse_cct('5600.1K'))
        with self.assertRaises(SyntaxError):
            cct = parse_cct('K5600')
        with self.assertRaises(SyntaxError):
            cct = parse_cct('5600-')
        with self.assertRaises(SyntaxError):
            cct = parse_cct('-5600')

    def test_tint_parser(self):
        self.assertEqual(5.0, parse_tint('5'))
        self.assertEqual(5.0, parse_tint('  5  '))
        self.assertEqual(5.0, parse_tint('5.'))
        self.assertEqual(5.0, parse_tint('5.0'))
        self.assertEqual(5.1, parse_tint('5.1'))
        self.assertEqual(-5.0, parse_tint('-5'))
        self.assertEqual(-5.0, parse_tint('-5.'))
        self.assertEqual(-5.0, parse_tint('-5.0'))
        self.assertEqual(-5.1, parse_tint('-5.1'))
        self.assertEqual(5.0, parse_tint('+5'))
        self.assertEqual(5.0, parse_tint('+5.'))
        self.assertEqual(5.0, parse_tint('+5.0'))
        self.assertEqual(5.1, parse_tint('+5.1'))
        with self.assertRaises(SyntaxError):
            cct = parse_tint('-5-')
        with self.assertRaises(SyntaxError):
            cct = parse_tint('5+')


if __name__ == '__main__':
    unittest.main()
