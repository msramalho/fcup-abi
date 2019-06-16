import io
import os
import unittest
import unittest.mock
from bioseq import Matrix


class TestMatrix(unittest.TestCase):
    def test_constructor_and_private_methods(self):
        m = Matrix(10, 10)
        self.assertEqual(10, len(m))
        self.assertEqual(10, len(m[0]))
        for r in m: self.assertTrue(all(c == 0 for c in r))
        m = Matrix(10, 10, True)
        for r in m: self.assertTrue(all(r))
        m = Matrix(10, 10, False)
        for r in m: self.assertFalse(any(r))
        self.assertEqual(len(m[3:10:2]), 4)
        self.assertEqual(len(m[3:10]), 7)
        # setitem
        m[0] = list(range(10))
        self.assertEqual(m[0], list(range(10)))
        self.assertRaises(Exception, m.__setitem__, 0, [1])
        # set_col
        m.set_col(0, list(range(10)))
        self.assertEqual(m[3][0], 3)
        self.assertRaises(Exception, m.set_col, 0, [1])
        # invalid matrix sizes
        self.assertRaises(Exception, Matrix, -1, 10)
        self.assertRaises(Exception, Matrix, 10, -1)
        self.assertRaises(Exception, Matrix, 0, 5)
        self.assertRaises(Exception, Matrix, 0, 0)
        self.assertRaises(Exception, Matrix, -10, -10)

    def test_sum(self):
        m = Matrix(10, 10, 3)
        self.assertEqual(300, m.sum())

    def test_add_val(self):
        m = Matrix(10, 10, 1)
        self.assertEqual(100, m.sum())
        m.add_val(5)
        self.assertEqual(600, m.sum())

    def test_mull_val(self):
        m = Matrix(10, 10, 1)
        self.assertEqual(100, m.sum())
        m.mul_val(5)
        self.assertEqual(500, m.sum())

    def test_max_min(self):
        m = Matrix(10, 10, 5)
        self.assertEqual(100, len(m.max()))
        self.assertEqual((5, 0, 0), m.max()[0])
        m[0][7] = 5005
        self.assertEqual([(5005, 0, 7)], m.max())
        self.assertEqual(99, len(m.min()))
        self.assertEqual((5, 0, 0), m.min()[0])
        m[7][0] = -5005
        self.assertEqual([(-5005, 7, 0)], m.min())

    def test_square(self):
        m = Matrix(10, 10, 3)
        self.assertTrue(m.square())
        m = Matrix(10, 5, 3)
        self.assertFalse(m.square())
        m = Matrix(10)
        self.assertTrue(m.square())

    def test_sum(self):
        m = Matrix(10, 10, 3)
        self.assertEqual(3, m.last())
        m[9][9] = 15
        self.assertEqual(15, m.last())

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def _test_str(self, m,l, mock_stdout):
        print(m)
        output = mock_stdout.getvalue()
        self.assertEqual(len(output), l)
        self.assertIn("12", output)

    def test_str(self):
        m = Matrix(10, 4, 12)
        self._test_str(m, 131)
        m.rows = ["-"]*10
        m.cols = ["-"]*4
        self._test_str(m, 166)

    def test_display(self):
        filename = "temp_matrix.png"
        m = Matrix(10, 4, 1)
        m[0][0] = 0
        m.display(filename)
        os.remove(filename)


if __name__ == '__main__':
    unittest.main()
