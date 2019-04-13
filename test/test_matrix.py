import unittest
from bioseq import Matrix


class TestMatrix(unittest.TestCase):
    def test_constructor_and_private_methods(self):
        m = Matrix(10, 10)
        self.assertEqual(10, len(m))
        self.assertEqual(10, len(m[0]))
        for r in m:
            self.assertTrue(all(c == 0 for c in r))
        m = Matrix(10, 10, True)
        for r in m:
            self.assertTrue(all(r))
        m = Matrix(10, 10, False)
        for r in m:
            self.assertFalse(any(r))
        self.assertEqual(len(m[3:10:2]), 4)
        self.assertEqual(len(m[3:10]), 7)
        self.assertRaises(Exception, Matrix, -1, 10)
        self.assertRaises(Exception, Matrix, 10, -1)
        self.assertRaises(Exception, Matrix, 0, 5)
        self.assertRaises(Exception, Matrix, 5, 0)
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
        self.assertEqual(5, m.max())
        m[0][7] = 5005
        self.assertEqual(5005, m.max())
        self.assertEqual(5, m.min())
        m[7][0] = -5005
        self.assertEqual(-5005, m.min())

    def test_square(self):
        m = Matrix(10, 10, 3)
        self.assertTrue(m.square())
        m = Matrix(10, 5, 3)
        self.assertFalse(m.square())


if __name__ == '__main__':
    unittest.main()