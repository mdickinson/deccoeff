from deccoeff import Deccoeff as D, LIMB_DIGITS, MAX_DIGITS
from test.support import run_unittest
from operator import add, mul, sub, floordiv, mod, pow
import operator
import unittest
from random import randrange

opname = {
    add: '+',
    sub: '-',
    mul: '*',
    floordiv: '/',
    mod: '%',
}

class DecccoeffTest(unittest.TestCase):
    def compare_op(self, op, a, b):
        """Check that Deccoeff and Python integers give the same
        result for the given arithmetic operation."""

        if op is sub and a < b:
            self.assertRaises(OverflowError, op, D(a), D(b))
        elif op in (floordiv, mod) and b == 0:
            self.assertRaises(ZeroDivisionError, op, D(a), D(b))
        else:
            int_result = D(op(a, b))
            deccoeff_result = op(D(a), D(b))
            if int_result != deccoeff_result:
                self.fail('integer and deccoeff give different results for: '
                          '%s %s %s' % (a, opname[op], b))

    def testRandomOps(self):
        sizes = [10**i for i in
                   [0, 1, 2, 4, 5, 9, 10, 18, 19, 20, 50, 100]]
        trials = 100
        for op in opname.keys():
            for a_max in sizes:
                for b_max in sizes:
                    for n in range(trials):
                        self.compare_op(op, randrange(a_max), randrange(b_max))

    def testAdd(self):
        # additional tests for addition
        for i in range(100):
            for j in range(100):
                self.compare_op(add, i, j)
        test_values = [
            (0, 123456789),
            (0, 999999999),
            (0, 10**9),
            (1, 10**9-1),
            (1, 10**9),
            (1, 10**18-1),
            (123456789, 876543210),
            (123456789, 876543211),
            (123456789, 999999999876543210),
            (123456789, 999999999876543211),
            (123456789, 999999999876543212),
            ]
        for i, j in test_values:
            self.compare_op(add, i, j)

    def testMixed(self):
        class MyInteger(object):
            def __init__(self, n):
                self.n = n
            def __index__(self):
                return self.n

        # integer types should be automatically converted to Deccoeffs
        self.assertEqual(D(123) + 47, D(123+47))
        self.assertEqual(True * D(123), D(True * 123))
        self.assertEqual(MyInteger(12345) % D(123), D(12345 % 123))

        # but an attempt to convert a negative integer type
        # should cause an error, even if the operation result is nonnegative
        self.assertRaises(OverflowError, add, D(47), -12)
        self.assertRaises(OverflowError, add, -12, D(47))

def test_main():
    run_unittest(DecccoeffTest)

if __name__ == '__main__':
    test_main()
