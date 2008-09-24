# Timings for multiplication

import timeit
from random import randrange
from deccoeff import Deccoeff as D

def getvals(N):
    a = randrange(10**N)
    b = randrange(10**N)
    return a, b

# test multiplication of two random n-digit numbers

N = 4000

T = timeit.Timer('c = a*b',
                 'from __main__ import getvals, D\n'
                 'a, b = getvals(%s)\n'
                 'a, b = D(a), D(b)' % N)

T2 = timeit.Timer('c = a*b',
                  'from __main__ import getvals\n'
                  'a, b = getvals(%s)' % N)

count = 10000
times = T.repeat(repeat = 5, number = count)
print("%s-digit by %s-digit multiply, decimal: %.9f" % (N, N, min(times)/count))
times = T2.repeat(repeat = 5, number = count)
print("%s-digit by %s-digit multiply, binary:  %.9f" % (N, N, min(times)/count))
