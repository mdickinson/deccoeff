import decimal
D = decimal.Decimal
DI = decimal.Deccoeff

# test divmod (and also, incidentally, multiplication and addition)
from random import randrange, choice

def test_divmod():
    # random lengths for deccoeffs
    alen = randrange(10)
    blen = randrange(10)

    a = DI(''.join(choice('0123456789') for i in range(alen)))
    b = DI(''.join(choice('0123456789') for i in range(blen)))
    if (not b):
        return
    q, r = divmod(a, b)
    assert b*q+r == a and 0 <= r < b

    # also check that results are same as for longs
    qq, rr = divmod(int(str(a) or '0'), int(str(b) or '0'))
    assert qq == int(str(q) or '0') and rr == int(str(r) or '0')



x = DI('9000000000000000')
print(x << 3)
