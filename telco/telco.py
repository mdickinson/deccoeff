""" Telco Benchmark for measuring the performance of decimal calculations

http://www2.hursley.ibm.com/decimal/telco.html
http://www2.hursley.ibm.com/decimal/telcoSpec.html

A call type indicator, c, is set from the bottom (least significant) bit of the duration (hence c is 0 or 1).
A r, r, is determined from the call type. Those calls with c=0 have a low r: 0.0013; the remainder (‘distance calls’) have a ‘premium’ r: 0.00894. (The rates are, very roughly, in Euros or dollarates per second.)
A price, p, for the call is then calculated (p=r*n). This is rounded to exactly 2 fractional digits using round-half-even (Banker’s round to nearest).
A basic tax, b, is calculated: b=p*0.0675 (6.75%). This is truncated to exactly 2 fractional digits (round-down), and the total basic tax variable is then incremented (sumB=sumB+b).
For distance calls: a distance tax, d, is calculated: d=p*0.0341 (3.41%). This is truncated to exactly 2 fractional digits (round-down), and then the total distance tax variable is incremented (sumD=sumD+d).
The total price, t, is calculated (t=p+b, and, if a distance call, t=t+d).
The total prices variable is incremented (sumT=sumT+t).
The total price, t, is converted to a string, s.

"""

from struct import unpack
from time import clock as time
from decimal import *
import sys, os

# To run the full test with 1,000,000 entries:  python telco.py full
test = 'full' not in ' '.join(sys.argv[1:]).lower()

if test:
    filename = "telco.testb"
    expected = map(Decimal, "8.91 0.50 0.22".split())
    print("  Time Rate | Price   Btax   Dtax  | Output")
    print("------------+----------------------+--------")
else:
    filename = "expon180.1e6b"
    if not os.access(filename, os.F_OK):
        print("You must download and unzip the test file from: http://www2.hursley.ibm.com/decimal/expon180-1e6b.zip")
        sys.exit(-1)
    expected = map(Decimal, "1004737.58 57628.30 25042.17".split())

getcontext().rounding = ROUND_DOWN
rates = list(map(Decimal, ('0.0013', '0.00894')))
twodig = Decimal('0.01')
Banker = Context(rounding=ROUND_HALF_EVEN)
basictax = Decimal("0.0675")
disttax = Decimal("0.0341")

infil = open(filename, "rb")
outfil = open("telco.out", "w")
start = time()

sumT = Decimal("0")   # sum of total prices
sumB = Decimal("0")   # sum of basic tax
sumD = Decimal("0")   # sum of 'distance' tax

while 1:
    datum = infil.read(8)
    if datum == b'': break
    n, =  unpack('>Q', datum)

    calltype = n & 1
    r = rates[calltype]

    p = Banker.quantize(r * n, twodig)

    b = p * basictax
    b = b.quantize(twodig)
    sumB += b

    t = p + b

    if calltype:
        d = p * disttax
        d = d.quantize(twodig)
        sumD += d
        t += d

    sumT += t
    print(t, file=outfil)

    if test:
        print('%6d   %1s  |%6s %6s %6s  |%6s' % (n, 'LD'[calltype], p, b, (not calltype and " " or d), t))

infil.close()
outfil.close()
end = time()

print('\nControl totals:')
print('Actual  ', list(map(str, (sumT, sumB, sumD))))
print('Expected', list(map(str, expected)))
print('Elapsed time:', end-start)

