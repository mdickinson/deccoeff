from deccoeff import _Decimal, Deccoeff
from decimal import Decimal

D = Deccoeff
x = D(3) + D(-4)
print(x)

class MyDecimal(Decimal): pass

x = MyDecimal('123')
