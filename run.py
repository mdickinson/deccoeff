from deccoeff import *
from decimal import *
D = Deccoeff

x = 228276521107663040363579649908999403271918866928960310856180122057271452991526637158444248473393759326399284721107971662136571255380685903109012621859563582388276575305208
y = 998796156233724397059690401224054249226343226019830987451098900118970645644339408309927492930375453732044542023199544202935189362859382717779953303760505818737894532587115

print(D(x*y))
print(D(x) * D(y))

assert D(x*y) == D(x) * D(y)


from random import randrange

B = 10**LIMB_DIGITS

def test(m, n):
    for i in range(1000):
        if m > 0:
            x = randrange(B**(m-1), B**m)
        else:
            x = 0
        if n > 0:
            y = randrange(B**(n-1), B**n)
        else:
            y = 0
        if not D(x*y) == D(x) * D(y):
            print("Failed equality: x = %s, y = %s" % (x, y))
            assert D(x*y) == D(x) * D(y)

for m in range(0, 20):
    for n in range(0, 20):
        test(m, n)
