from cProfile import run
from test_decimal import test_main
run("test_main(arith=True)", sort=1)
