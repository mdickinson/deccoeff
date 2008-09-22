PYTHON = ../py3k/python.exe

all: deccoeff.so
module:
	$(PYTHON) setup.py build
	cp build/lib*/deccoeff.so .
clean:
	-rm -rf build/
	-rm deccoeff.o deccoeff.so
	-rm *.pyc
	-rm config.status config.log autom4te.cache

test: deccoeff.so
	$(PYTHON) test_decimal.py

run: deccoeff.so
	$(PYTHON) -i run.py

deccoeff.so: deccoeff.c
	$(PYTHON) setup.py build
	cp build/lib*/deccoeff.so .
