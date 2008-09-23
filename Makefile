PYTHON = python3.0

all: deccoeff.so
module:
	$(PYTHON) setup.py build
	cp build/lib*/deccoeff.so .
clean:
	-rm -fr build/
	-rm deccoeff.o deccoeff.so
	-rm *.pyc
	-rm config.status config.log config.h
	-rm -fr autom4te.cache/

distclean: clean
	-rm deccoeff_config.h.in
	-rm configure

test: deccoeff.so
	$(PYTHON) test_decimal.py

run: deccoeff.so
	$(PYTHON) -i run.py

install: all
	$(PYTHON) setup.py install

deccoeff.so: deccoeff.c
	$(PYTHON) setup.py build
	cp build/lib*/deccoeff.so .
