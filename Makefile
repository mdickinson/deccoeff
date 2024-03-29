# Replace this with something that points to your local Python 3 installation.
PYTHON = python3

all: deccoeff.so
module:
	$(PYTHON) setup.py build
	cp build/lib*/deccoeff.so .
clean:
	-rm -fr build/
	-rm deccoeff.o deccoeff.so
	-rm *.pyc
	-rm config.status config.log deccoeff_config.h
	-rm -fr autom4te.cache/
	-rm telco/telco.out

distclean: clean
	-rm deccoeff_config.h.in
	-rm configure

test: deccoeff.so
	$(PYTHON) test_deccoeff.py
	$(PYTHON) test_decimal.py

run: deccoeff.so
	$(PYTHON) -i run.py

profile: deccoeff.so
	$(PYTHON) profile.py

install: all
	$(PYTHON) setup.py install

deccoeff.so: deccoeff.c
	$(PYTHON) setup.py build
	cp build/lib*/deccoeff.so .
