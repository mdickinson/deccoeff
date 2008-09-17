PYTHON = ../py3k/python.exe

all: _decimal.so
module:
	$(PYTHON) setup.py build
	cp build/lib*/_decimal.so .
clean:
	-rm -rf build/
	-rm _decimal.o _decimal.so limb.o
	-rm decimal.pyc test_decimal.pyc
test: _decimal.so
	$(PYTHON) test_decimal.py

run: _decimal.so
	$(PYTHON) -i run.py
profile: _decimal.so
	$(PYTHON) -i profile.py

debug: _decimal.so
	gdb $(PYTHON)

_decimal.so: _decimal.c limb.c
	$(PYTHON) setup.py build
	cp build/lib*/_decimal.so .
