.PHONY: help all lib install test clean

all: lib install test

help:
	@echo "make all"
	@echo "       install pyssw and dependencies into current python environment"
	@echo "make lib"
	@echo "       prepare development environment, for debug only"
	@echo "make install"
	@echo "       install pyssw"
	@echo "make test"
	@echo "       unit test, run after installation"
	@echo "make clean"
	@echo "       clean python cache files"

lib:
	cd striped_smith_waterman && make libssw.so
	cp striped_smith_waterman/libssw.so ssw/libssw.so

install:
	python3 setup.py install

test:
	python3 setup.py build_ext --inplace
	python3 setup.py test

clean:
	python3 setup.py clean --all
