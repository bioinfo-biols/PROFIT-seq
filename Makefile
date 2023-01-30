.PHONY: help all lib

all: lib

help:
	@echo "make all"
	@echo "       install pyssw and dependencies into current python environment"
	@echo "make lib"
	@echo "       prepare development environment, for debug only"

lib:
    pip install -r scripts/requirements.txt
    cd vendor/pyssw & make
    cd vendor/Porechop && python setup.py install

