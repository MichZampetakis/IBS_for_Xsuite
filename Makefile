# Copyright 2019 Felix Soubelet <felix.soubelet@cern.ch>
# MIT License

# Documentation for most of what you will see here can be found at the following links:
# for the GNU make special targets: https://www.gnu.org/software/make/manual/html_node/Special-Targets.html
# for python packaging: https://docs.python.org/3/distutils/introduction.html

# ANSI escape sequences for colors
# In order, B=bold, C=cyan, D=dark blue, E=end, P=pink, R=red and Y=yellow
B=\033[1m
C=\033[96m
D=\033[34m
E=\033[0m
P=\033[95m
R=\033[31m
Y=\033[33m

.PHONY : help build clean docs format install lines lint typing tests

all: install

help:
	@echo "Please use 'make $(R)<target>$(E)' where $(R)<target>$(E) is one of:"
	@echo "  $(R) build $(E)  \t  to build wheel and source distribution with $(P)Hatch$(E)."
	@echo "  $(R) clean $(E)  \t  to recursively remove build, run and bitecode files/dirs."
	@echo "  $(R) docs $(E)  \t  to build the documentation for the package with $(P)Sphinx$(E)."
	@echo "  $(R) format $(E)  \t  to recursively apply PEP8 formatting through the $(P)Black$(E) and $(P)isort$(E) cli tools."
	@echo "  $(R) install $(E)  \t  to $(C)pip install$(E) this package into the current environment."
	@echo "  $(R) lines $(E)  \t  to count lines of code in the package folder with the $(P)tokei$(E) tool."
	@echo "  $(R) lint $(E)  \t  to lint the packages' code though $(P)Ruff$(E)."
	@echo "  $(R) typing $(E)  \t  to run type checking on the codebase with $(P)MyPy$(E)."
	@echo "  $(R) tests $(E)  \t  to run the test suite with $(P)pytest$(E)."


# ----- Dev Tools Targets ----- #

build:
	@echo "Re-building wheel and sdist"
	@echo "Cleaning up package builds and distutils remains."
	@find . -type d -name "*build" -exec rm -rf {} +
	@find . -type d -name "*dist" -exec rm -rf {} +
	@hatch build --clean
	@echo "Created build is located in the $(C)dist$(E) folder."

clean:
	@echo "Cleaning up documentation pages."
	@rm -rf doc_build
	@echo "Cleaning up sphinx-gallery build artifacts."
	@rm -rf docs/gallery
	@rm -rf docs/gen_modules
	@echo "Cleaning up package builds and distutils remains."
	@find . -type d -name "*build" -exec rm -rf {} +
	@find . -type d -name "*dist" -exec rm -rf {} +
	@rm -rf xibs.egg-info
	@rm -rf .eggs
	@echo "Cleaning up bitecode files and python cache."
	@find . -type f -name '*.py[co]' -delete -o -type d -name __pycache__ -delete
	@echo "Cleaning up pytest cache & test artifacts."
	@find . -type d -name '*.pytest_cache' -exec rm -rf {} + -o -type f -name '*.pytest_cache' -exec rm -rf {} +
	@echo "Cleaning up mypy and ruff caches."
	@find . -type d -name "*.mypy_cache" -exec rm -rf {} +
	@find . -type d -name "*.ruff_cache" -exec rm -rf {} +
	@echo "Cleaning up ipython notebooks cache."
	@find . -type d -name "*.ipynb_checkpoints" -exec rm -rf {} +
	@echo "Cleaning up coverage reports."
	@find . -type f -name '.coverage*' -exec rm -rf {} + -o -type f -name 'coverage.xml' -delete
	@echo "All cleaned up!\n"

docs:
	@echo "Building static pages with $(D)Sphinx$(E)."
	@python -m sphinx -v -b html docs doc_build -d doc_build

format:
	@echo "Formatting code to PEP8 with $(P)isort$(E) and $(P)Black$(E). Max line length is 110 characters."
	@python -m isort . && black .

install: clean
	@echo "Installing with $(D)pip$(E) in the current environment."
	@python -m pip install . -v

lines: format
	@tokei xibs

lint: format
	@echo "Linting code with $(P)Pylint$(E)."
	@ruff check xibs/

typing: format
	@echo "Checking code typing with $(P)mypy$(E)."
	@python -m mypy xibs


# ----- Tests Targets ----- #

tests:  # all tests not involving pyhdtoolkit.cpymadtools
	@python -m pytest -v

# Catch-all unknow targets without returning an error. This is a POSIX-compliant syntax.
.DEFAULT:
	@echo "Make caught an invalid target! See help output below for available targets."
	@make help
