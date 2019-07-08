# Fcup ABI
Python Module for biological sequence handling developed in the course: Algorithms for Bionformatics

[![Build Status](https://travis-ci.org/msramalho/fcup-abi.svg?branch=master)](https://travis-ci.org/msramalho/fcup-abi)
[![Coverage Status](https://coveralls.io/repos/github/msramalho/fcup-abi/badge.svg)](https://coveralls.io/github/msramalho/fcup-abi)

 * The 1st delivery report can be found [here](report/part1/report.pdf).
 * The 2nd delivery report can be found [here](report/part2/report.pdf).
 * The 3rd delivery report can be found [here](report/part3/report.pdf).

## Install
`pip install git+https://github.com/msramalho/fcup-abi`

## Import

To use it simply do `import bioseq` or `from bioseq import DNASeq` (for instance).


## Demo
You can run any of the _demos_ in [demo](demo/) with `python DEMO_NAME.py` to see some examples of what this tool does and you also have a [Jupyter Notebook](README.ipynb) demonstrating the 3rd part. 

<p align="center"><img alt="some prints" src="https://i.imgur.com/HC5kmRL.gif"/></p>


## Tests
There are tests that guarantee 100% coverage (see [html version of coverage report here](htmlcov/index.html)) and there is human friendly documentation [here](docs/_build/html/index.html).

* Run tests: `python -m unittest discover`
* Run with [coverage](https://coverage.readthedocs.io/): `coverage run -m unittest discover`
* Generate (HTML) report: `coverage report` | `coverage html`

## Docs
To generate the docs again, `cd ./docs` and:
 * Windows: `make.bat html`
 * Linux: `make html`
