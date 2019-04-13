# Fcup ABI
Algorithms for Bionformatics

[![Bioseq Module Build Status](https://travis-ci.com/msramalho/fcup-abi.svg?token=peexEWPZnZuwQpsysvdo&branch=master)](https://travis-ci.com/msramalho/fcup-abi)

This project contains a python module that is capable of handling biological sequences and perform operations on them, with special focus on DNA, RNA and Protein sequences.

 * The delivery 1 report can be found [here](report/report_1.pdf).

This file should be seen with a Markdown reader _or a human that can read markdown :)_.

To use it simply do `import bioseq` or `from bioseq import DNASeq` (for instance).

You can run `python run_me.py` to see some examples of what this tool does. 

Additionally, there are tests that guarantee 100% coverage (see [html version of coverage report here](htmlcov/index.html)) and there is human friendly documentation [here](docs/_build/html/index.html).

## Tests
* Run tests: `python -m unittest discover`
* Run with [coverage](https://coverage.readthedocs.io/): `coverage run -m unittest discover`
* Generate (HTML) report: `coverage report` | `coverage html`

![](https://i.imgur.com/uI7VdLS.png)

## Docs
To generate the docs again, `cd ./docs` and:
 * Windows: `make.bat html`
 * Linux: `make html`
