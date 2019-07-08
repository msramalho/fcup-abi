from setuptools import setup
import os

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(name='bioseq',
      version='1.0',
      description='Biological Sequences Handling',
      url='https://github.com/msramalho/fcup-abi',
      author='msramalho',
      packages=['bioseq'],
      install_requires=required,
      zip_safe=False)
