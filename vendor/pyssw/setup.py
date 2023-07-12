#! /usr/bin/env python
# -*- encoding:utf-8 -*-
import os
import sys

import codecs

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from ssw.version import __version__


def read(infile):
    return codecs.open(os.path.join(os.path.dirname(__file__), infile)).read()


setup(
    name='pyssw',
    version=__version__,
    url='https://github.com/Kevinzjy/pyssw',
    description='Python binding for stripped smith waterman library',
    long_description=read('README.md'),
    long_description_content_type='text/markdown',
    author='Jinyang Zhang',
    author_email='zhangjinyang@biols.ac.cn',
    maintainer='Jinyang Zhang',
    maintainer_email='zhangjinyang@biols.ac.cn',
    license='MIT',
    keywords='circRNA',
    packages=find_packages(exclude=['doc', 'tests']),
    # entry_points={
    #     'console_scripts': [
    #         'CIRI-long=CIRI.main:main',
    #     ]
    # },
    include_package_data=True,
    zip_safe=False,
    test_suite="nose.collector",
    tests_require=['nose==1.3.7'],
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Intended Audience :: Developers',
        'Operating System :: OS Independent',
    ],
    cmdclass={
        'build_ext': build_ext,
    }
)