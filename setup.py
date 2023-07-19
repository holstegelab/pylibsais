#! /usr/bin/env python
from setuptools import setup, Extension
import numpy as np

numpy_nodepr_api = dict(define_macros=[("NPY_NO_DEPRECATED_API",
                                        "NPY_1_9_API_VERSION")])

pylibsais = Extension(
    'pylibsais',
    sources=['libsais.c', 'pylibsais.c'], **numpy_nodepr_api
)

setup(
    name='pylibsais',
    version='0.1',
    description='A Python module wrapper for the SA-IS algorithm.',
    author='Marc Hulsman',
    author_email='m.hulsman@tudelft.nl',
    long_description='A Python module wrapper for the SA-IS algorithm implementation by Ilya Grebnov.',
    include_dirs=[np.get_include() + '/numpy', 'tests'],
    ext_modules=[pylibsais],
)
