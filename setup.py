from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy as np

setup(ext_modules=cythonize("lib.pyx"), include_dirs=[np.get_include()])
