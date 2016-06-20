from setuptools import setup, Extension

from Cython.Build import cythonize
import numpy as np

extensions = cythonize([
    Extension("DREAM11.lib", 
              ["DREAM11/lib.pyx", ],
              include_dirs=[np.get_include()]),
])

config = {
    'include_package_data': True,
    'ext_modules': extensions, 
    'description': 'DREAM11',
    'author': 'Nathan Boley',
    'url': 'NA',
    'download_url': 'http://github.com/nboley/DREAM_tfbinding_challenge/',
    'author_email': 'npboley@gmail.com',
    'version': '0.1.1',
    'packages': ['DREAM11', ],
    'setup_requires': [],
    'install_requires': [ 'scipy', 'numpy' ],
    'scripts': [],
    'name': 'DREAM11'
}

if __name__== '__main__':
    setup(**config)
