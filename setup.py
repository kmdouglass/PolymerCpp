# © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE,
# Switzerland, Laboratory of Experimental Biophysics, 2017
# See the LICENSE.txt file for more details.

from distutils.core import setup, Extension
from os import listdir

# Find all source files
source_ext = '.cpp'
source_dir = 'PolymerCpp/core/'
source_files = ['{0:s}{1:s}'.format(source_dir, f)
                for f in listdir(source_dir)
                if f.endswith(source_ext)]

module = Extension('PolymerCppCore',
                    define_macros = [('MAJOR_VERSION', '0'),
                                     ('MINOR_VERSION', '0')],
                    include_dirs = ['include/'],
                    sources = source_files,
                    extra_compile_args = ['-std=c++11','-O2', '-fPIC'])

config={'name': 'PolymerCpp',
        'version': '0.0.1',
        'description': 'A 3D wormlike chain algorithm written in C++',
        'author': 'Marcel Stefko, Kyle M. Douglass',
        'author_email': 'marcel.stefko@epfl.ch, kyle.m.douglass@gmail.com',
        'url': 'https://github.com/MStefko/PolymerCpp',
        'packages': ['PolymerCpp'],
        'ext_modules': [module]}

setup(**config)
