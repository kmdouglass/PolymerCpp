# © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE,
# Switzerland, Laboratory of Experimental Biophysics, 2017
# See the LICENSE.txt file for more details.

from setuptools import setup, Extension
from os.path import join
from os import listdir

# Find all source files
source_ext = '.cpp'
source_dir = join('PolymerCpp', 'core')
source_files = [join(source_dir, f)
                for f in listdir(source_dir)
                if f.endswith(source_ext)]

module = Extension('PolymerCppCore',
                    define_macros = [('MAJOR_VERSION', '0'),
                                     ('MINOR_VERSION', '1')],
                    include_dirs = ['include'],
                    sources = source_files,
                    extra_compile_args = ['-std=c++11','-O2', '-fPIC'])

config={'name': 'PolymerCpp',
        'version': '0.1.3',
        'description': '2D and 3D wormlike chain generator for Python and written in C++',
        'author': 'Kyle M. Douglass, Marcel Stefko',
        'author_email': 'kyle.m.douglass@gmail.com',
        'url': 'https://github.com/kmdouglass/PolymerCpp',
        'download_url': 'https://github.com/kmdouglass/PolymerCpp/archive/v0.1.3.tar.gz',
        'keywords': ['polymer', 'wormlike chain', 'random walk'],
        'classifiers': ['Topic :: Scientific/Engineering :: Physics',
                        'Programming Language :: Python :: 3',
                        'Programming Language :: Python :: 3.5'],
        'license': 'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'packages': ['PolymerCpp'],
        'ext_modules': [module],
        'install_requires': ['numpy>=1.11.0', 'matplotlib>=2.0.0']}

setup(**config)
