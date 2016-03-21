from distutils.core import setup, Extension

""" Install guide:
  1. Copy src/Eigen directory to /usr/local/include
  2. run command in src directory: "sudo python setup.py install"
  3. You can now "import PolymerCpp"!
"""


module1 = Extension('PolymerCpp',
                    define_macros = [('MAJOR_VERSION', '1'),
                                     ('MINOR_VERSION', '0')],
                    include_dirs = ['/usr/local/include'],
                    sources = ['PyPolymerCpp.cpp','WLC.cpp',
                    'Path.cpp','RgDict.cpp','Misc.cpp',
                    'SAWLC.cpp','SAWLC_Rosenbluth.cpp',
                    'Stdafx.cpp','Stopwatch.cpp','PyUtils.cpp'],
                    extra_compile_args = ['-std=c++11','-O2'])

setup (name = 'PolymerCpp',
       version = '1.0',
       description = 'Python wrapper for PolymerCpp',
       author = 'Marcel Stefko',
       author_email = 'marcel.stefko@epfl.ch',
       url = 'https://github.com/MStefko/PolymerCpp',
       long_description = '''
Hacked together in a few hours.
''',
       ext_modules = [module1])