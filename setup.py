#!/usr/bin/env python3
import os
import re
import sys
import sysconfig
import platform
import subprocess
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext


# Refer to "https://www.benjack.io/2017/06/12/python-cpp-tests.html"
# I changed minor parts of the code written in the above site.
class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError(
              "CMake must be installed to build the following extensions: " +
              ", ".join(e.name for e in self.extensions))

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(
                 os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable,
                      '-Wno-dev']
        build_args = ['--', '-j1']
        env = os.environ.copy()
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args,
                              cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args,
                              cwd=self.build_temp)

# download trng4 and pybind11 libraries from github and install them at here
module_name = 'install_libraries'
setup (
    name = module_name,
    version = '0.1',
    author = 'Dongkyu Kim',
    author_email = 'dkkim1005@gmail.com',
    description = 'Python binding for libraries',
    long_description = '',
    platforms = ['linux'],
    # add extension module
    ext_modules = [CMakeExtension('cpp_extension')],
    # add custom build_ext command
    cmdclass = dict(build_ext=CMakeBuild),
    zip_safe = False
)


# build main parts
module_name   = '_cpp_module'
include_dirs  = ['./pybind11/include', './trng4']
libraries     = ['openblas', 'lapack', 'fftw3', 'trng4']
library_dirs  = ['./trng4.build/trng']
define_macros = [('MAJOR_VERSION', '1'), ('MINOR_VERSION', '0'),
                 ('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION'), ('USE_TRNG4', None)]

module = Extension(module_name,
        define_macros = define_macros,
        include_dirs = include_dirs,
        libraries = libraries,
        library_dirs = library_dirs,
        sources = ['./src/etc.cpp', './src/hirschfye_qmc.cpp', './src/green_tools.cpp', './src/pymodule_hfqmc.cpp'],
        language = 'c++',
        extra_compile_args = ['-std=c++11', '-O3'])

setup (name = module_name,
       version = '1.0',
       description = '...',
       author = 'Dongkyu Kim',
       author_email = 'dkkim1005@gmail.com',
       url = '...',
       long_description = "HF-QMC",
       platforms = ['linux'],
       ext_modules = [module])

print ('\n@Note! Type the following command to include a search path for dynamic linker: export LD_LIBRARY_PATH=$(pwd)/trng4.build/trng:$LD_LIBRARY_PATH')
