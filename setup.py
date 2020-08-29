#!/usr/bin/env python3
from distutils.core import setup, Extension
module_name = 'cpp_module'

include_dirs = ['./pybind11/include']
libraries    = ['stdc++', 'blas', 'lapack', 'trng4']

module = Extension(module_name,
        define_macros = [('MAJOR_VERSION', '1'), ('MINOR_VERSION', '0'),
                         ('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
        include_dirs = include_dirs,
        libraries = libraries,
        sources = ['hfqmc_pymodule.cpp'],
        language = 'c++',
        extra_compile_args = ['-std=c++11', '-O3'])

setup (name = module_name,
       version = '1.0',
       description = '...',
       author = 'Dongkyu Kim',
       author_email = 'dkkim1005@gmail.com',
       url = 'https://docs.python.org/extending/building',
       long_description = "HF-QMC",
       platforms = ['macos'],
       ext_modules = [module])
