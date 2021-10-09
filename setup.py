from distutils.core import setup, Extension
from Cython.Build import cythonize

ext = Extension("_flash", sources=["_flash.pyx"], extra_compile_args=["-O3"], language="c++")
setup(name="_flash", ext_modules=cythonize([ext]))
