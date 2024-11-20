from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy

#ext = Extension("vlmMethod1", ["vlmMethod1.pyx"],
    #include_dirs = [numpy.get_include()])

#setup(ext_modules=[ext], cmdclass = {'build_ext': build_ext})

setup(
  name = 'vlmMethod1',
  ext_modules=[
    Extension('vlmMethod1',
              sources=['vlmMethod1.pyx'],
              extra_compile_args=['-O2'],
              language='c')
    ],
  cmdclass = {'build_ext': build_ext}
)