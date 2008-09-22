from distutils.core import setup, Extension

files = ['deccoeff.c']

libraries = []
includes = []
setup(name = "deccoeff", version = "0.2",
      ext_modules = [Extension("deccoeff", files,
                               libraries = libraries,
                               include_dirs =  includes,
                               )
                     ],
      )
