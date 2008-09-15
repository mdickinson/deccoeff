from distutils.core import setup, Extension

files = ['_decimal.c', 'limb.c']

libraries = []
includes = []
setup(name = "_decimal", version = "0.1",
      ext_modules = [Extension("_decimal", files,
                               libraries = libraries,
                               include_dirs =  includes,
                               )
                     ],
      )
