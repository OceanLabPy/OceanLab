# -*- coding: utf-8 -*-
from distutils.core import setup
from distutils import util

def readme():
    with open('README.md') as f:
        return f.read()

if __name__=='__main__':
    setup(name='ocean',
             version='0.2',
             description='Python functions for Physical Oceanograhy',
             long_description=readme(),
             url='https://github.com/iuryt/ocean',
             author='Iury Tércio Simões de Sousa',
             author_email='simoesiury@gmail.com',
             license='MIT',
             py_modules=['ADCP','AO','CTD','EOF_fillgap','MOORING','SEAPLOT'],
             install_requires=['pandas','numpy','scipy','matplotlib','basemap',
                               'gsw','seawater','netCDF4'],
             zip_safe=False)
