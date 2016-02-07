# -*- coding: utf-8 -*-
#from distutils.core import setup
#from distutils import util
from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

if __name__=='__main__':
    setup(name='OceanLab',
             version='0.1',
             description='Python functions for Physical Oceanograhy',
             long_description=readme(),
             package_dir={'OceanLab':'OceanLab'},
             url='https://github.com/iuryt/ocean',
             author='Iury Tércio Simões de Sousa',
             author_email='simoesiury@gmail.com',
             license='MIT',
             py_modules=['ADCP','AO','CTD','EOF_fillgap','MOORING','SEAPLOT'],
             install_requires=['pandas','numpy','scipy','matplotlib',
                               'gsw','seawater','netCDF4'])
