from setuptools import setup

def read(file):
    return open(file, 'r').read()

LONG_DESCRIPTION = read('README.md')
LICENSE = read('LICENSE.txt')

setup(
    name='OceanLab',
    version='0.1.0',
    packages=['OceanLab'],
    include_package_data=True,
    description='Python functions for Physical Oceanography',
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    download_url = 'https://pypi.python.org/pypi/OceanLab',
    url='https://github.com/iuryt/OceanLab',
    author='Iury T. Simoes-Sousa',
    author_email='simoesiury@gmail.com',
    license=LICENSE,
    py_modules=['OA','EOF','DYN'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires = [
        'seawater ~= 3.3',
        'numpy ~= 1.18',
        'scipy ~= 1.6',
        'xarray ~= 0.18',
        'dask ~= 2021.06',
        'dask[distributed] ~= 2021.06'
    ],
)
